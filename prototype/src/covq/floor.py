"""The information-floor compiler.

The compiler's contract is not "reproduce this QFIM". It is:

    min_Pi  Cost_C(Pi)      subject to      A^T F_Pi A  >=  G_req ,

where ``A`` (m x d) names the downstream parameter combinations the caller
actually cares about and ``G_req >= 0`` is the Fisher information floor those
combinations must clear.  Exact matching ``F_Pi = F_star`` survives as one mode
of this, not as the reason the compiler exists.

For pair-entangled labelled programs on a hardware graph ``H``, Theorem 1 makes
the feasible body exact, so the contract becomes a concrete convex program:

    min   sum_{e in E} c_e t_e
    s.t.  A^T F A >= G_req                     (Loewner, semidefinite)
          diag F = 1
          F_ij = 0                             for ij not in E
          -t_ij <= F_ij <= t_ij                for ij in E
          t in MATCH(H)                        (degree + blossom)

Two things make this useful rather than merely well-posed.

*The objective is a circuit count.*  ``MATCH(H)`` is down-closed, so at a
non-negative-cost optimum ``t_e = |F_e|``, and any matching decomposition
``|F| = sum_r p_r chi_{M_r}`` gives expected per-shot native cost
``sum_r p_r sum_{e in M_r} c_e = sum_e c_e |F_e|``.  The objective *is* the
expected Bell-pair bill.

*It is solved by cutting planes, so it certifies itself.*  The only
non-polyhedral constraint is the Loewner one, and a violated instance separates
by an eigenvector: if ``v`` is a minimal eigenvector of ``A^T F A - G_req`` then

    sum_e F_e * 2 (a_i . v)(a_j . v)  >=  v^T G_req v - v^T A^T A v

is a valid linear cut (``a_i`` is row ``i`` of ``A``).  Every LP relaxation
along the way is a valid *lower bound* on the true cost, and an infeasible
relaxation is a proof that the requested floor is unreachable at width two on
``H``.  The claim to make is convex, certificate-producing optimisation to a
prescribed precision -- not a strongly polynomial algorithm.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linprog

from .polytope import pairs
from .width import _matching_column_generation


@dataclass
class FloorResult:
    """One of Theorem 3's three outcomes, with the evidence for it."""

    status: str  # "feasible" | "infeasible" | "iteration_limit"
    m: int
    F: np.ndarray | None = None
    t: np.ndarray | None = None
    cost: float | None = None
    lower_bound: float | None = None
    n_cuts: int = 0
    slack_min_eig: float | None = None
    certificate: dict | None = None
    dual: "DualCertificate | None" = None
    info: dict = field(default_factory=dict)

    @property
    def gap(self) -> float | None:
        if self.cost is None or self.lower_bound is None:
            return None
        return float(self.cost - self.lower_bound)


def _odd_set_rows(m: int, edge_index: dict, max_enumerate: int = 16):
    """Blossom constraints, enumerated (exact) for the sizes in scope."""
    rows, rhs = [], []
    if m > max_enumerate:
        raise ValueError("odd-set enumeration capped; use Padberg-Rao separation")
    for size in range(3, m + 1, 2):
        for S in itertools.combinations(range(m), size):
            row = np.zeros(len(edge_index))
            hit = False
            for i, j in itertools.combinations(S, 2):
                k = edge_index.get((i, j))
                if k is not None:
                    row[k] = 1.0
                    hit = True
            if hit:
                rows.append(row)
                rhs.append((size - 1) / 2.0)
    return rows, rhs


def information_floor_compile(A: np.ndarray, G_req: np.ndarray, m: int, edges,
                              costs=None, tol: float = 1e-8,
                              max_cuts: int = 300) -> FloorResult:
    """Least-cost pair-entangled labelled program meeting a Fisher floor."""
    A = np.asarray(A, dtype=float).reshape(m, -1)
    G_req = np.atleast_2d(np.asarray(G_req, dtype=float))
    d = A.shape[1]
    if G_req.shape != (d, d):
        raise ValueError(f"G_req must be {d}x{d} to match A's {d} columns")

    E = sorted({tuple(sorted((int(a), int(b)))) for a, b in edges})
    if not E:
        raise ValueError("hardware graph has no edges")
    edge_index = {e: k for k, e in enumerate(E)}
    ne = len(E)
    c_e = np.ones(ne) if costs is None else np.array(
        [float(costs[e] if isinstance(costs, dict) else costs[k])
         for k, e in enumerate(E)], dtype=float)
    if np.any(c_e < 0):
        raise ValueError("negative edge costs break the t_e = |F_e| reduction")

    # variables x = [F_e (ne), t_e (ne)];  objective sums c_e t_e
    obj = np.concatenate([np.zeros(ne), c_e])
    ub_rows, ub_rhs = [], []
    for k in range(ne):                       # F_e - t_e <= 0 ; -F_e - t_e <= 0
        r = np.zeros(2 * ne); r[k] = 1.0; r[ne + k] = -1.0
        ub_rows.append(r); ub_rhs.append(0.0)
        r = np.zeros(2 * ne); r[k] = -1.0; r[ne + k] = -1.0
        ub_rows.append(r); ub_rhs.append(0.0)
    for i in range(m):                        # degree: sum_{e ni i} t_e <= 1
        r = np.zeros(2 * ne)
        for e, k in edge_index.items():
            if i in e:
                r[ne + k] = 1.0
        if r.any():
            ub_rows.append(r); ub_rhs.append(1.0)
    odd_rows, odd_rhs = _odd_set_rows(m, edge_index)
    for row, rhs in zip(odd_rows, odd_rhs):
        r = np.zeros(2 * ne); r[ne:] = row
        ub_rows.append(r); ub_rhs.append(rhs)
    n_structural = len(ub_rows)

    AtA = A.T @ A
    cuts: list[tuple[np.ndarray, float]] = []
    cut_vs: list[np.ndarray] = []
    bounds = [(-1.0, 1.0)] * ne + [(0.0, 1.0)] * ne
    F = None
    lower_bound = None

    for it in range(max_cuts):
        rows = list(ub_rows) + [r for r, _ in cuts]
        rhs = list(ub_rhs) + [b for _, b in cuts]
        res = linprog(obj, A_ub=np.array(rows), b_ub=np.array(rhs),
                      bounds=bounds, method="highs")
        if not res.success:
            cert = {"kind": "lp_infeasible_under_valid_cuts",
                    "n_loewner_cuts": len(cuts),
                    "n_structural_rows": n_structural,
                    "reason": ("every added cut is valid for the true feasible set, so "
                               "an infeasible relaxation proves the floor is unreachable "
                               "at width two on this graph")}
            dual = _dual_lp(A, G_req, m, E, cut_vs, c_e, homogeneous=True)
            if dual is not None and dual[3] > 1e-9:
                Y, mu, nu, val = dual
                dc = DualCertificate("farkas", Y, mu, nu, val, A, G_req, E, c_e)
                cert["kind"] = "farkas_dual_ray"
                cert["dual"] = {"objective": val, "verification": dc.verify()}
                return FloorResult("infeasible", m, n_cuts=len(cuts),
                                   certificate=cert, dual=dc, info={"iterations": it})
            return FloorResult("infeasible", m, n_cuts=len(cuts),
                               certificate=cert, info={"iterations": it})
        x = np.asarray(res.x, dtype=float)
        lower_bound = float(res.fun)          # relaxation optimum is a valid lower bound
        F = np.eye(m)
        for e, k in edge_index.items():
            F[e[0], e[1]] = F[e[1], e[0]] = x[k]

        M = A.T @ F @ A - G_req
        w, V = np.linalg.eigh(M)
        if w[0] >= -tol:
            t = np.abs(np.array([x[k] for k in range(ne)]))
            cert = {"kind": "primal_feasible_plus_relaxation_bound",
                    "min_eigenvalue_of_slack": float(w[0])}
            dc = None
            dual = _dual_lp(A, G_req, m, E, cut_vs, c_e, homogeneous=False)
            if dual is not None:
                Y, mu, nu, val = dual
                dc = DualCertificate("bound", Y, mu, nu, val, A, G_req, E, c_e)
                cert["dual_bound"] = val
                cert["dual_verification"] = dc.verify()
            return FloorResult(
                "feasible", m, F=F, t=t,
                cost=float(c_e @ t), lower_bound=lower_bound,
                n_cuts=len(cuts), slack_min_eig=float(w[0]),
                certificate=cert, dual=dc,
                info={"iterations": it, "n_edges": ne})
        v = V[:, 0]
        cut_vs.append(v.copy())
        coef = np.zeros(2 * ne)
        for e, k in edge_index.items():
            coef[k] = -2.0 * float(A[e[0]] @ v) * float(A[e[1]] @ v)
        cuts.append((coef, float(v @ AtA @ v - v @ G_req @ v)))

    return FloorResult("iteration_limit", m, F=F, lower_bound=lower_bound,
                       n_cuts=len(cuts), info={"max_cuts": max_cuts})


def emit_floor_program(result: FloorResult, edges=None):
    """Turn a feasible floor solution into a labelled signed-Bell-pair schedule."""
    from .programs import labelled_schedule_from_width

    if result.status != "feasible" or result.F is None:
        raise ValueError("no feasible F to emit")
    F = result.F
    m = result.m
    x = np.abs(F).copy()
    np.fill_diagonal(x, 0.0)
    mask = np.zeros((m, m))
    if edges is None:
        mask[:] = 1.0
        np.fill_diagonal(mask, 0.0)
    else:
        for a, b in edges:
            mask[a, b] = mask[b, a] = 1.0
    dec = _matching_column_generation(F, x, mask, tol=1e-9, max_rounds=500)
    prog = labelled_schedule_from_width(dec)

    # The cost identity: expected per-shot native edge usage equals sum_e |F_e|.
    expected_edges = sum(w * sum(1 for b in blocks if len(b) == 2)
                         for w, blocks, _ in dec.branches)
    return prog, dec, {
        "n_branches": dec.n_branches,
        "caratheodory_bound": int(np.count_nonzero(mask) // 2) + 1,
        "expected_active_edges_per_shot": float(expected_edges),
        "sum_abs_offdiag": float(x.sum() / 2.0),
        "cost_identity_residual": float(abs(expected_edges - x.sum() / 2.0)),
    }


def common_mode_contract(m: int, gamma: float):
    """The analytic unit test: ``A = 1/sqrt(m)``, ``G_req = [gamma]``.

    ``u^T F u = 1 + (2/m) sum_{i<j} F_ij``, so the contract needs
    ``sum_{i<j} F_ij >= m(gamma-1)/2`` and, at unit edge costs, the optimum is

        cost* = m (gamma - 1) / 2 ,

    attainable exactly for ``1 <= gamma <= 1 + 2 floor(m/2) / m`` -- the upper
    limit being the k=2 support-function value.  A product probe reaches only
    ``gamma = 1``.
    """
    u = np.ones((m, 1)) / np.sqrt(m)
    return u, np.array([[float(gamma)]]), max(0.0, m * (gamma - 1.0) / 2.0)


# ----------------------------------------------------------------------
# The explicit conic dual, and the matching separation oracle
# ----------------------------------------------------------------------
#
# Lagrangian dual of Theorem 2.  Write F = I + sum_e F_e (E_ij + E_ji) and
# S_e = a_i a_j^T + a_j a_i^T with a_i the i-th row of A, so that
# <Y, S_e> = 2 a_i^T Y a_j.  Dualising the Loewner constraint with Y >= 0, the
# degree constraints with mu >= 0 and the blossom constraints with nu >= 0:
#
#   max   <Y, G_req - A^T A>  -  sum_i mu_i  -  sum_S nu_S (|S|-1)/2
#   s.t.  |2 a_i^T Y a_j|  <=  c_ij + mu_i + mu_j + sum_{S superset ij} nu_S
#         Y >= 0,  mu >= 0,  nu >= 0
#
# Reading: Y prices the information requirement, mu prices qubit congestion,
# nu prices the odd-set (blossom) obstructions, and the constraint says the
# marginal information value of an edge cannot exceed its cost plus the
# congestion it consumes.
#
# Weak duality, in three lines.  For any primal-feasible (F, t) and any
# dual-feasible (Y, mu, nu),
#
#   <Y, G_req - A^T A>  <=  sum_e F_e * 2 a_i^T Y a_j            (Y, slack >= 0)
#                       <=  sum_e t_e (c_e + mu_i + mu_j + sum nu_S)
#                       <=  sum_e c_e t_e + sum_i mu_i + sum_S nu_S (|S|-1)/2 ,
#
# the last step using the degree and blossom inequalities on t.  With c = 0 the
# same chain proves *infeasibility* whenever the dual objective is positive:
# no primal-feasible point can exist.

@dataclass
class DualCertificate:
    """A checkable dual point: a cost lower bound, or a proof of infeasibility."""

    kind: str                      # "bound" | "farkas"
    Y: np.ndarray
    mu: np.ndarray
    nu: dict                       # odd set (tuple) -> price
    value: float
    A: np.ndarray
    G_req: np.ndarray
    edges: list
    costs: np.ndarray

    def verify(self, tol: float = 1e-7) -> dict:
        """Re-derive the certificate from scratch; never trust the solver."""
        Y = np.asarray(self.Y, dtype=float)
        sym = float(np.abs(Y - Y.T).max())
        psd = float(np.linalg.eigvalsh((Y + Y.T) / 2).min())
        c = np.zeros(len(self.edges)) if self.kind == "farkas" else self.costs
        worst = -np.inf
        for k, (i, j) in enumerate(self.edges):
            lhs = abs(2.0 * float(self.A[i] @ Y @ self.A[j]))
            rhs = c[k] + self.mu[i] + self.mu[j] + sum(
                p for S, p in self.nu.items() if i in S and j in S)
            worst = max(worst, lhs - rhs)
        obj = (float(np.sum(Y * (self.G_req - self.A.T @ self.A)))
               - float(self.mu.sum())
               - sum(p * (len(S) - 1) / 2.0 for S, p in self.nu.items()))
        ok = (sym <= tol and psd >= -tol and worst <= tol
              and self.mu.min(initial=0.0) >= -tol
              and (min(self.nu.values()) if self.nu else 0.0) >= -tol)
        if self.kind == "farkas":
            ok = ok and obj > tol
        return {"valid": bool(ok), "Y_symmetry_defect": sym, "Y_min_eigenvalue": psd,
                "worst_edge_violation": float(worst), "recomputed_objective": obj,
                "reported_value": self.value}


def matching_separation_oracle(t: np.ndarray, m: int, edges, tol: float = 1e-9):
    """Separation for ``MATCH(H)``: return a violated Edmonds inequality, or None.

    Returns ``("degree", (i,), lhs, 1.0)`` or ``("odd_set", S, lhs, (|S|-1)/2)``.

    This is the explicit oracle interface the complexity claim rests on.  The
    implementation enumerates odd sets, which is exact but exponential; Padberg-Rao
    gives the same answers in polynomial time and slots in here without changing a
    single caller.  Nothing else in the compiler knows how separation is done.
    """
    E = sorted({tuple(sorted((int(a), int(b)))) for a, b in edges})
    idx = {e: k for k, e in enumerate(E)}
    for i in range(m):
        lhs = float(sum(t[idx[e]] for e in E if i in e))
        if lhs > 1.0 + tol:
            return ("degree", (i,), lhs, 1.0)
    best = None
    for size in range(3, m + 1, 2):
        for S in itertools.combinations(range(m), size):
            lhs = float(sum(t[idx[e]] for e in E
                            if e[0] in S and e[1] in S))
            rhs = (size - 1) / 2.0
            if lhs > rhs + tol and (best is None or lhs - rhs > best[2] - best[3]):
                best = ("odd_set", S, lhs, rhs)
    return best


def _dual_lp(A, G_req, m, E, cut_vs, costs, homogeneous: bool):
    """Maximise the dual objective over (eta >= 0, mu >= 0, nu >= 0).

    ``Y`` is restricted to the cone spanned by the cut eigenvectors,
    ``Y = sum_k eta_k v_k v_k^T``, which is automatically PSD and is exactly the
    aggregation the cutting-plane LP already performs.
    """
    odd = [S for size in range(3, m + 1, 2) for S in itertools.combinations(range(m), size)]
    nk, nm, nn = len(cut_vs), m, len(odd)
    if nk == 0:
        return None
    AtA = A.T @ A
    # objective coefficients (we maximise, so negate for linprog)
    obj_eta = np.array([float(v @ (G_req - AtA) @ v) for v in cut_vs])
    obj = -np.concatenate([obj_eta, -np.ones(nm), -np.array([(len(S) - 1) / 2.0 for S in odd])])
    rows, rhs = [], []
    c = np.zeros(len(E)) if homogeneous else costs
    for k, (i, j) in enumerate(E):
        pe = np.array([2.0 * float(A[i] @ np.outer(v, v) @ A[j]) for v in cut_vs])
        base = np.zeros(nm + nn)
        base[i] -= 1.0
        base[j] -= 1.0
        for si, S in enumerate(odd):
            if i in S and j in S:
                base[nm + si] -= 1.0
        for sign in (1.0, -1.0):                      # |.| <= rhs
            rows.append(np.concatenate([sign * pe, base]))
            rhs.append(float(c[k]))
    if homogeneous:
        # the ray is a direction, so normalise it; in the non-homogeneous case weak
        # duality already bounds the objective by the primal optimum, and normalising
        # would only make the bound loose.
        rows.append(np.ones(nk + nm + nn))
        rhs.append(1.0)
    res = linprog(obj, A_ub=np.array(rows), b_ub=np.array(rhs),
                  bounds=(0.0, None), method="highs")
    if not res.success:
        return None
    x = np.asarray(res.x, dtype=float)
    eta, mu = x[:nk], x[nk:nk + nm]
    nu = {S: float(x[nk + nm + si]) for si, S in enumerate(odd)
          if x[nk + nm + si] > 1e-12}
    Y = sum(e * np.outer(v, v) for e, v in zip(eta, cut_vs))
    Y = np.zeros_like(AtA) if np.isscalar(Y) else Y
    return Y, mu, nu, -float(res.fun)


# ----------------------------------------------------------------------
# Shot-scaled information-floor compilation (manuscript Problem 8.1)
# ----------------------------------------------------------------------
#
# The normalised per-shot formulation above needs a downstream map A to escape
# the Loewner collapse.  The shot-scaled formulation does not: it allocates a
# non-negative number of uses t_{M,sigma} to each (matching, sign) branch and
# compares the *total* information against the floor directly.
#
#   min  sum_{M,sigma} c(M) t_{M,sigma}
#   s.t. sum_{M,sigma} t_{M,sigma} B_{M,sigma}  >=  G_star ,   t >= 0
#   with c(M) = c0 + sum_{e in M} c_e  and  B_{M,sigma} = I + sum_{ij in M} sigma_ij (E_ij + E_ji).
#
# Always feasible: the empty matching gives B = I, so t = lambda_max(G_star)
# shots of the product probe suffice, whence OPT <= c0 * lambda_max(G_star).
#
# The conic dual has one PSD variable and exponentially many branch constraints,
#
#   max <G_star, Y>   s.t.  Y >= 0,  <B_{M,sigma}, Y> <= c(M)  for all M, sigma,
#
# and separates in polynomial time, because signs are chosen edgewise and the
# remaining problem is a maximum-weight matching:
#
#   max_{M,sigma} ( <B_{M,sigma},Y> - c(M) )
#        = tr Y - c0 + max_{M} sum_{ij in M} ( 2|Y_ij| - c_ij ) .
#
# That identity is what makes the oracle-polynomial claim work, and it is
# checked against brute-force enumeration in the tests.

def branch_matrix(m: int, matching, signs=None) -> np.ndarray:
    """``B_{M,sigma} = I + sum_{ij in M} sigma_ij (E_ij + E_ji)``."""
    B = np.eye(m)
    for k, (i, j) in enumerate(matching):
        s = 1.0 if signs is None else float(signs[k])
        B[i, j] = B[j, i] = s
    return B


def price_branch(Y: np.ndarray, m: int, edges, c0: float, costs: dict):
    """Most-violated branch constraint, via Eq. (65). Returns (matching, signs, value)."""
    from .width import max_weight_matching

    W = np.zeros((m, m))
    for (i, j) in edges:
        w = 2.0 * abs(float(Y[i, j])) - float(costs[(i, j)])
        W[i, j] = W[j, i] = w
    M, best = max_weight_matching(W)
    M = [tuple(sorted(e)) for e in M]
    signs = [1.0 if Y[i, j] >= 0 else -1.0 for i, j in M]
    return M, signs, float(np.trace(Y)) - c0 + best


def shot_scaled_floor_compile(G_star: np.ndarray, m: int, edges, c0: float = 1.0,
                              costs=None, tol: float = 1e-8,
                              max_rounds: int = 200, max_cuts: int = 400) -> dict:
    """Solve Problem 8.1 by column generation over branches with a matching oracle.

    The restricted master is itself semidefinite, and is solved by the same
    Loewner cutting-plane scheme; its cut multipliers aggregate into the dual
    ``Y = sum_k eta_k v_k v_k^T``, which is PSD by construction and is exactly
    the dual variable the pricing step needs.
    """
    G_star = np.atleast_2d(np.asarray(G_star, dtype=float))
    E = sorted({tuple(sorted((int(a), int(b)))) for a, b in edges})
    ce = {e: 1.0 for e in E} if costs is None else {
        e: float(costs[e] if isinstance(costs, dict) else costs[k]) for k, e in enumerate(E)}

    cols = [(([], []), c0)]                      # the empty matching: B = I
    Y = np.zeros((m, m))
    history = []
    n_solved = 1                                 # columns present when t was last solved
    hit_limit = True
    for rnd in range(max_rounds):
        Bs = [branch_matrix(m, M, s) for (M, s), _ in cols]
        cvec = np.array([c for _, c in cols])
        cuts_v: list[np.ndarray] = []
        t = None
        for _ in range(max_cuts):               # master, by Loewner cutting planes
            if cuts_v:
                rows = np.array([[-float(v @ B @ v) for B in Bs] for v in cuts_v])
                rhs = np.array([-float(v @ G_star @ v) for v in cuts_v])
            else:
                rows, rhs = np.zeros((1, len(Bs))), np.array([0.0])
            res = linprog(cvec, A_ub=rows, b_ub=rhs, bounds=(0.0, None), method="highs")
            if not res.success:
                return {"status": "master_failed", "round": rnd}
            t = np.asarray(res.x, dtype=float)
            S = sum(ti * B for ti, B in zip(t, Bs)) - G_star
            w, V = np.linalg.eigh(S)
            if w[0] >= -tol:
                break
            cuts_v.append(V[:, 0].copy())
        n_solved = len(cols)
        eta = -np.asarray(res.ineqlin.marginals, dtype=float) if cuts_v else None
        Y = (sum(e * np.outer(v, v) for e, v in zip(eta, cuts_v))
             if cuts_v else np.zeros((m, m)))
        Y = np.asarray(Y, dtype=float)
        cost = float(cvec @ t)
        M, signs, reduced = price_branch(Y, m, E, c0, ce)
        history.append({"round": rnd, "cost": cost, "dual": float(np.sum(G_star * Y)),
                        "reduced_cost": reduced, "n_branches": len(cols)})
        if reduced <= tol:
            hit_limit = False
            break
        cols.append(((M, signs), c0 + sum(ce[e] for e in M)))

    cols = cols[:n_solved]                       # keep t and cols in step
    total = sum(ti * branch_matrix(m, M, s) for ti, ((M, s), _) in zip(t, cols))
    return {
        "status": "round_limit" if hit_limit else "solved",
        "cost": float(np.array([c for _, c in cols]) @ t),
        "dual_bound": float(np.sum(G_star * Y)),
        "Y": Y,
        "shots": t,
        "branches": [(M, s) for (M, s), _ in cols],
        "total_information": total,
        "floor_slack_min_eig": float(np.linalg.eigvalsh(total - G_star).min()),
        "product_probe_bound": float(c0 * np.linalg.eigvalsh(G_star).max()),
        "n_branches_used": int((t > 1e-10).sum()),
        "n_entangled_branches_used": int(sum(
            1 for ti, ((M, _), _) in zip(t, cols) if ti > 1e-10 and len(M) > 0)),
        "uses_entanglement": bool(any(
            ti > 1e-10 and len(M) > 0 for ti, ((M, _), _) in zip(t, cols))),
        "history": history,
    }
