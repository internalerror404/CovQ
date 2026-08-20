"""The sign-correlation polytope Q_m = conv{ s s^T : s in {+-1}^m }.

Everything here is about the *Z frame* target: a unit-diagonal, zero-mean
target QFIM for generators ``Z_1..Z_m``.

Two facts fix the scope of this module, and both are inherited rather than new:

* ``F`` with unit diagonal is realisable by a zero-mean pure state iff
  ``F in Q_m`` (Pitowsky's correlation polytope, transported to the +-1 sign
  convention);
* deciding membership is NP-hard in general, so every exact routine here is
  explicitly exponential and explicitly capped by ``m``.

The one piece with a genuinely quantum reading is the certificate family:
for integer ``b`` with ``sum_i b_i`` odd, every realisable ``F`` obeys

    b^T F b = Var(sum_i b_i P_i) >= 1,

because ``sum_i b_i P_i`` has integer spectrum of the parity of ``sum_i b_i``.
A violated inequality is therefore an *uncertainty-relation* obstruction, not
merely a separating hyperplane.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linprog

DEFAULT_ATOL = 1e-9


# ----------------------------------------------------------------------
# Basic geometry
# ----------------------------------------------------------------------

def pairs(m: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(m) for j in range(i + 1, m)]


def vertices(m: int) -> np.ndarray:
    """All sign vectors with ``s_0 = +1`` (quotienting ``s ~ -s``)."""
    if m > 24:
        raise ValueError(f"refusing to enumerate 2^{m - 1} vertices")
    k = m - 1
    bits = ((np.arange(1 << k)[:, None] >> np.arange(k)) & 1).astype(np.int8)
    s = np.ones((1 << k, m), dtype=np.int8)
    s[:, 1:] = 1 - 2 * bits
    return s


def vertex_pair_matrix(m: int, s: np.ndarray | None = None) -> np.ndarray:
    """``A[p, v] = s_v[i] * s_v[j]`` for pair ``p = (i, j)``."""
    s = vertices(m) if s is None else s
    idx = pairs(m)
    ii = np.array([i for i, _ in idx])
    jj = np.array([j for _, j in idx])
    return (s[:, ii] * s[:, jj]).T.astype(float)


def offdiag(F: np.ndarray) -> np.ndarray:
    m = F.shape[0]
    idx = pairs(m)
    return np.array([F[i, j] for i, j in idx], dtype=float)


def from_offdiag(vec: np.ndarray, m: int) -> np.ndarray:
    F = np.eye(m)
    for (i, j), v in zip(pairs(m), vec):
        F[i, j] = F[j, i] = v
    return F


def is_unit_diagonal_psd(F: np.ndarray, atol: float = 1e-8) -> tuple[bool, str]:
    F = np.asarray(F, dtype=float)
    if F.ndim != 2 or F.shape[0] != F.shape[1]:
        return False, "not square"
    if not np.allclose(F, F.T, atol=atol):
        return False, "not symmetric"
    if not np.allclose(np.diag(F), 1.0, atol=atol):
        return False, "diagonal is not unit"
    w = np.linalg.eigvalsh(F)
    if w.min() < -atol:
        return False, f"not PSD (min eig {w.min():.3e})"
    return True, "ok"


# ----------------------------------------------------------------------
# Certificates of infeasibility
# ----------------------------------------------------------------------

@dataclass
class Certificate:
    """A checkable proof that ``F`` is outside ``Q_m``."""

    kind: str  # "hypermetric" | "lp_hyperplane"
    data: dict
    margin: float

    def verify(self, F: np.ndarray, s: np.ndarray | None = None) -> bool:
        F = np.asarray(F, dtype=float)
        m = F.shape[0]
        if self.kind == "hypermetric":
            b = np.asarray(self.data["b"], dtype=float)
            if int(round(b.sum())) % 2 == 0:
                return False
            return float(b @ F @ b) < 1.0 - 1e-12
        if self.kind == "lp_hyperplane":
            y = np.asarray(self.data["y"], dtype=float)
            t = float(self.data["t"])
            s = vertices(m) if s is None else s
            A = vertex_pair_matrix(m, s)
            if float(np.max(y @ A)) > t + 1e-9:
                return False
            return float(y @ offdiag(F)) > t + 1e-12
        return False


def hypermetric_certificate(F: np.ndarray, max_entry: int = 1,
                            rng: np.random.Generator | None = None,
                            local_search_restarts: int = 0) -> Certificate | None:
    """Search for integer ``b`` with odd sum and ``b^T F b < 1``.

    Exhaustive over ``b in {-max_entry..max_entry}^m`` when that is small
    enough, plus optional randomised local search for larger ``m``.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    best: tuple[float, np.ndarray] | None = None

    span = list(range(-max_entry, max_entry + 1))
    if (2 * max_entry + 1) ** m <= 2_000_000:
        for b in itertools.product(span, repeat=m):
            if sum(b) % 2 == 0:
                continue
            bv = np.array(b, dtype=float)
            val = float(bv @ F @ bv)
            if best is None or val < best[0]:
                best = (val, bv)
    if local_search_restarts and rng is not None:
        for _ in range(local_search_restarts):
            bv = rng.integers(-max_entry, max_entry + 1, size=m).astype(float)
            if int(bv.sum()) % 2 == 0:
                bv[int(rng.integers(0, m))] += 1
            improved = True
            while improved:
                improved = False
                cur = float(bv @ F @ bv)
                for i in range(m):
                    for delta in (-2, 2):  # preserve parity of sum(b)
                        cand = bv.copy()
                        cand[i] += delta
                        if np.abs(cand).max() > max_entry + 2:
                            continue
                        val = float(cand @ F @ cand)
                        if val < cur - 1e-12:
                            bv, cur, improved = cand, val, True
            val = float(bv @ F @ bv)
            if best is None or val < best[0]:
                best = (val, bv)

    if best is None or best[0] >= 1.0 - 1e-12:
        return None
    return Certificate("hypermetric", {"b": best[1].astype(int).tolist()}, margin=1.0 - best[0])


def lp_separating_hyperplane(F: np.ndarray, s: np.ndarray | None = None) -> Certificate | None:
    """Exact LP separation against the full (quotiented) vertex list.

    Solves ``max <y, F> - t`` s.t. ``<y, V> <= t`` for every vertex and
    ``|y| <= 1``.  A positive optimum is a complete proof of infeasibility.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    s = vertices(m) if s is None else s
    A = vertex_pair_matrix(m, s)  # (n_pairs, n_vertices)
    npair, nvert = A.shape
    # variables: y (npair), t (1);  maximise y.f - t  ->  minimise -y.f + t
    c = np.concatenate([-offdiag(F), [1.0]])
    A_ub = np.hstack([A.T, -np.ones((nvert, 1))])
    b_ub = np.zeros(nvert)
    bounds = [(-1.0, 1.0)] * npair + [(None, None)]
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    if not res.success:
        return None
    value = -float(res.fun)
    if value <= 1e-9:
        return None
    y = np.asarray(res.x[:npair], dtype=float)
    t = float(res.x[npair])
    return Certificate("lp_hyperplane", {"y": y.tolist(), "t": t}, margin=value)


# ----------------------------------------------------------------------
# Exact feasibility and decomposition
# ----------------------------------------------------------------------

@dataclass
class Decomposition:
    """``F = sum_r weights[r] * outer(signs[r], signs[r])`` within ``residual``."""

    signs: np.ndarray  # (q, m) int8
    weights: np.ndarray  # (q,)
    residual: float
    feasible: bool
    certificate: Certificate | None = None
    method: str = ""
    info: dict = field(default_factory=dict)

    @property
    def support(self) -> int:
        return int(self.signs.shape[0])

    def matrix(self) -> np.ndarray:
        m = self.signs.shape[1]
        if self.support == 0:
            return np.zeros((m, m))
        S = self.signs.astype(float)
        return np.einsum("r,ri,rj->ij", self.weights, S, S)


def exact_decompose(F: np.ndarray, atol: float = DEFAULT_ATOL,
                    certify: bool = True) -> Decomposition:
    """Exact LP feasibility and decomposition over the full vertex list.

    HiGHS returns a basic solution, whose support is automatically bounded by
    the number of equality rows, ``C(m,2) + 1`` -- i.e. Caratheodory sparsity
    comes for free.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    s = vertices(m)
    A = vertex_pair_matrix(m, s)
    A_eq = np.vstack([A, np.ones((1, A.shape[1]))])
    b_eq = np.concatenate([offdiag(F), [1.0]])
    res = linprog(np.zeros(A.shape[1]), A_eq=A_eq, b_eq=b_eq,
                  bounds=(0.0, None), method="highs")
    if res.success:
        lam = np.asarray(res.x, dtype=float)
        keep = lam > 1e-12
        d = Decomposition(s[keep], lam[keep], 0.0, True, None, "exact_lp")
        d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
        d.feasible = d.residual <= atol
        d.info = {"n_vertices_enumerated": int(s.shape[0])}
        return d
    cert = None
    if certify:
        cert = hypermetric_certificate(F, max_entry=1) or lp_separating_hyperplane(F, s)
    d = Decomposition(np.zeros((0, m), dtype=np.int8), np.zeros(0), float("inf"),
                      False, cert, "exact_lp")
    d.info = {"n_vertices_enumerated": int(s.shape[0])}
    return d


def caratheodory_reduce(dec: Decomposition, atol: float = 1e-10) -> Decomposition:
    """Reduce support to at most ``C(m,2) + 1`` atoms without changing ``F``."""
    from scipy.linalg import null_space

    S = dec.signs.astype(float).copy()
    lam = dec.weights.astype(float).copy()
    m = S.shape[1]
    limit = m * (m - 1) // 2 + 1
    while S.shape[0] > limit:
        A = np.vstack([vertex_pair_matrix(m, S.astype(np.int8)), np.ones((1, S.shape[0]))])
        ns = null_space(A)
        if ns.size == 0:
            break
        w = ns[:, 0]
        pos = w > atol
        neg = w < -atol
        if not pos.any() and not neg.any():
            break
        ratios = np.where(pos, lam / np.where(pos, w, 1.0), np.inf)
        theta_pos = ratios.min() if pos.any() else np.inf
        ratios_n = np.where(neg, lam / np.where(neg, -w, 1.0), np.inf)
        theta_neg = ratios_n.min() if neg.any() else np.inf
        if theta_pos <= theta_neg:
            lam = lam - theta_pos * w
        else:
            lam = lam + theta_neg * w
        keep = lam > atol
        S, lam = S[keep], lam[keep]
    out = Decomposition(S.astype(np.int8), lam, 0.0, dec.feasible, dec.certificate,
                        dec.method + "+caratheodory", dict(dec.info))
    out.residual = float(np.linalg.norm(out.matrix() - dec.matrix(), "fro"))
    return out


def min_support_decompose(F: np.ndarray, max_m: int = 8,
                          time_limit: float = 60.0) -> Decomposition:
    """Minimum-support decomposition by MILP (exact, small ``m`` only)."""
    from scipy.optimize import Bounds, LinearConstraint, milp

    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    if m > max_m:
        raise ValueError(f"min_support_decompose capped at m={max_m}")
    s = vertices(m)
    A = vertex_pair_matrix(m, s)
    nv = A.shape[1]
    # variables [lambda (nv), y (nv binary)]
    c = np.concatenate([np.zeros(nv), np.ones(nv)])
    eq = np.hstack([np.vstack([A, np.ones((1, nv))]), np.zeros((A.shape[0] + 1, nv))])
    b = np.concatenate([offdiag(F), [1.0]])
    link = np.hstack([np.eye(nv), -np.eye(nv)])  # lambda - y <= 0
    cons = [
        LinearConstraint(eq, b, b),
        LinearConstraint(link, -np.inf, 0.0),
    ]
    integrality = np.concatenate([np.zeros(nv), np.ones(nv)])
    res = milp(c, constraints=cons, integrality=integrality,
               bounds=Bounds(lb=np.zeros(2 * nv), ub=np.ones(2 * nv)),
               options={"time_limit": time_limit})
    if not res.success or res.x is None:
        return Decomposition(np.zeros((0, m), dtype=np.int8), np.zeros(0),
                             float("inf"), False, None, "min_support_milp")
    lam = np.asarray(res.x[:nv], dtype=float)
    keep = lam > 1e-9
    d = Decomposition(s[keep], lam[keep], 0.0, True, None, "min_support_milp")
    d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
    d.feasible = d.residual <= 1e-8
    d.info = {"milp_status": str(res.status), "objective": float(res.fun)}
    return d


def column_generation(F: np.ndarray, max_rounds: int = 200,
                      init: int = 8, rng: np.random.Generator | None = None,
                      atol: float = DEFAULT_ATOL) -> Decomposition:
    """LP over a growing vertex subset; pricing is exact but exponential.

    The pricing subproblem is ``max_s s^T Y s``, an unconstrained binary
    quadratic (MaxCut-type) problem.  Its hardness is precisely why the exact
    arm cannot be pushed past small ``m``: this is a re-derivation of the known
    hardness, not a new obstruction.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    rng = np.random.default_rng(0) if rng is None else rng
    npair = m * (m - 1) // 2

    cols = [np.ones(m, dtype=np.int8)]
    for _ in range(init - 1):
        s = rng.choice(np.array([-1, 1], dtype=np.int8), size=m)
        s = s * s[0]
        cols.append(s.astype(np.int8))
    rounds = 0
    while rounds < max_rounds:
        rounds += 1
        S = np.array(cols, dtype=np.int8)
        A = vertex_pair_matrix(m, S)
        nv = A.shape[1]
        # min ||residual||_1 : lambda>=0, sum lambda = 1
        c = np.concatenate([np.zeros(nv), np.ones(2 * npair)])
        A_eq = np.zeros((npair + 1, nv + 2 * npair))
        A_eq[:npair, :nv] = A
        A_eq[:npair, nv:nv + npair] = np.eye(npair)
        A_eq[:npair, nv + npair:] = -np.eye(npair)
        A_eq[npair, :nv] = 1.0
        b_eq = np.concatenate([offdiag(F), [1.0]])
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
        if not res.success:
            break
        lam = np.asarray(res.x[:nv], dtype=float)
        obj = float(res.fun)
        if obj <= atol:
            keep = lam > 1e-12
            d = Decomposition(S[keep], lam[keep], 0.0, True, None, "column_generation")
            d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
            d.feasible = d.residual <= 1e-8
            d.info = {"rounds": rounds, "columns": int(S.shape[0])}
            return d
        # Reduced cost of a column v is -(sum_p y_p A[p,v] + y0); it enters when
        # that is negative, i.e. when sum_{i<j} y_ij s_i s_j + y0 > 0.  Note
        # s^T Y s = 2 sum_{i<j} y_ij s_i s_j for Y with zero diagonal.
        y = np.asarray(res.eqlin.marginals[:npair], dtype=float)
        y0 = float(res.eqlin.marginals[npair])
        Y = from_offdiag(y, m) - np.eye(m)
        best_s, best_val = _max_quadratic(Y, m)
        if 0.5 * best_val + y0 <= 1e-9:
            break
        if any(np.array_equal(best_s, c) for c in cols):
            break  # pricing returned a column already present: degenerate, stop
        cols.append(best_s.astype(np.int8))
    S = np.array(cols, dtype=np.int8)
    d = exact_decompose_over(F, S)
    d.method = "column_generation(exhausted)"
    d.info = {"rounds": rounds, "columns": int(S.shape[0])}
    return d


def _max_quadratic(Y: np.ndarray, m: int) -> tuple[np.ndarray, float]:
    """Exact ``max_s s^T Y s`` by enumeration (capped at m <= 20)."""
    if m > 20:
        raise ValueError("exact pricing capped at m=20")
    s = vertices(m).astype(float)
    vals = np.einsum("vi,ij,vj->v", s, Y, s)
    k = int(np.argmax(vals))
    return s[k].astype(np.int8), float(vals[k])


def exact_decompose_over(F: np.ndarray, S: np.ndarray,
                         atol: float = DEFAULT_ATOL) -> Decomposition:
    """Best l1 fit of ``F`` by a convex combination of the given sign vectors."""
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    npair = m * (m - 1) // 2
    A = vertex_pair_matrix(m, S)
    nv = A.shape[1]
    c = np.concatenate([np.zeros(nv), np.ones(2 * npair)])
    A_eq = np.zeros((npair + 1, nv + 2 * npair))
    A_eq[:npair, :nv] = A
    A_eq[:npair, nv:nv + npair] = np.eye(npair)
    A_eq[:npair, nv + npair:] = -np.eye(npair)
    A_eq[npair, :nv] = 1.0
    b_eq = np.concatenate([offdiag(F), [1.0]])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
    if not res.success:
        return Decomposition(np.zeros((0, m), dtype=np.int8), np.zeros(0),
                             float("inf"), False, None, "restricted_lp")
    lam = np.asarray(res.x[:nv], dtype=float)
    keep = lam > 1e-12
    d = Decomposition(S[keep], lam[keep], 0.0, False, None, "restricted_lp")
    d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
    d.feasible = d.residual <= atol
    return d


def frank_wolfe(F: np.ndarray, iters: int = 200, away_steps: bool = True,
                tol: float = 1e-12) -> Decomposition:
    """Fedorov-Wynn / Frank-Wolfe sparse approximation with a duality gap.

    This is the classical exchange algorithm of optimal experimental design,
    applied to the design space of signed cat settings; it is reported as such
    and is not claimed as new.  Away steps are on by default -- without them the
    O(1/t) rate makes the arm look far weaker than it is, which would be an
    unfair comparison in the compiler's own favour.

    The returned ``info['final_duality_gap']`` is a genuine certificate: it upper
    bounds the suboptimality of the reported approximation.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    atoms = [np.ones(m, dtype=np.int8)]
    weights = np.array([1.0])
    gap_hist: list[float] = []
    n_away = 0
    for _ in range(iters):
        S = np.array(atoms, dtype=float)
        X = np.einsum("r,ri,rj->ij", weights, S, S)
        G = X - F
        s_fw, _ = _max_quadratic(-G, m)
        V_fw = np.outer(s_fw.astype(float), s_fw.astype(float))
        gap_fw = float(np.sum(G * (X - V_fw)))
        gap_hist.append(gap_fw)
        if gap_fw <= tol:
            break

        use_away = False
        if away_steps and len(atoms) > 1:
            vals = [float(np.sum(G * np.outer(a.astype(float), a.astype(float))))
                    for a in atoms]
            a_idx = int(np.argmax(vals))
            V_aw = np.outer(atoms[a_idx].astype(float), atoms[a_idx].astype(float))
            gap_aw = float(np.sum(G * (V_aw - X)))
            if gap_aw > gap_fw:
                use_away = True

        if use_away:
            D = X - V_aw
            wa = float(weights[a_idx])
            gmax = wa / (1.0 - wa) if wa < 1.0 - 1e-15 else 0.0
        else:
            D = V_fw - X
            gmax = 1.0
        dd = float(np.sum(D * D))
        if dd <= 1e-30 or gmax <= 0.0:
            break
        gamma = float(np.clip(np.sum((F - X) * D) / dd, 0.0, gmax))
        if gamma <= 0.0:
            break

        if use_away:
            n_away += 1
            weights = weights * (1.0 + gamma)
            weights[a_idx] -= gamma
        else:
            match = [i for i, a in enumerate(atoms) if np.array_equal(a, s_fw)]
            weights = weights * (1.0 - gamma)
            if match:
                weights[match[0]] += gamma
            else:
                weights = np.append(weights, gamma)
                atoms.append(s_fw)
        keep = weights > 1e-14
        atoms = [a for a, k in zip(atoms, keep) if k]
        weights = weights[keep]
        weights = weights / weights.sum()

    S = np.array(atoms, dtype=np.int8)
    d = Decomposition(S, weights, 0.0, False, None,
                      "frank_wolfe_away" if away_steps else "frank_wolfe")
    d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
    d.feasible = d.residual <= 1e-8
    d.info = {"iterations": len(gap_hist),
              "away_steps": n_away,
              "final_duality_gap": gap_hist[-1] if gap_hist else 0.0,
              "support": int(S.shape[0]),
              "operator_norm_error": float(np.linalg.norm(d.matrix() - F, 2))}
    return d
