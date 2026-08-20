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
    bounds = [(-1.0, 1.0)] * ne + [(0.0, 1.0)] * ne
    F = None
    lower_bound = None

    for it in range(max_cuts):
        rows = list(ub_rows) + [r for r, _ in cuts]
        rhs = list(ub_rhs) + [b for _, b in cuts]
        res = linprog(obj, A_ub=np.array(rows), b_ub=np.array(rhs),
                      bounds=bounds, method="highs")
        if not res.success:
            return FloorResult(
                "infeasible", m, n_cuts=len(cuts),
                certificate={
                    "kind": "lp_infeasible_under_valid_cuts",
                    "n_loewner_cuts": len(cuts),
                    "n_structural_rows": n_structural,
                    "reason": ("every added cut is valid for the true feasible set, so "
                               "an infeasible relaxation proves the floor is unreachable "
                               "at width two on this graph"),
                },
                info={"iterations": it})
        x = np.asarray(res.x, dtype=float)
        lower_bound = float(res.fun)          # relaxation optimum is a valid lower bound
        F = np.eye(m)
        for e, k in edge_index.items():
            F[e[0], e[1]] = F[e[1], e[0]] = x[k]

        M = A.T @ F @ A - G_req
        w, V = np.linalg.eigh(M)
        if w[0] >= -tol:
            t = np.abs(np.array([x[k] for k in range(ne)]))
            return FloorResult(
                "feasible", m, F=F, t=t,
                cost=float(c_e @ t), lower_bound=lower_bound,
                n_cuts=len(cuts), slack_min_eig=float(w[0]),
                certificate={"kind": "primal_feasible_plus_relaxation_bound",
                             "min_eigenvalue_of_slack": float(w[0])},
                info={"iterations": it, "n_edges": ne})
        v = V[:, 0]
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
