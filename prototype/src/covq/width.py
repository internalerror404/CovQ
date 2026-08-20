"""The bounded branch-entanglement-width hierarchy Q^prog_{m,k}.

``Q^prog_{m,k}`` is the convex hull of block-diagonal sign-correlation matrices
whose blocks have size at most ``k``:

    Q^prog_{m,k} = conv{ M(pi, s) : max_B |B| <= k },
    M(pi, s)_ij = s_i s_j if i ~_pi j, else 0;  M(pi, s)_ii = 1.

Two levels are settled here.

k = 1
    ``Q^prog_{m,1} = {I}``: a width-1 branch is a product state, whose
    covariance is diagonal, and every branch has unit diagonal.

k = 2  (the result this module exists for)
    ``F in Q^prog_{m,2}``  iff  the matrix of off-diagonal *magnitudes*
    ``x_ij = |F_ij|`` lies in the matching polytope of ``K_m``, i.e. iff

        sum_{j != i} |F_ij| <= 1                       for every i,
        sum_{{i,j} subset S} |F_ij| <= (|S| - 1) / 2   for every odd |S| >= 3.

    Necessity: a width-2 branch has off-diagonal support on a matching with
    entries in [-1, 1], so |F| is dominated by a convex combination of matching
    incidence vectors, and the matching polytope is down-closed.
    Sufficiency is constructive: Edmonds' theorem writes ``x`` as a convex
    combination of matchings, and matching ``M`` is realised by a branch that
    puts a two-qubit signed cat on each edge of ``M`` and |+> on every unmatched
    qubit.

    Consequences that matter for the paper: membership at ``k = 2`` is
    polynomial-time decidable (Edmonds separation) while membership at ``k = m``
    is NP-hard; and every feasible width-2 target is realised by at most
    ``O(m^2)`` branches whose deepest primitive is a single CX.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import linprog

from .polytope import from_offdiag, offdiag, pairs


# ----------------------------------------------------------------------
# Matching machinery
# ----------------------------------------------------------------------

def max_weight_matching(W: np.ndarray) -> tuple[list[tuple[int, int]], float]:
    """Exact maximum-weight matching of ``K_m`` by subset DP (O(2^m m^2)).

    Negative weights are dropped (a matching may leave any vertex uncovered).
    """
    m = W.shape[0]
    if m > 22:
        raise ValueError("subset-DP matching capped at m=22")
    best = np.zeros(1 << m)
    choice: list[tuple[int, int] | None] = [None] * (1 << m)
    for S in range(1, 1 << m):
        i = (S & -S).bit_length() - 1
        rest = S & ~(1 << i)
        best[S] = best[rest]
        choice[S] = None
        for j in range(m):
            if j == i or not (rest >> j) & 1:
                continue
            w = W[i, j]
            if w <= 0:
                continue
            cand = w + best[rest & ~(1 << j)]
            if cand > best[S] + 1e-15:
                best[S] = cand
                choice[S] = (i, j)
    S = (1 << m) - 1
    matching: list[tuple[int, int]] = []
    while S:
        c = choice[S]
        i = (S & -S).bit_length() - 1
        if c is None:
            S &= ~(1 << i)
        else:
            matching.append(c)
            S &= ~((1 << c[0]) | (1 << c[1]))
    return sorted(matching), float(best[(1 << m) - 1])


@dataclass
class MatchingViolation:
    """A violated facet of the matching polytope: a width-2 infeasibility proof."""

    kind: str  # "degree" | "odd_set"
    subset: tuple[int, ...]
    lhs: float
    rhs: float

    @property
    def margin(self) -> float:
        return self.lhs - self.rhs

    def verify(self, F: np.ndarray) -> bool:
        x = np.abs(np.asarray(F, dtype=float))
        if self.kind == "degree":
            i = self.subset[0]
            lhs = float(x[i].sum() - x[i, i])
            return lhs > self.rhs + 1e-12
        S = list(self.subset)
        lhs = float(sum(x[i, j] for i, j in itertools.combinations(S, 2)))
        return lhs > self.rhs + 1e-12


def matching_polytope_violation(F: np.ndarray, tol: float = 1e-9,
                                max_enumerate: int = 22) -> MatchingViolation | None:
    """Most-violated Edmonds constraint for ``x = |offdiag(F)|``, if any.

    Odd sets are enumerated exhaustively (exact, exponential).  Polynomial-time
    separation exists (Padberg-Rao) and is what the complexity claim rests on;
    enumeration is used here because ``m <= 12`` in scope and an exhaustive
    search cannot silently miss a facet.
    """
    x = np.abs(np.asarray(F, dtype=float)).copy()
    m = x.shape[0]
    np.fill_diagonal(x, 0.0)
    worst: MatchingViolation | None = None

    for i in range(m):
        lhs = float(x[i].sum())
        if lhs > 1.0 + tol:
            v = MatchingViolation("degree", (i,), lhs, 1.0)
            if worst is None or v.margin > worst.margin:
                worst = v
    if m > max_enumerate:
        raise ValueError("odd-set enumeration capped; use Padberg-Rao separation")
    for size in range(3, m + 1, 2):
        for S in itertools.combinations(range(m), size):
            lhs = float(sum(x[i, j] for i, j in itertools.combinations(S, 2)))
            rhs = (size - 1) / 2.0
            if lhs > rhs + tol:
                v = MatchingViolation("odd_set", S, lhs, rhs)
                if worst is None or v.margin > worst.margin:
                    worst = v
    return worst


# ----------------------------------------------------------------------
# Width-2 membership and constructive decomposition
# ----------------------------------------------------------------------

@dataclass
class WidthDecomposition:
    """A labelled width-bounded schedule.

    ``branches[r] = (weight, blocks, signs)`` where ``blocks`` is a partition of
    ``range(m)`` with every block of size <= ``width``, and ``signs`` is a
    +-1 vector; the branch QFIM is ``M(blocks, signs)``.
    """

    m: int
    width: int
    branches: list[tuple[float, list[tuple[int, ...]], np.ndarray]]
    residual: float
    feasible: bool
    violation: MatchingViolation | None = None
    method: str = ""
    info: dict = field(default_factory=dict)

    @property
    def n_branches(self) -> int:
        return len(self.branches)

    def matrix(self) -> np.ndarray:
        F = np.zeros((self.m, self.m))
        for w, blocks, s in self.branches:
            B = np.zeros((self.m, self.m))
            for blk in blocks:
                for i in blk:
                    for j in blk:
                        B[i, j] = float(s[i]) * float(s[j])
            np.fill_diagonal(B, 1.0)
            F += w * B
        return F


def decompose_width2(F: np.ndarray, tol: float = 1e-9,
                     max_rounds: int = 500) -> WidthDecomposition:
    """Decide ``F in Q^prog_{m,2}`` and, if so, emit the schedule.

    Uses column generation on the matching polytope with an exact max-weight
    matching pricing oracle, so both outcomes are certified: a decomposition, or
    a violated Edmonds facet.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    viol = matching_polytope_violation(F, tol=tol)
    if viol is not None:
        return WidthDecomposition(m, 2, [], float("inf"), False, viol, "edmonds_separation")

    x = np.abs(F).copy()
    np.fill_diagonal(x, 0.0)
    sgn = np.sign(F)
    target = np.array([x[i, j] for i, j in pairs(m)], dtype=float)
    npair = target.size

    cols: list[np.ndarray] = [np.zeros(npair)]  # the empty matching
    matchings: list[list[tuple[int, int]]] = [[]]
    pair_pos = {p: k for k, p in enumerate(pairs(m))}

    rounds = 0
    lam = None
    while rounds < max_rounds:
        rounds += 1
        A = np.array(cols).T  # (npair, ncols)
        nv = A.shape[1]
        c = np.concatenate([np.zeros(nv), np.ones(2 * npair)])
        A_eq = np.zeros((npair + 1, nv + 2 * npair))
        A_eq[:npair, :nv] = A
        A_eq[:npair, nv:nv + npair] = np.eye(npair)
        A_eq[:npair, nv + npair:] = -np.eye(npair)
        A_eq[npair, :nv] = 1.0
        b_eq = np.concatenate([target, [1.0]])
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
        if not res.success:
            break
        lam = np.asarray(res.x[:nv], dtype=float)
        if float(res.fun) <= tol:
            break
        # A matching column has incidence entries in {0,1}, so its reduced cost
        # is -(sum_{e in M} y_e + y0) and it enters when that sum exceeds -y0.
        y = np.asarray(res.eqlin.marginals[:npair], dtype=float)
        y0 = float(res.eqlin.marginals[npair])
        W = np.zeros((m, m))
        for (i, j), val in zip(pairs(m), y):
            W[i, j] = W[j, i] = val
        M, best = max_weight_matching(W)
        if best + y0 <= 1e-9:
            break
        if any(M == prev for prev in matchings):
            break
        col = np.zeros(npair)
        for e in M:
            col[pair_pos[tuple(sorted(e))]] = 1.0
        cols.append(col)
        matchings.append(M)

    if lam is None:
        return WidthDecomposition(m, 2, [], float("inf"), False, None, "column_generation_failed")

    branches: list[tuple[float, list[tuple[int, ...]], np.ndarray]] = []
    for w, M in zip(lam, matchings):
        if w <= 1e-12:
            continue
        s = np.ones(m, dtype=np.int8)
        blocks: list[tuple[int, ...]] = []
        covered: set[int] = set()
        for i, j in M:
            blocks.append((i, j))
            covered.update((i, j))
            s[j] = np.int8(1 if sgn[i, j] >= 0 else -1)
        for i in range(m):
            if i not in covered:
                blocks.append((i,))
        branches.append((float(w), sorted(blocks), s))

    d = WidthDecomposition(m, 2, branches, 0.0, True, None, "edmonds_column_generation")
    d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
    d.feasible = d.residual <= 1e-8
    d.info = {"rounds": rounds, "columns": len(cols)}
    return d


# ----------------------------------------------------------------------
# General k (exact, small m)
# ----------------------------------------------------------------------

def partitions_with_max_block(m: int, k: int):
    """All set partitions of ``range(m)`` with every block of size <= ``k``."""

    def rec(i: int, blocks: list[list[int]]):
        if i == m:
            yield [tuple(b) for b in blocks]
            return
        for b in blocks:
            if len(b) < k:
                b.append(i)
                yield from rec(i + 1, blocks)
                b.pop()
        blocks.append([i])
        yield from rec(i + 1, blocks)
        blocks.pop()

    yield from rec(0, [])


def width_k_generators(m: int, k: int, cap: int = 400_000) -> tuple[np.ndarray, list]:
    """Extreme points of ``Q^prog_{m,k}`` as off-diagonal coordinate vectors."""
    pos = {p: idx for idx, p in enumerate(pairs(m))}
    rows: list[np.ndarray] = []
    meta: list[tuple[list[tuple[int, ...]], np.ndarray]] = []
    for part in partitions_with_max_block(m, k):
        sign_choices = []
        for blk in part:
            opts = []
            for mask in range(1 << (len(blk) - 1)):  # quotient the global block sign
                s = np.ones(len(blk), dtype=np.int8)
                for b in range(len(blk) - 1):
                    if (mask >> b) & 1:
                        s[b + 1] = -1
                opts.append(s)
            sign_choices.append(opts)
        for combo in itertools.product(*sign_choices):
            s = np.ones(m, dtype=np.int8)
            for blk, sv in zip(part, combo):
                for q, val in zip(blk, sv):
                    s[q] = val
            vec = np.zeros(len(pos))
            for blk in part:
                for i, j in itertools.combinations(sorted(blk), 2):
                    vec[pos[(i, j)]] = float(s[i]) * float(s[j])
            rows.append(vec)
            meta.append((list(part), s.copy()))
            if len(rows) > cap:
                raise ValueError(f"Q^prog_{{{m},{k}}} generator count exceeds cap {cap}")
    return np.array(rows).T, meta


def decompose_width_k(F: np.ndarray, k: int, max_m: int = 8, tol: float = 1e-9,
                      force_enumeration: bool = False) -> WidthDecomposition:
    """Exact LP membership in ``Q^prog_{m,k}`` by full generator enumeration.

    ``force_enumeration`` keeps the brute-force path even at ``k = 2``, which is
    what makes the matching-polytope characterisation *checkable*: two
    independent code paths must return the same verdict on every instance.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    if k == 2 and not force_enumeration:
        return decompose_width2(F, tol=tol)
    if m > max_m:
        raise ValueError(f"exact width-{k} membership capped at m={max_m}")
    A, meta = width_k_generators(m, k)
    nv = A.shape[1]
    npair = A.shape[0]
    c = np.concatenate([np.zeros(nv), np.ones(2 * npair)])
    A_eq = np.zeros((npair + 1, nv + 2 * npair))
    A_eq[:npair, :nv] = A
    A_eq[:npair, nv:nv + npair] = np.eye(npair)
    A_eq[:npair, nv + npair:] = -np.eye(npair)
    A_eq[npair, :nv] = 1.0
    b_eq = np.concatenate([offdiag(F), [1.0]])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
    if not res.success:
        return WidthDecomposition(m, k, [], float("inf"), False, None, "width_k_lp_failed")
    lam = np.asarray(res.x[:nv], dtype=float)
    branches = [(float(w), meta[i][0], meta[i][1]) for i, w in enumerate(lam) if w > 1e-12]
    d = WidthDecomposition(m, k, branches, 0.0, float(res.fun) <= tol, None, "width_k_exact_lp")
    d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
    d.feasible = d.residual <= 1e-8
    d.info = {"l1_distance": float(res.fun), "n_generators": nv}
    return d


def min_width(F: np.ndarray, max_m_exact: int = 8, tol: float = 1e-9) -> dict:
    """Smallest ``k`` with ``F in Q^prog_{m,k}`` (the compiler resource k_0)."""
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    off = np.abs(offdiag(F))
    if off.max(initial=0.0) <= tol:
        return {"k": 1, "method": "identity", "certified": True}
    d2 = decompose_width2(F, tol=tol)
    if d2.feasible:
        return {"k": 2, "method": d2.method, "certified": True,
                "n_branches": d2.n_branches, "violation": None}
    viol = None if d2.violation is None else {
        "kind": d2.violation.kind, "subset": list(d2.violation.subset),
        "lhs": d2.violation.lhs, "rhs": d2.violation.rhs}
    if m > max_m_exact:
        return {"k": None, "method": "width2_refuted_only", "certified": False,
                "violation": viol,
                "note": f"k >= 3; exact k>=3 membership capped at m={max_m_exact}"}
    for k in range(3, m + 1):
        dk = decompose_width_k(F, k, max_m=max_m_exact, tol=tol)
        if dk.feasible:
            return {"k": k, "method": dk.method, "certified": True,
                    "n_branches": dk.n_branches, "violation": viol}
    return {"k": None, "method": "infeasible_at_every_width", "certified": True,
            "violation": viol}


# ----------------------------------------------------------------------
# Hardware-native width 2:  Q^lab_{H,2}
# ----------------------------------------------------------------------

def decompose_width2_hardware(F: np.ndarray, edges, tol: float = 1e-9,
                              max_rounds: int = 500) -> WidthDecomposition:
    """Decide ``F in Q^lab_{H,2}`` for a hardware graph ``H = (V, E)``.

        Q^lab_{H,2} = { F : diag F = 1,
                        F_ij = 0 for ij not in E,
                        (|F_ij|)_{ij in E} in MATCH(H) }

    This is the form the compiler actually needs: it maps a target QFIM directly
    onto a schedule of *parallel Bell pairs on physical edges*, with no routing
    inside the primitive backend.  Routing is a separate backend that enlarges
    the admissible edge set at an explicit SWAP cost, and must be priced as
    such rather than assumed free.
    """
    F = np.asarray(F, dtype=float)
    m = F.shape[0]
    E = {tuple(sorted((int(a), int(b)))) for a, b in edges}

    offgraph = [(i, j) for i, j in pairs(m)
                if (i, j) not in E and abs(F[i, j]) > tol]
    if offgraph:
        i, j = max(offgraph, key=lambda p: abs(F[p[0], p[1]]))
        v = MatchingViolation("off_graph_support", (i, j), abs(float(F[i, j])), 0.0)
        return WidthDecomposition(m, 2, [], float("inf"), False, v, "hardware_support")

    x = np.abs(F).copy()
    np.fill_diagonal(x, 0.0)
    for i, j in pairs(m):
        if (i, j) not in E:
            x[i, j] = x[j, i] = 0.0

    # Edmonds separation restricted to H
    worst: MatchingViolation | None = None
    for i in range(m):
        lhs = float(x[i].sum())
        if lhs > 1.0 + tol:
            v = MatchingViolation("degree", (i,), lhs, 1.0)
            if worst is None or v.margin > worst.margin:
                worst = v
    for size in range(3, m + 1, 2):
        for S in itertools.combinations(range(m), size):
            lhs = float(sum(x[i, j] for i, j in itertools.combinations(S, 2)
                            if (i, j) in E))
            rhs = (size - 1) / 2.0
            if lhs > rhs + tol:
                v = MatchingViolation("odd_set", S, lhs, rhs)
                if worst is None or v.margin > worst.margin:
                    worst = v
    if worst is not None:
        return WidthDecomposition(m, 2, [], float("inf"), False, worst,
                                  "edmonds_separation_hardware")

    mask = np.zeros((m, m))
    for i, j in E:
        if i < m and j < m:
            mask[i, j] = mask[j, i] = 1.0
    d = _matching_column_generation(F, x, mask, tol=tol, max_rounds=max_rounds)
    d.method = "edmonds_column_generation_hardware"
    d.info["n_edges"] = len(E)
    return d


def _matching_column_generation(F: np.ndarray, x: np.ndarray, mask: np.ndarray,
                                tol: float, max_rounds: int) -> WidthDecomposition:
    """Shared column-generation core; ``mask`` restricts the admissible edges."""
    m = F.shape[0]
    sgn = np.sign(F)
    target = np.array([x[i, j] for i, j in pairs(m)], dtype=float)
    npair = target.size
    pair_pos = {p: k for k, p in enumerate(pairs(m))}

    cols: list[np.ndarray] = [np.zeros(npair)]
    matchings: list[list[tuple[int, int]]] = [[]]
    lam = None
    rounds = 0
    while rounds < max_rounds:
        rounds += 1
        A = np.array(cols).T
        nv = A.shape[1]
        c = np.concatenate([np.zeros(nv), np.ones(2 * npair)])
        A_eq = np.zeros((npair + 1, nv + 2 * npair))
        A_eq[:npair, :nv] = A
        A_eq[:npair, nv:nv + npair] = np.eye(npair)
        A_eq[:npair, nv + npair:] = -np.eye(npair)
        A_eq[npair, :nv] = 1.0
        b_eq = np.concatenate([target, [1.0]])
        res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
        if not res.success:
            break
        lam = np.asarray(res.x[:nv], dtype=float)
        if float(res.fun) <= tol:
            break
        y = np.asarray(res.eqlin.marginals[:npair], dtype=float)
        y0 = float(res.eqlin.marginals[npair])
        W = np.zeros((m, m))
        for (i, j), val in zip(pairs(m), y):
            W[i, j] = W[j, i] = val
        W = W * mask
        M, best = max_weight_matching(W)
        if best + y0 <= 1e-9 or any(M == prev for prev in matchings):
            break
        col = np.zeros(npair)
        for e in M:
            col[pair_pos[tuple(sorted(e))]] = 1.0
        cols.append(col)
        matchings.append(M)

    if lam is None:
        return WidthDecomposition(m, 2, [], float("inf"), False, None, "cg_failed")
    branches = []
    for w, M in zip(lam, matchings):
        if w <= 1e-12:
            continue
        s = np.ones(m, dtype=np.int8)
        blocks: list[tuple[int, ...]] = []
        covered: set[int] = set()
        for i, j in M:
            blocks.append((i, j))
            covered.update((i, j))
            s[j] = np.int8(1 if sgn[i, j] >= 0 else -1)
        for i in range(m):
            if i not in covered:
                blocks.append((i,))
        branches.append((float(w), sorted(blocks), s))
    d = WidthDecomposition(m, 2, branches, 0.0, True, None, "cg")
    d.residual = float(np.linalg.norm(d.matrix() - F, "fro"))
    d.feasible = d.residual <= 1e-8
    d.info = {"rounds": rounds, "columns": len(cols),
              "parallel_bell_depth": max((sum(1 for b in bl if len(b) == 2) > 0)
                                         for _, bl, _ in branches) if branches else 0}
    return d
