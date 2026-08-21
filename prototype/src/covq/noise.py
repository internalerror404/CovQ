"""Block-local noise, and whether it preserves matching-based compilation.

Equation (111) in its reformulated shape is a *noisy operational information
floor*, not a weighted norm-to-target heuristic:

    min  C_total(Pi, n, R)   s.t.   A^T F_C^N(Pi, n, R) A >= G_req,

with `F_C^N = sum_b n_b F_{C,b}^{N,R_b}` the classical Fisher matrix of the
**declared readout** under the noise model, `n_b` the exposure of branch `b`,
and `R_b` its pilot and production readout.  The mixed-state QFIM enters only
as a reference upper bound `F_C^N <= F_Q^N`; the generator covariance is a
diagnostic and nothing more -- the pinching proposition shows it can be exactly
constant while the program is statistically dead.

The central structural question this module exists to answer:

    **Does block-local noise preserve matching-based pricing?**

If the channel and the readout are block-local, every branch template splits as

    F_{C,M,sigma}^N  =  F_0^N  +  sum_{e in M} DeltaF_{e,sigma_e}^N,

because the state, the noise and the readout all factorise over blocks.  Then
for a dual matrix `Q = A Y A^T` the pricing subproblem collapses to

    max_{M,sigma} <Q, F^N_{C,M,sigma}> - c(M)
      = <Q, F_0^N> - c_0
        + max_{M in MATCH(H)} sum_{e in M} [ max_sigma <Q, DeltaF_{e,sigma}^N> - c_e ],

which is still a maximum-weight matching -- the *same* oracle as the ideal
Eq (65), with reweighted edges.  The noise changes the templates and the edge
weights; it does not change the combinatorics.  Everything here is computed on
one- and two-qubit subsystems rather than on the full `2^m` state, which is
exactly the computational payoff of block-locality.

Scope, declared up front.  This holds for noise that acts *within* blocks.  It
is not claimed for crosstalk, correlated branch noise, route collisions, or
coherent errors that couple distinct blocks; those break the factorisation that
the argument rests on, and are excluded from the first campaign.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

import numpy as np

from .measurement import (block_quadrature, classical_fisher_matrix, equatorial_basis_change,
                          mixed_state_qfim)
from .paulis import PauliSet, apply_pauli, z_generators
from .width import max_weight_matching

_I2 = np.eye(2, dtype=complex)
_PAULIS = {
    "I": _I2,
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _kron(*ops):
    out = np.ones((1, 1), dtype=complex)
    for o in ops:
        out = np.kron(o, out)          # little-endian, matching covq.sim
    return out


@dataclass
class BlockLocalNoise:
    """A noise model that respects block boundaries.

    ``dephasing`` and ``depolarizing`` are per-qubit and may be heterogeneous;
    ``edge_depolarizing`` is two-qubit depolarization applied after a pair
    preparation, so it is charged only to edges that a branch actually
    activates; ``idle_dephasing`` is charged to qubits left unmatched while the
    pairs are prepared; ``readout_confusion`` maps qubit -> ``(e0, e1)``, the
    asymmetric probabilities of misreporting each outcome.
    """

    dephasing: dict[int, float] = field(default_factory=dict)
    depolarizing: dict[int, float] = field(default_factory=dict)
    edge_depolarizing: dict[tuple[int, int], float] = field(default_factory=dict)
    idle_dephasing: float = 0.0
    readout_confusion: dict[int, tuple[float, float]] = field(default_factory=dict)

    def q_dephase(self, q: int) -> float:
        return float(self.dephasing.get(q, 0.0))

    def q_depol(self, q: int) -> float:
        return float(self.depolarizing.get(q, 0.0))

    def e_depol(self, e) -> float:
        return float(self.edge_depolarizing.get(tuple(sorted(e)), 0.0))

    def confusion(self, q: int) -> tuple[float, float]:
        return tuple(self.readout_confusion.get(q, (0.0, 0.0)))

    @property
    def is_symmetric_readout(self) -> bool:
        return all(abs(a - b) < 1e-15 for a, b in self.readout_confusion.values())


# ----------------------------------------------------------------------
# Block channels
# ----------------------------------------------------------------------

def _apply_1q_channel(rho: np.ndarray, k: int, nq: int, p_deph: float, p_depol: float):
    if p_deph > 0:
        z = _embed(_PAULIS["Z"], k, nq)
        rho = (1 - p_deph) * rho + p_deph * (z @ rho @ z.conj().T)
    if p_depol > 0:
        acc = (1 - p_depol) * rho
        for name in ("X", "Y", "Z"):
            u = _embed(_PAULIS[name], k, nq)
            acc = acc + (p_depol / 3.0) * (u @ rho @ u.conj().T)
        rho = acc
    return rho


def _embed(op: np.ndarray, k: int, nq: int) -> np.ndarray:
    return _kron(*[op if i == k else _I2 for i in range(nq)])


def _two_qubit_depolarize(rho: np.ndarray, q: float) -> np.ndarray:
    if q <= 0:
        return rho
    return (1 - q) * rho + q * np.trace(rho).real * np.eye(4, dtype=complex) / 4.0


def _confusion_matrix(qubits, noise: BlockLocalNoise) -> np.ndarray:
    """``C[obs, true]`` for independent asymmetric readout error."""
    mats = []
    for q in qubits:
        e0, e1 = noise.confusion(q)
        mats.append(np.array([[1 - e0, e1], [e0, 1 - e1]]))
    out = np.ones((1, 1))
    for mtx in mats:
        out = np.kron(mtx, out)
        
    return out


def _block_cfi(rho: np.ndarray, ps: PauliSet, alphas, qubits, noise: BlockLocalNoise):
    """CFI of the declared product readout, with readout confusion applied."""
    dim = rho.shape[0]
    basis = np.eye(dim, dtype=complex)
    nq = int(np.log2(dim))
    b = np.eye(dim, dtype=complex)
    for k, a in enumerate(alphas):
        b = _apply_basis(b, equatorial_basis_change(float(a)), k, nq)
    rho_b = b @ rho @ b.conj().T
    p = np.real(np.diag(rho_b)).copy()
    dp = np.empty((ps.m, p.size))
    for i in range(ps.m):
        col = np.column_stack([apply_pauli(basis[:, k], ps, i) for k in range(dim)])
        pr = b @ (col @ rho) @ b.conj().T
        dp[i] = np.imag(np.diag(pr))
    conf = _confusion_matrix(qubits, noise)
    return classical_fisher_matrix(conf @ p, dp @ conf.T)


def _apply_basis(mat: np.ndarray, u: np.ndarray, k: int, nq: int) -> np.ndarray:
    return _embed(u, k, nq) @ mat


def singleton_channel(theta: float, q: int, noise: BlockLocalNoise, idle: bool):
    """``|+>`` under ``Z`` rotation and local noise."""
    psi = np.array([1.0, 1.0], dtype=complex) / math.sqrt(2.0)
    psi = np.array([np.exp(-1j * theta / 2) * psi[0], np.exp(1j * theta / 2) * psi[1]])
    rho = np.outer(psi, psi.conj())
    p_deph = noise.q_dephase(q) + (noise.idle_dephasing if idle else 0.0)
    return _apply_1q_channel(rho, 0, 1, min(p_deph, 0.5), noise.q_depol(q))


def pair_channel(thetas, e, sign: int, noise: BlockLocalNoise):
    """Signed Bell block under ``Z`` rotations, pair depolarization, local noise."""
    psi = np.zeros(4, dtype=complex)
    # |s> and |-s> for s = (+1, sign).  Qubit q holds bit (1 - s_q)/2 and the
    # index is little-endian, so s = (+1,+1) -> {0b00, 0b11} while
    # s = (+1,-1) -> {0b10, 0b01}.  Getting this wrong yields a state that is
    # not a cat at all and silently produces a rank-one, asymmetric block.
    lo = 0b00 if sign > 0 else 0b10
    hi = 0b11 if sign > 0 else 0b01
    psi[lo] = 1 / math.sqrt(2.0)
    psi[hi] = 1 / math.sqrt(2.0)
    phase = np.array([np.exp(-1j * (thetas[0] * (1 - 2 * ((k >> 0) & 1))
                                    + thetas[1] * (1 - 2 * ((k >> 1) & 1))) / 2)
                      for k in range(4)])
    psi = psi * phase
    rho = np.outer(psi, psi.conj())
    rho = _two_qubit_depolarize(rho, noise.e_depol(e))
    for k, q in enumerate(e):
        rho = _apply_1q_channel(rho, k, 2, noise.q_dephase(q), noise.q_depol(q))
    return rho


def singleton_template(theta: float, q: int, noise: BlockLocalNoise,
                       idle: bool = True, three_point: bool | None = None) -> float:
    rho = singleton_channel(theta, q, noise, idle)
    ps = z_generators(1)
    tp = (not noise.is_symmetric_readout) if three_point is None else three_point
    rec = block_quadrature(rho, [0], 1, three_point=tp)
    return float(_block_cfi(rho, ps, [rec["alpha"]], [q], noise)[0, 0])


def pair_template(thetas, e, sign: int, noise: BlockLocalNoise,
                  three_point: bool | None = None) -> np.ndarray:
    rho = pair_channel(thetas, e, sign, noise)
    ps = z_generators(2)
    tp = (not noise.is_symmetric_readout) if three_point is None else three_point
    rec = block_quadrature(rho, [0, 1], 2, three_point=tp)
    alphas = [rec["alpha"], 0.0]
    return _block_cfi(rho, ps, alphas, list(e), noise)


def branch_template(matching, signs, m: int, theta, noise: BlockLocalNoise) -> np.ndarray:
    """Full ``m x m`` declared-readout CFI of one noisy branch."""
    F = np.zeros((m, m))
    matched = set()
    for e, s in zip(matching, signs):
        i, j = int(e[0]), int(e[1])
        blk = pair_template((theta[i], theta[j]), (i, j), int(s), noise)
        F[np.ix_([i, j], [i, j])] += blk
        matched.update((i, j))
    for q in range(m):
        if q not in matched:
            F[q, q] += singleton_template(theta[q], q, noise, idle=bool(matching))
    return F


def product_template(m: int, theta, noise: BlockLocalNoise,
                     idle: bool = False) -> np.ndarray:
    """``F_0^N``: every qubit a singleton.

    ``idle`` selects the reference.  The pairless branch genuinely does not wait
    for a pair preparation, so it is idle-free; but the reference used by the
    *pricing* decomposition must idle, because every non-empty pair-width branch
    has two-qubit depth exactly one and so imposes the same wait on its
    unmatched qubits.  Using the idle-free reference inside pricing silently
    breaks edge additivity as soon as ``idle_dephasing > 0``.
    """
    F = np.zeros((m, m))
    for q in range(m):
        F[q, q] = singleton_template(theta[q], q, noise, idle=idle)
    return F


def delta_edge(e, sign: int, m: int, theta, noise: BlockLocalNoise) -> np.ndarray:
    """``DeltaF_{e,sigma}^N``: activating one edge, against the product template."""
    i, j = int(e[0]), int(e[1])
    D = np.zeros((m, m))
    D[np.ix_([i, j], [i, j])] += pair_template((theta[i], theta[j]), (i, j), sign, noise)
    D[i, i] -= singleton_template(theta[i], i, noise, idle=True)
    D[j, j] -= singleton_template(theta[j], j, noise, idle=True)
    return D


# ----------------------------------------------------------------------
# The pricing theorem
# ----------------------------------------------------------------------

def edge_weight(Q: np.ndarray, e, m: int, theta, noise: BlockLocalNoise,
                cost: float) -> tuple[float, int]:
    """``max_sigma <Q, DeltaF_{e,sigma}^N> - c_e`` and the winning sign."""
    best_val, best_sign = -np.inf, 1
    for sign in (1, -1):
        val = float(np.sum(Q * delta_edge(e, sign, m, theta, noise))) - cost
        if val > best_val:
            best_val, best_sign = val, sign
    return best_val, best_sign


def price_branch_noisy(Q: np.ndarray, m: int, edges, theta, noise: BlockLocalNoise,
                       costs=None, c0: float = 0.0) -> dict:
    """Noisy analogue of the Eq (65) pricing step.

    Reduces to a maximum-weight matching on reweighted edges, which is the
    content of the block-local pricing theorem.  The reference template is the
    *idling* product branch, because every pair-width branch has two-qubit
    depth exactly one and therefore imposes the same idle wait on its unmatched
    qubits; the pairless branch is priced separately, since it does not wait.
    """
    costs = {} if costs is None else costs
    base = product_template(m, theta, noise, idle=True)
    W = np.zeros((m, m))
    signs: dict[tuple[int, int], int] = {}
    for e in edges:
        i, j = (int(e[0]), int(e[1])) if e[0] < e[1] else (int(e[1]), int(e[0]))
        w, s = edge_weight(Q, (i, j), m, theta, noise, float(costs.get((i, j), 1.0)))
        W[i, j] = W[j, i] = max(w, 0.0)
        signs[(i, j)] = s
    matching, _ = max_weight_matching(W)
    chosen = [e for e in matching if W[e[0], e[1]] > 0.0]
    paired = float(np.sum(Q * base)) - c0 + sum(W[i, j] for i, j in chosen)

    # The pairless branch is a *separate* column, not the zero-edge case of the
    # matching problem.  It activates no pair, so it never waits through a pair
    # preparation and carries no idle dephasing, which puts it on a different
    # template from the idling reference the edge weights are measured against.
    # With idle_dephasing = 0 the two references coincide and this reduces to
    # one matching problem; with idle noise present, dropping this candidate
    # loses the optimum whenever idling costs more than the best edge gains.
    pairless = float(np.sum(Q * product_template(m, theta, noise, idle=False))) - c0

    if pairless > paired:
        return {"value": pairless, "matching": [], "signs": [],
                "base_value": pairless, "branch_kind": "pairless_idle_free"}
    return {"value": paired, "matching": chosen,
            "signs": [signs[tuple(sorted(e))] for e in chosen],
            "base_value": float(np.sum(Q * base)) - c0,
            "branch_kind": "matching_from_idling_reference"}


def price_branch_bruteforce(Q: np.ndarray, m: int, edges, theta,
                            noise: BlockLocalNoise, costs=None,
                            c0: float = 0.0) -> dict:
    """Exhaustive enumeration of signed matchings, for the N5 gate."""
    costs = {} if costs is None else costs
    norm = [tuple(sorted((int(a), int(b)))) for a, b in edges]
    best = {"value": -np.inf, "matching": [], "signs": []}
    for r in range(len(norm) + 1):
        for sub in itertools.combinations(norm, r):
            seen: set[int] = set()
            if any(q in seen or seen.add(q) for e in sub for q in e):
                continue
            for sgn in itertools.product((1, -1), repeat=r):
                F = branch_template(list(sub), list(sgn), m, theta, noise)
                val = float(np.sum(Q * F)) - c0 - sum(costs.get(e, 1.0) for e in sub)
                if val > best["value"]:
                    best = {"value": val, "matching": list(sub), "signs": list(sgn)}
    return best


# ----------------------------------------------------------------------
# Eq (111): the noisy operational information-floor compiler
# ----------------------------------------------------------------------

def noise_aware_floor_compile(G_req: np.ndarray, m: int, edges, theta,
                              noise: BlockLocalNoise, A=None, c0: float = 1.0,
                              costs=None, tol: float = 1e-8,
                              max_rounds: int = 120, max_cuts: int = 300) -> dict:
    """Minimum-cost exposure meeting a floor on the **declared readout** CFI.

        min_n  sum_b n_b c_b   s.t.   A^T (sum_b n_b F_{C,b}^N) A >= G_req,  n >= 0.

    Convex, and solved by exactly the ideal machinery: a Loewner cutting-plane
    master whose aggregated multipliers form the PSD dual `Y`, and a pricing
    step that is a maximum-weight matching on noise-reweighted edges.  That the
    oracle survives the move to noisy templates is the block-local pricing
    theorem; that it stays *convex* is why the dual bound below is meaningful.

    Per-use costs only.  A fixed charge per distinct setting is a cardinality
    term, which would make this a mixed-integer conic program and void the
    oracle-polynomial claim; that variant is handled separately in
    :func:`setting_cost_crossover` and is explicitly not solved here.
    """
    from scipy.optimize import linprog

    G_req = np.atleast_2d(np.asarray(G_req, dtype=float))
    A = np.eye(m) if A is None else np.asarray(A, dtype=float)
    E = sorted({tuple(sorted((int(a), int(b)))) for a, b in edges})
    ce = {e: 1.0 for e in E} if costs is None else {
        e: float(costs[e] if isinstance(costs, dict) else costs[k]) for k, e in enumerate(E)}

    def template(M, s):
        return (product_template(m, theta, noise, idle=False) if not M
                else branch_template(list(M), list(s), m, theta, noise))

    cols = [(([], []), c0)]
    Y = np.zeros(G_req.shape)
    history, n_solved, hit_limit, t = [], 1, True, None
    for rnd in range(max_rounds):
        Bs = [A.T @ template(M, s) @ A for (M, s), _ in cols]
        cvec = np.array([c for _, c in cols])
        cuts_v: list[np.ndarray] = []
        for _ in range(max_cuts):
            if cuts_v:
                rows = np.array([[-float(v @ B @ v) for B in Bs] for v in cuts_v])
                rhs = np.array([-float(v @ G_req @ v) for v in cuts_v])
            else:
                rows, rhs = np.zeros((1, len(Bs))), np.array([0.0])
            res = linprog(cvec, A_ub=rows, b_ub=rhs, bounds=(0.0, None), method="highs")
            if not res.success:
                return {"status": "master_infeasible", "round": rnd,
                        "n_branches": len(cols)}
            t = np.asarray(res.x, dtype=float)
            w, V = np.linalg.eigh(sum(ti * B for ti, B in zip(t, Bs)) - G_req)
            if w[0] >= -tol:
                break
            cuts_v.append(V[:, 0].copy())
        n_solved = len(cols)
        eta = -np.asarray(res.ineqlin.marginals, dtype=float) if cuts_v else None
        Y = np.asarray(sum(e * np.outer(v, v) for e, v in zip(eta, cuts_v))
                       if cuts_v else np.zeros(G_req.shape), dtype=float)
        priced = price_branch_noisy(A @ Y @ A.T, m, E, theta, noise, ce, c0)
        reduced = priced["value"]
        history.append({"round": rnd, "cost": float(cvec @ t),
                        "dual": float(np.sum(G_req * Y)), "reduced_cost": reduced,
                        "n_branches": len(cols)})
        if reduced <= tol:
            hit_limit = False
            break
        key = (tuple(map(tuple, priced["matching"])), tuple(priced["signs"]))
        if any(key == (tuple(map(tuple, M)), tuple(s)) for (M, s), _ in cols):
            hit_limit = False           # oracle repeated a column: nothing left to add
            break
        cols.append(((priced["matching"], priced["signs"]),
                     c0 + sum(ce[tuple(sorted(e))] for e in priced["matching"])))

    cols = cols[:n_solved]
    t = t[:n_solved]
    total = sum(ti * template(M, s) for ti, ((M, s), _) in zip(t, cols))
    cost = float(np.array([c for _, c in cols]) @ t)
    dual = float(np.sum(G_req * Y))
    return {
        "status": "round_limit" if hit_limit else "solved",
        "cost": cost, "dual_bound": dual,
        "relative_gap": abs(cost - dual) / max(abs(cost), 1e-12),
        "Y": Y, "exposures": t,
        "branches": [(M, s) for (M, s), _ in cols],
        "total_information": total,
        "floor_slack_min_eig": float(np.linalg.eigvalsh(A.T @ total @ A - G_req).min()),
        "n_settings_used": int((t > 1e-10).sum()),
        "n_entangled_settings_used": int(sum(1 for ti, ((M, _), _) in zip(t, cols)
                                             if ti > 1e-10 and len(M) > 0)),
        "history": history,
    }


def setting_cost_crossover(cost_covq: float, q_covq: int,
                           cost_quest: float, q_quest: int,
                           c_setup: float) -> dict:
    """Amortization boundary between a many-setting schedule and a one-state probe.

    With batched execution the number of switches is ``q - 1``, not the number
    of shots, so the totals are affine in ``N``:

        C(N) = q c_setup + N c_bar,

    and when the per-shot cost of the one-state probe exceeds the schedule's,

        N* = (q_C - q_Q) c_setup / (c_Q - c_bar_C)

    is the repetition budget beyond which the schedule wins.  This converts the
    unpriced-settings caveat into a measurable boundary rather than leaving it
    as a hedge.
    """
    denom = cost_quest - cost_covq
    if denom <= 0:
        return {"crossover_shots": None,
                "verdict": "schedule never cheaper: its per-shot cost is not lower"}
    n_star = (q_covq - q_quest) * c_setup / denom
    return {"crossover_shots": float(max(n_star, 0.0)),
            "per_shot_saving": float(denom),
            "setup_penalty": float((q_covq - q_quest) * c_setup),
            "verdict": ("schedule cheaper for every N" if n_star <= 0
                        else f"schedule cheaper once N > {n_star:.1f}")}
