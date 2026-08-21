"""QUEST baseline: expectation-value targeting by depth-adaptive Pauli rotations.

Two implementations live here, and only one of them is QUEST.

:func:`quest_published` implements the algorithm of arXiv:2605.02367 --
Mahapatra and Kadiri, *Quantum State Engineering Under Multiple
Expectation-Value Constraints*; QUEST expands to "Quantum Unitary Engineering
of States to Target".  Every published variant runs **two phases per
iteration**: insert one Pauli rotation, then **jointly reoptimise all
accumulated angles** by classical optimisation (L-BFGS in the paper's numerical
study).  ``variant='tE'`` inserts only at the terminal position;
``variant='bE'`` searches every insertion position.  Written from the published
method description rather than from released source.

:func:`quest` is the *earlier* routine and is **not QUEST**.  It appends only at
the terminal position and never reoptimises earlier angles, so it is a terminal
greedy Pauli-path expectation-targeting baseline inspired by QUEST.  It was
originally described here as a faithful reimplementation; that was wrong, and
the fidelity audit in ``docs/audits/QUEST_FIDELITY_AUDIT_v0.4.md`` is what
caught it.  The omission is not cosmetic: on the registered ``m = 4`` contract
target the greedy routine needs 50 rotations and 44 two-qubit gates where
QUEST-tE needs 4 and 3.  Any resource comparison built on the greedy routine
understates the baseline by an order of magnitude, so it is retained only as a
labelled ablation and must never be reported as QUEST.

Why this is the right comparator, and only here.  CovQ's *exact-target* mode
asks for a state whose first and second generator moments hit prescribed
values.  That is an instance of expectation-value targeting, which is why the
K1a specification collision is real.  It is **not** an information-floor
optimizer: QUEST is given a single target point, whereas the floor compiler is
given a Loewner constraint `A^T F A >= G_req` and searches over every `F` that
satisfies it.  Per the v0.3 handoff, QUEST must not be forced into that primary
role without a separately specified outer optimization over targets.  This
module therefore compares on exact targets only.

The line search is exact rather than iterative.  For a Pauli rotation
`R(t) = exp(-i t P / 2)` and any observable `O`,

    <R(t)^dag O R(t)> = a + b cos t + c sin t,

identically, because `P^2 = I`.  Three evaluations determine `(a, b, c)` per
observable, so the residual cost is a known degree-two trigonometric polynomial
in `t` and its minimizer is found on a dense grid plus a parabolic refinement.
No gradients, no barren plateau, and no optimizer state -- which is the
character of the published method.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .circuits import Circuit, lower_to_cx, resources
from .paulis import PauliSet, apply_pauli
from .qfim import qfim_from_statevector
from .sim import simulate

_PAULI = {
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


@dataclass(frozen=True)
class PauliRotation:
    """``exp(-i t/2 * prod_q sigma_{axis[q]}^{(q)})``."""

    axes: tuple[tuple[int, str], ...]

    @property
    def weight(self) -> int:
        return len(self.axes)

    def label(self) -> str:
        return "".join(f"{a}{q}" for q, a in self.axes)


def _dense(axes, n: int) -> np.ndarray:
    op = np.ones((1, 1), dtype=complex)
    table = {q: a for q, a in axes}
    for q in range(n):
        op = np.kron(_PAULI[table[q]] if q in table else np.eye(2, dtype=complex), op)
    return op


def _apply_rotation(psi: np.ndarray, dense: np.ndarray, t: float) -> np.ndarray:
    return math.cos(t / 2.0) * psi - 1j * math.sin(t / 2.0) * (dense @ psi)


def default_pool(n: int, edges=None) -> list[PauliRotation]:
    """Single-qubit rotations plus two-qubit Pauli rotations on the coupling graph."""
    pool = [PauliRotation(((q, a),)) for q in range(n) for a in ("X", "Y", "Z")]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)] if edges is None else edges
    for i, j in edges:
        i, j = (int(i), int(j)) if i < j else (int(j), int(i))
        for a in ("X", "Y", "Z"):
            for b in ("X", "Y", "Z"):
                pool.append(PauliRotation(((i, a), (j, b))))
    return pool


@dataclass
class QuestResult:
    converged: bool
    residual: float
    rotations: list[tuple[PauliRotation, float]]
    state: np.ndarray
    history: list[float] = field(default_factory=list)
    stop_reason: str = ""

    @property
    def depth_adaptive_length(self) -> int:
        return len(self.rotations)

    @property
    def two_qubit_rotations(self) -> int:
        return sum(1 for r, _ in self.rotations if r.weight == 2)


def moment_constraints(ps: PauliSet, F: np.ndarray) -> list[tuple[np.ndarray, float]]:
    """Exact first/second generator moments for a zero-mean unit-diagonal target.

    ``<P_i> = 0`` and ``<P_i P_j> = F_ij``: the constraint set CovQ's exact mode
    and QUEST are both asked to satisfy, stated identically for both.
    """
    m, n = ps.m, ps.n
    out = []
    basis = np.eye(1 << n, dtype=complex)
    mats = []
    for i in range(m):
        mats.append(np.column_stack([apply_pauli(basis[:, k], ps, i) for k in range(1 << n)]))
    for i in range(m):
        out.append((mats[i], 0.0))
    for i in range(m):
        for j in range(i + 1, m):
            out.append((mats[i] @ mats[j], float(F[i, j])))
    return out


def _residual(psi: np.ndarray, constraints) -> float:
    return float(sum((np.real(np.vdot(psi, O @ psi)) - tau) ** 2 for O, tau in constraints))


def quest(constraints, n: int, pool=None, max_depth: int = 64,
          tol: float = 1e-12, psi0: np.ndarray | None = None,
          grid: int = 256) -> QuestResult:
    """Depth-adaptive greedy Pauli-rotation synthesis for expectation targeting."""
    pool = default_pool(n) if pool is None else pool
    dense_pool = [(r, _dense(r.axes, n)) for r in pool]
    dense_obs = [(np.asarray(O), tau) for O, tau in constraints]

    # Default start is |+>^n, not |0..0>.  For zero-mean unit-diagonal moment
    # targets |0..0> is the worst possible seed: it is a joint Z eigenstate, so
    # every <Z_i> is maximally wrong and every Z-type rotation in the pool is a
    # no-op on it.  |+>^n sits at the origin of the moment space -- all first
    # and second generator moments vanish -- and the greedy descent converges
    # from there on targets where it stalls from |0..0>.  Reporting the |0..0>
    # numbers as this method's performance would understate it badly.
    if psi0 is None:
        psi = np.ones(1 << n, dtype=complex) / math.sqrt(1 << n)
    else:
        psi = psi0.astype(complex).copy()
    ts = np.linspace(-math.pi, math.pi, grid, endpoint=False)
    cos_t, sin_t = np.cos(ts), np.sin(ts)

    rotations: list[tuple[PauliRotation, float]] = []
    history = [_residual(psi, dense_obs)]
    stop = "max_depth"
    for _ in range(max_depth):
        if history[-1] <= tol:
            stop = "converged"
            break
        best = None
        for rot, dense in dense_pool:
            # Exact three-point trigonometric fit, per observable.
            p0, pp, pm = psi, _apply_rotation(psi, dense, math.pi / 2), \
                _apply_rotation(psi, dense, -math.pi / 2)
            cost = np.zeros(grid)
            for O, tau in dense_obs:
                v0 = float(np.real(np.vdot(p0, O @ p0)))
                vp = float(np.real(np.vdot(pp, O @ pp)))
                vm = float(np.real(np.vdot(pm, O @ pm)))
                a = 0.5 * (vp + vm)
                b = v0 - a
                c = 0.5 * (vp - vm)
                cost += (a + b * cos_t + c * sin_t - tau) ** 2
            k = int(np.argmin(cost))
            t_star, c_star = _refine(ts, cost, k, dense_obs, psi, dense)
            if best is None or c_star < best[2]:
                best = (rot, t_star, c_star, dense)
        rot, t_star, c_star, dense = best
        if c_star >= history[-1] - 1e-15:
            stop = "stalled"
            break
        psi = _apply_rotation(psi, dense, t_star)
        rotations.append((rot, float(t_star)))
        history.append(c_star)

    return QuestResult(history[-1] <= tol, history[-1], rotations, psi, history,
                       "converged" if history[-1] <= tol else stop)


def _refine(ts, cost, k, dense_obs, psi, dense, iters: int = 40):
    """Golden-section refinement inside the winning grid cell."""
    n = ts.size
    lo, hi = ts[(k - 1) % n], ts[(k + 1) % n]
    if hi < lo:
        hi += 2 * math.pi
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc = _residual(_apply_rotation(psi, dense, c), dense_obs)
    fd = _residual(_apply_rotation(psi, dense, d), dense_obs)
    for _ in range(iters):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = _residual(_apply_rotation(psi, dense, c), dense_obs)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = _residual(_apply_rotation(psi, dense, d), dense_obs)
    t = c if fc < fd else d
    return float(t), float(min(fc, fd))


def canonicalise_rotations(rotations, angle_tol: float = 1e-8):
    """Apply the two lowering rules that both backends must share.

    *Near-zero-angle removal*: a rotation with ``|theta| <= angle_tol`` is the
    identity and must not be billed for two CX gates.  *Adjacent-rotation
    combination*: consecutive rotations about the same Pauli commute into one.

    Both are no-ops on the registered QUEST solutions -- no angle falls below
    tolerance and no two adjacent rotations share a Pauli -- and CovQ's
    preparation contains no parameterised rotations at all.  They are
    implemented anyway so the resource audit can state that identical
    simplification was available to both sides, rather than leaving it silently
    absent from one.
    """
    out: list = []
    for rot, theta in rotations:
        if out and out[-1][0].axes == rot.axes:
            merged = out[-1][1] + theta
            if abs(merged) <= angle_tol:
                out.pop()
            else:
                out[-1] = (rot, merged)
            continue
        if abs(theta) <= angle_tol:
            continue
        out.append((rot, float(theta)))
    return out


def quest_circuit(result: QuestResult, n: int, angle_tol: float = 1e-8) -> Circuit:
    """Emit the rotation sequence so resources are parsed, never estimated.

    Each two-qubit Pauli rotation lowers to basis changes, two CX, and one RZ;
    single-qubit rotations lower to basis changes and one RZ.  Gate 9 of the
    handoff requires resource claims to come from emitted artifacts, so the
    comparison below counts these, not a formula.
    """
    c = Circuit(n)
    # The |+>^n preparation is part of the emitted program and must be in the
    # circuit.  Without it the artifact is not self-contained: simulating it
    # starts from |0..0> and produces a different state, and its resource count
    # omits the preparation layer.  It adds no two-qubit gate, so the CX-based
    # comparisons are unchanged -- but a circuit that does not prepare its own
    # input is not an emitted program.
    for q in range(n):
        c.h(q)
    for rot, t in canonicalise_rotations(result.rotations, angle_tol):
        pre = []
        for q, a in rot.axes:
            if a == "X":
                c.h(q); pre.append(("h", q))
            elif a == "Y":
                c.sdg(q); c.h(q); pre.append(("y", q))
        qs = [q for q, _ in rot.axes]
        if len(qs) == 2:
            c.cx(qs[0], qs[1])
            c.rz(qs[1], t)
            c.cx(qs[0], qs[1])
        else:
            c.rz(qs[0], t)
        for kind, q in reversed(pre):
            if kind == "h":
                c.h(q)
            else:
                c.h(q); c.s(q)
    return c


# ----------------------------------------------------------------------
# The published algorithm: insertion + joint reoptimisation
# ----------------------------------------------------------------------

def _forward(rots, angles, psi0):
    """States after each prefix: ``phi[k] = R_k ... R_1 |psi0>``, with ``phi[0] = psi0``."""
    out = [psi0]
    cur = psi0
    for (_, dense), t in zip(rots, angles):
        cur = _apply_rotation(cur, dense, t)
        out.append(cur)
    return out


def _cost_and_grad(angles, rots, psi0, obs):
    """Sum-of-squared-residuals cost and its exact gradient.

    The gradient is adjoint, not finite-difference and not parameter-shift.
    Writing ``W = sum_c 2 (<O_c> - tau_c) O_c``, the chain rule collapses every
    observable into one operator, so a single backward sweep gives every angle
    derivative:

        dC/dtheta_k = 2 Re <chi_k| (-i P_k / 2) |phi_k>,

    with ``phi_k`` the forward prefix state and ``chi_k`` obtained by pulling
    ``W|psi>`` back through the later rotations.  Cost is one forward and one
    backward pass regardless of how many observables or angles there are, which
    is what makes joint L-BFGS over all accumulated angles affordable at every
    insertion.
    """
    phis = _forward(rots, angles, psi0)
    psi = phis[-1]
    resid = []
    w = np.zeros_like(psi)
    for O, tau in obs:
        Opsi = O @ psi
        v = float(np.real(np.vdot(psi, Opsi)))
        resid.append((v - tau) ** 2)
        w = w + 2.0 * (v - tau) * Opsi
    cost = float(sum(resid))

    grad = np.zeros(len(angles))
    chi = w
    for k in range(len(angles) - 1, -1, -1):
        _, dense = rots[k]
        # Undo R_k on the adjoint state, then contract with the generator.
        chi = _apply_rotation(chi, dense, -angles[k])
        dpsi = -0.5j * (dense @ phis[k])
        grad[k] = 2.0 * float(np.real(np.vdot(chi, dpsi)))
    return cost, grad


def _best_single_angle(rots, angles, psi0, obs, cand_dense, position, grid):
    """Exact one-angle minimum for inserting ``cand_dense`` at ``position``."""
    trial_rots = rots[:position] + [(None, cand_dense)] + rots[position:]
    base = list(angles[:position]) + [0.0] + list(angles[position:])
    ts = np.linspace(-math.pi, math.pi, grid, endpoint=False)
    best_t, best_c = 0.0, np.inf
    for t in ts:
        base[position] = float(t)
        c, _ = _cost_and_grad(base, trial_rots, psi0, obs)
        if c < best_c:
            best_t, best_c = float(t), c
    lo, hi = best_t - math.pi / grid, best_t + math.pi / grid
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c1, d1 = b - phi * (b - a), a + phi * (b - a)

    def f(t):
        base[position] = float(t)
        return _cost_and_grad(base, trial_rots, psi0, obs)[0]

    fc, fd = f(c1), f(d1)
    for _ in range(40):
        if fc < fd:
            b, d1, fd = d1, c1, fc
            c1 = b - phi * (b - a)
            fc = f(c1)
        else:
            a, c1, fc = c1, d1, fd
            d1 = a + phi * (b - a)
            fd = f(d1)
    t = c1 if fc < fd else d1
    return float(t), float(min(fc, fd))


def quest_published(constraints, n: int, pool=None, variant: str = "tE",
                    max_depth: int = 64, tol: float = 1e-12,
                    psi0: np.ndarray | None = None, grid: int = 64,
                    maxiter: int = 200, position_window: int | None = None) -> QuestResult:
    """QUEST as published (arXiv:2605.02367): insert, then jointly reoptimise.

    Each iteration has the paper's two phases:

    1. **insert one Pauli rotation**, chosen by the exact one-angle minimum of
       the sum-of-squared-residuals cost.  ``variant='tE'`` searches only the
       terminal position; ``variant='bE'`` searches every insertion position,
       which the paper reports as the faster exact variant.
    2. **jointly reoptimise every accumulated angle** by L-BFGS.

    Phase 2 is the load-bearing difference from a terminal greedy Pauli-path
    method, and omitting it is what made the earlier implementation in this file
    a QUEST-*inspired* baseline rather than QUEST.  Later rotations routinely
    make earlier angles suboptimal; without reoptimisation the depth needed to
    hit a target is an overestimate, and any resource comparison built on it
    understates the baseline.
    """
    from scipy.optimize import minimize

    if variant not in ("tE", "bE"):
        raise ValueError("variant must be 'tE' (terminal exact) or 'bE' (best-position)")
    if position_window is not None and variant != "bE":
        raise ValueError("position_window only applies to the bE variant")
    pool = default_pool(n) if pool is None else pool
    dense_pool = [(r, _dense(r.axes, n)) for r in pool]
    obs = [(np.asarray(O), float(tau)) for O, tau in constraints]

    if psi0 is None:
        psi0 = np.ones(1 << n, dtype=complex) / math.sqrt(1 << n)
    else:
        psi0 = psi0.astype(complex)

    rots: list = []
    angles: list[float] = []
    history = [float(sum((np.real(np.vdot(psi0, O @ psi0)) - tau) ** 2 for O, tau in obs))]
    stop = "max_depth"

    for _ in range(max_depth):
        if history[-1] <= tol:
            stop = "converged"
            break
        if variant == "tE":
            positions = [len(rots)]
        elif position_window is None:
            positions = list(range(len(rots) + 1))
        else:
            # Windowed best-position insertion: search only the last W slots.
            # W = 1 recovers tE and W = inf recovers bE, so the *cost*
            # interpolates between the two published variants, linear rather
            # than quadratic in depth.
            #
            # The solution quality does NOT interpolate, and it is worth saying
            # so plainly: on path_m4, W = 2 needs 152 emitted CX and fails to
            # converge where W = 1 converges at 24 and W = 4 at 20.  Greedy
            # insertion is not monotone in its candidate set -- widening the
            # window changes which rotation wins the first few rounds and can
            # steer the whole trajectory into a worse basin.  A windowed run is
            # therefore evidence about the instance, never a bound on bE, and
            # must never be reported as bE itself.
            lo = max(0, len(rots) + 1 - int(position_window))
            positions = list(range(lo, len(rots) + 1))
        best = None
        for rot, dense in dense_pool:
            for pos in positions:
                t, c = _best_single_angle(rots, angles, psi0, obs, dense, pos, grid)
                if best is None or c < best[0]:
                    best = (c, rot, dense, pos, t)
        _, rot, dense, pos, t = best
        rots.insert(pos, (rot, dense))
        angles.insert(pos, t)

        # Phase 2: joint reoptimisation of every accumulated angle.
        res = minimize(lambda a: _cost_and_grad(list(a), rots, psi0, obs), np.array(angles),
                       jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
        if res.fun <= _cost_and_grad(angles, rots, psi0, obs)[0]:
            angles = [float(x) for x in res.x]
        cost = _cost_and_grad(angles, rots, psi0, obs)[0]
        if cost >= history[-1] - 1e-15:
            history.append(cost)
            stop = "stalled"
            break
        history.append(cost)

    psi = _forward(rots, angles, psi0)[-1]
    return QuestResult(history[-1] <= tol, history[-1],
                       [(r, float(t)) for (r, _), t in zip(rots, angles)],
                       psi, history,
                       "converged" if history[-1] <= tol else stop)
