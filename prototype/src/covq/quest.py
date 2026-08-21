"""QUEST baseline: expectation-value targeting by depth-adaptive Pauli rotations.

Reimplementation of the method of arXiv:2605.02367, *Quantum State Engineering
Under Multiple Expectation-Value Constraints* (QUEST, "Quantum Unitary
Engineering of States to Target"), from its published description: the
engineered state is built as a **depth-adaptive sequence of Pauli rotations**,
each rotation chosen to descend a **sum-of-squared-residuals** cost, with no
high-dimensional classical optimizer.  Written from the paper's method
description rather than from released source; it is a faithful reimplementation
of the stated algorithm, not a port, and should be read as such.

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


def quest_circuit(result: QuestResult, n: int) -> Circuit:
    """Emit the rotation sequence so resources are parsed, never estimated.

    Each two-qubit Pauli rotation lowers to basis changes, two CX, and one RZ;
    single-qubit rotations lower to basis changes and one RZ.  Gate 9 of the
    handoff requires resource claims to come from emitted artifacts, so the
    comparison below counts these, not a formula.
    """
    c = Circuit(n)
    for rot, t in result.rotations:
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
