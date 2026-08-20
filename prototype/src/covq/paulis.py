"""Commuting Pauli generator programs and simultaneous Clifford normalisation.

A Hermitian Pauli is stored symplectically as ``P = sign * i^(x.z) X^x Z^z``
with ``x, z`` in F_2^n and ``sign`` in {+1, -1}.

The diagonalising Clifford is computed from the symplectic part alone; the
resulting signs are then read off *numerically* rather than by propagating a
phase rule, which removes an entire class of silent sign bugs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .circuits import Circuit
from .sim import simulate


def _popcount(a: np.ndarray) -> np.ndarray:
    try:
        return np.bitwise_count(a)
    except AttributeError:  # pragma: no cover - numpy < 2.0
        out = np.zeros_like(a)
        v = a.copy()
        while np.any(v):
            out += (v & 1).astype(out.dtype)
            v >>= 1
        return out


@dataclass
class PauliSet:
    """An ordered list of ``m`` Hermitian Pauli strings on ``n`` qubits."""

    x: np.ndarray  # (m, n) bool
    z: np.ndarray  # (m, n) bool
    sign: np.ndarray  # (m,) int8 in {+1, -1}

    def __post_init__(self) -> None:
        self.x = np.asarray(self.x, dtype=bool)
        self.z = np.asarray(self.z, dtype=bool)
        self.sign = np.asarray(self.sign, dtype=np.int8)
        if self.x.shape != self.z.shape:
            raise ValueError("x and z blocks must have the same shape")
        if self.sign.shape != (self.x.shape[0],):
            raise ValueError("one sign per generator")

    @property
    def m(self) -> int:
        return self.x.shape[0]

    @property
    def n(self) -> int:
        return self.x.shape[1]

    def label(self, i: int) -> str:
        chars = []
        for q in range(self.n):
            xi, zi = bool(self.x[i, q]), bool(self.z[i, q])
            chars.append({(0, 0): "I", (1, 0): "X", (0, 1): "Z", (1, 1): "Y"}[(int(xi), int(zi))])
        return ("-" if self.sign[i] < 0 else "+") + "".join(reversed(chars))

    def as_ints(self, i: int) -> tuple[int, int]:
        xi = int(np.dot(self.x[i].astype(np.int64), 1 << np.arange(self.n)))
        zi = int(np.dot(self.z[i].astype(np.int64), 1 << np.arange(self.n)))
        return xi, zi


def z_generators(m: int, n: int | None = None) -> PauliSet:
    """The physical-Z generator class: ``P_i = Z_i``."""
    n = m if n is None else n
    x = np.zeros((m, n), dtype=bool)
    z = np.zeros((m, n), dtype=bool)
    for i in range(m):
        z[i, i] = True
    return PauliSet(x, z, np.ones(m, dtype=np.int8))


def commutes(ps: PauliSet, i: int, j: int) -> bool:
    v = np.logical_and(ps.x[i], ps.z[j]).sum() + np.logical_and(ps.z[i], ps.x[j]).sum()
    return int(v) % 2 == 0


def all_commute(ps: PauliSet) -> bool:
    return all(commutes(ps, i, j) for i in range(ps.m) for j in range(i + 1, ps.m))


def _f2_rank(mat: np.ndarray) -> int:
    a = np.asarray(mat, dtype=np.uint8).copy() % 2
    rows, cols = a.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if a[i, c]:
                piv = i
                break
        if piv is None:
            continue
        a[[r, piv]] = a[[piv, r]]
        for i in range(rows):
            if i != r and a[i, c]:
                a[i] ^= a[r]
        r += 1
        if r == rows:
            break
    return r


def is_independent(ps: PauliSet) -> bool:
    return _f2_rank(np.hstack([ps.x, ps.z])) == ps.m


def apply_pauli(psi: np.ndarray, ps: PauliSet, i: int) -> np.ndarray:
    """Return ``P_i |psi>``."""
    n = ps.n
    xi, zi = ps.as_ints(i)
    idx = np.arange(1 << n, dtype=np.int64)
    phase = 1.0 - 2.0 * (_popcount(idx & zi) & 1)
    xz = int(np.logical_and(ps.x[i], ps.z[i]).sum())
    coef = complex(int(ps.sign[i])) * (1j) ** xz
    return coef * (phase * psi)[idx ^ xi]


# ----------------------------------------------------------------------
# Symplectic conjugation (phase-free)
# ----------------------------------------------------------------------

def _conj_h(x, z, q):
    x[:, q], z[:, q] = z[:, q].copy(), x[:, q].copy()


def _conj_s(x, z, q):
    z[:, q] ^= x[:, q]


def _conj_cx(x, z, c, t):
    x[:, t] ^= x[:, c]
    z[:, c] ^= z[:, t]


def _conj_cz(x, z, a, b):
    z[:, a] ^= x[:, b]
    z[:, b] ^= x[:, a]


def conjugate_symplectic(ps: PauliSet, circ: Circuit) -> tuple[np.ndarray, np.ndarray]:
    """Apply ``P -> U P U^dagger`` for a Clifford circuit, ignoring signs."""
    x, z = ps.x.copy(), ps.z.copy()
    for g in circ.gates:
        if g.name == "h":
            _conj_h(x, z, g.qubits[0])
        elif g.name == "s":
            _conj_s(x, z, g.qubits[0])
        elif g.name == "sdg":
            _conj_s(x, z, g.qubits[0])
        elif g.name == "x" or g.name == "z" or g.name == "y":
            pass  # Pauli conjugation only changes signs
        elif g.name == "cx":
            _conj_cx(x, z, *g.qubits)
        elif g.name == "cz":
            _conj_cz(x, z, *g.qubits)
        elif g.name == "swap":
            a, b = g.qubits
            x[:, [a, b]] = x[:, [b, a]]
            z[:, [a, b]] = z[:, [b, a]]
        else:
            raise ValueError(f"{g.name!r} is not a Clifford gate handled here")
    return x, z


@dataclass
class Normalisation:
    """Result of simultaneous diagonalisation.

    ``U P_i U^dagger = sign_out[i] * Z_{perm[i]}`` where ``U`` is the unitary of
    :attr:`circuit` (gates applied left to right).
    """

    circuit: Circuit
    perm: list[int]
    sign_out: np.ndarray
    cx_count: int

    def z_frame_target(self, F: np.ndarray) -> np.ndarray:
        """Map a physical-frame target QFIM into the Z frame.

        ``F_ij`` refers to generators ``P_i, P_j``; in the Z frame the same
        information sits at rows ``perm[i], perm[j]`` with a sign
        ``sign_out[i]*sign_out[j]``.
        """
        m = F.shape[0]
        G = np.empty_like(F)
        d = self.sign_out.astype(float)
        for i in range(m):
            for j in range(m):
                G[self.perm[i], self.perm[j]] = d[i] * d[j] * F[i, j]
        return G


def diagonalizing_clifford(ps: PauliSet) -> Normalisation:
    """Find a Clifford mapping independent commuting Paulis to signed ``Z_i``.

    Cost is O(m*n) two-qubit gates per generator in the worst case, i.e. the
    O(m^2) two-qubit overhead that any 'minimum entanglement' claim for a
    Clifford-conjugated generator class has to be reported against.
    """
    if not all_commute(ps):
        raise ValueError("generators do not pairwise commute")
    if not is_independent(ps):
        raise ValueError("generators are not independent")

    x, z = ps.x.copy(), ps.z.copy()
    m, n = ps.m, ps.n
    circ = Circuit(n)

    # Phase 1: make every row Z-type.
    for i in range(m):
        if not x[i].any():
            continue
        q = int(np.flatnonzero(x[i])[0])
        for q2 in np.flatnonzero(x[i]):
            q2 = int(q2)
            if q2 == q:
                continue
            circ.cx(q, q2)
            _conj_cx(x, z, q, q2)
        if z[i, q]:
            circ.s(q)
            _conj_s(x, z, q)
        for q2 in np.flatnonzero(z[i]):
            q2 = int(q2)
            if q2 == q:
                continue
            circ.cz(q, q2)
            _conj_cz(x, z, q, q2)
        assert not z[i].any() and x[i].sum() == 1 and x[i, q]
        circ.h(q)
        _conj_h(x, z, q)
    assert not x.any(), "phase 1 failed to remove all X components"

    # Phase 2: column-reduce the Z block so each row becomes a single Z.
    perm: list[int] = []
    for i in range(m):
        cand = [c for c in np.flatnonzero(z[i]) if int(c) not in perm]
        if not cand:
            raise RuntimeError("generators were not independent after phase 1")
        c = int(cand[0])
        for c2 in np.flatnonzero(z[i]):
            c2 = int(c2)
            if c2 == c:
                continue
            circ.cx(c2, c)  # column c2 += column c
            _conj_cx(x, z, c2, c)
        perm.append(c)
    for i in range(m):
        assert z[i].sum() == 1 and z[i, perm[i]]

    # Signs, read off numerically: sigma_i = <chi|P_i|chi> with |chi> = U^dag|0>.
    chi = simulate(circ.inverse())
    sign_out = np.empty(m, dtype=np.int8)
    for i in range(m):
        val = complex(np.vdot(chi, apply_pauli(chi, ps, i)))
        if abs(abs(val) - 1.0) > 1e-8 or abs(val.imag) > 1e-8:
            raise RuntimeError(f"generator {i} did not normalise to a signed Z ({val})")
        sign_out[i] = 1 if val.real > 0 else -1

    cx_count = sum(1 for g in circ.gates if g.name in {"cx", "cz"})
    return Normalisation(circ, perm, sign_out, cx_count)


def random_commuting_paulis(m: int, rng: np.random.Generator, n: int | None = None,
                            depth: int | None = None) -> PauliSet:
    """Random independent commuting Paulis: Clifford-conjugate ``Z_1..Z_m``."""
    n = m if n is None else n
    depth = 4 * n if depth is None else depth
    base = z_generators(m, n)
    c = Circuit(n)
    for _ in range(depth):
        kind = rng.integers(0, 3)
        if kind == 0:
            c.h(int(rng.integers(0, n)))
        elif kind == 1:
            c.s(int(rng.integers(0, n)))
        else:
            a, b = rng.choice(n, size=2, replace=False)
            c.cx(int(a), int(b))
    x, z = conjugate_symplectic(base, c)
    keep = np.any(x | z, axis=1)
    if not keep.all():  # pragma: no cover - conjugation is invertible
        raise RuntimeError("conjugation produced an identity generator")
    return PauliSet(x, z, np.ones(m, dtype=np.int8))
