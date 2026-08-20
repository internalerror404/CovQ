"""Dense statevector simulator for the CovQ circuit IR.

Little-endian: basis index ``b`` assigns qubit ``q`` the bit ``(b >> q) & 1``.
Exact for the widths in scope (m <= 12 data qubits plus flag/lowering ancillas).
"""

from __future__ import annotations

import math

import numpy as np

from .circuits import Circuit, Gate

_SQ = {
    "h": np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2.0),
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "z": np.array([[1, 0], [0, -1]], dtype=complex),
    "s": np.array([[1, 0], [0, 1j]], dtype=complex),
    "sdg": np.array([[1, 0], [0, -1j]], dtype=complex),
    "t": np.array([[1, 0], [0, np.exp(1j * math.pi / 4)]], dtype=complex),
    "tdg": np.array([[1, 0], [0, np.exp(-1j * math.pi / 4)]], dtype=complex),
}


def _apply_1q(psi: np.ndarray, u: np.ndarray, q: int, n: int) -> np.ndarray:
    t = psi.reshape(2 ** (n - q - 1), 2, 2**q)
    return np.einsum("ij,ajb->aib", u, t, optimize=True).reshape(-1)


def _mask_indices(n: int, ones: tuple[int, ...], zero: int) -> np.ndarray:
    idx = np.arange(1 << n, dtype=np.int64)
    sel = (idx >> zero) & 1 == 0
    for c in ones:
        sel &= (idx >> c) & 1 == 1
    return idx[sel]


def apply_gate(psi: np.ndarray, g: Gate, n: int) -> np.ndarray:
    name = g.name
    if name in _SQ:
        return _apply_1q(psi, _SQ[name], g.qubits[0], n)
    if name == "ry":
        c, s = math.cos(g.param / 2.0), math.sin(g.param / 2.0)
        u = np.array([[c, -s], [s, c]], dtype=complex)
        return _apply_1q(psi, u, g.qubits[0], n)
    if name == "rz":
        e = np.exp(-1j * g.param / 2.0)
        u = np.array([[e, 0], [0, e.conjugate()]], dtype=complex)
        return _apply_1q(psi, u, g.qubits[0], n)
    if name == "cx":
        c, t = g.qubits
        i0 = _mask_indices(n, (c,), t)
        i1 = i0 | (1 << t)
        psi = psi.copy()
        psi[i0], psi[i1] = psi[i1].copy(), psi[i0].copy()
        return psi
    if name == "cz":
        a, b = g.qubits
        idx = np.arange(1 << n, dtype=np.int64)
        sel = ((idx >> a) & 1 == 1) & ((idx >> b) & 1 == 1)
        psi = psi.copy()
        psi[sel] *= -1.0
        return psi
    if name == "swap":
        a, b = g.qubits
        idx = np.arange(1 << n, dtype=np.int64)
        i0 = idx[(((idx >> a) & 1) == 1) & (((idx >> b) & 1) == 0)]
        i1 = (i0 ^ (1 << a)) | (1 << b)
        psi = psi.copy()
        psi[i0], psi[i1] = psi[i1].copy(), psi[i0].copy()
        return psi
    if name == "mcx":
        i0 = _mask_indices(n, g.controls, g.target)
        i1 = i0 | (1 << g.target)
        psi = psi.copy()
        psi[i0], psi[i1] = psi[i1].copy(), psi[i0].copy()
        return psi
    if name == "mcry":
        i0 = _mask_indices(n, g.controls, g.target)
        i1 = i0 | (1 << g.target)
        c, s = math.cos(g.param / 2.0), math.sin(g.param / 2.0)
        psi = psi.copy()
        a0, a1 = psi[i0].copy(), psi[i1].copy()
        psi[i0] = c * a0 - s * a1
        psi[i1] = s * a0 + c * a1
        return psi
    raise ValueError(f"unknown gate {name!r}")


def simulate(circ: Circuit, psi0: np.ndarray | None = None) -> np.ndarray:
    """Run ``circ`` on |0...0> (or ``psi0``) and return the full statevector."""
    n = circ.n_qubits
    if psi0 is None:
        psi = np.zeros(1 << n, dtype=complex)
        psi[0] = 1.0
    else:
        psi = np.asarray(psi0, dtype=complex).copy()
        if psi.size != 1 << n:
            raise ValueError("psi0 has the wrong dimension")
    for g in circ.gates:
        psi = apply_gate(psi, g, n)
    return psi


def data_statevector(circ: Circuit, atol: float = 1e-10) -> np.ndarray:
    """Statevector of the data register, asserting all ancillas end in |0>.

    Ancillas are the ``n_ancilla`` highest-index wires.  Raises if they carry
    residual amplitude, which is exactly the failure mode a buggy uncomputation
    would produce.
    """
    psi = simulate(circ)
    nd, na = circ.n_data, circ.n_ancilla
    if na == 0:
        return psi
    block = psi.reshape(1 << na, 1 << nd)
    leak = float(np.linalg.norm(block[1:]))
    if leak > atol:
        raise AssertionError(f"ancillas not returned to |0>: residual norm {leak:.3e}")
    out = block[0]
    return out


def marginal_partition_is_product(psi: np.ndarray, blocks) -> tuple[bool, float]:
    """Check |psi> = tensor product across ``blocks``.

    A pure state factorises across a partition iff every block's reduced state
    is pure, so we return the worst block purity defect.
    """
    n = int(round(math.log2(psi.size)))
    worst = 0.0
    for blk in blocks:
        blk = tuple(sorted(blk))
        if len(blk) == n:
            continue
        rest = tuple(q for q in range(n) if q not in blk)
        t = psi.reshape([2] * n)
        # numpy axis j corresponds to qubit n-1-j under little-endian reshape
        axes_blk = [n - 1 - q for q in blk]
        axes_rest = [n - 1 - q for q in rest]
        mat = np.transpose(t, axes_blk + axes_rest).reshape(1 << len(blk), 1 << len(rest))
        rho = mat @ mat.conj().T
        purity = float(np.real(np.trace(rho @ rho)))
        worst = max(worst, abs(1.0 - purity))
    return worst <= 1e-9, worst


def finest_product_partition(psi: np.ndarray, declared) -> list[tuple[int, ...]]:
    """Refine each declared block into the finest genuine product partition.

    Used by the width gate: it catches a compiler that *declares* width k while
    emitting a branch that is really width k' < k (over-reporting) as well as
    one that emits a branch that does not factorise at all (under-reporting is
    caught separately by :func:`marginal_partition_is_product`).
    """
    out: list[tuple[int, ...]] = []
    for blk in declared:
        blk = tuple(sorted(blk))
        out.extend(_refine(psi, blk))
    return sorted(out)


def _refine(psi: np.ndarray, blk: tuple[int, ...]) -> list[tuple[int, ...]]:
    k = len(blk)
    if k <= 1:
        return [blk]
    # try every proper bipartition of blk (k <= ~10 in scope)
    for mask in range(1, 1 << (k - 1)):
        left = tuple(blk[i] for i in range(k) if (mask >> i) & 1)
        right = tuple(q for q in blk if q not in left)
        ok, _ = marginal_partition_is_product(psi, [left, right])
        if ok:
            return _refine(psi, left) + _refine(psi, right)
    return [blk]
