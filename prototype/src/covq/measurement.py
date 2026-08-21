"""Measurement compilation: from a probe schedule to a certified readout.

Section 11 of the manuscript is a *scoping* section.  It states local
attainability (11.1), separates the three flag operations (11.2), warns that
``4 Cov_rho`` is not the QFIM off pure states (11.3), and then says plainly in
11.5 that v0.2 "does not yet solve minimum-cost measurement synthesis".  The
paper therefore uses "realizes the QFIM" in the *state-family* sense.

This module closes the operational gap for exactly the class the compiler
emits, and no wider.  Three things are established numerically rather than
asserted:

M1  For every branch of a pair-width (indeed any cat-block) schedule, a
    **product of single-qubit equatorial measurements** attains the branch
    QFIM exactly -- not asymptotically, not up to a constant.  Because the
    label is retained, the schedule's total classical Fisher matrix is then
    ``sum_r p_r CFI_r = sum_r p_r F_r = F``.  So for this backend the
    state-family statement and the operational statement coincide, and the
    word "certified" survives the move to a real readout.

M2  The attaining measurement **depends on the operating point**, and the
    dependence is not a technicality.  At ``theta = 0`` -- the natural place
    to sit -- the fixed all-X readout returns *identically zero* Fisher
    information on every branch.  Qualification (ii) of Section 11.1 is
    therefore load-bearing: a compiler that emits a probe schedule without
    emitting the matched readout has emitted nothing.

M3  Equation (110) is implemented and checked against Equation (109).  They
    agree on pure states (110 reduces to ``4 Cov``) and disagree on noisy
    ones by a factor that is quantified here.  Reading a noisy covariance
    table as a QFIM *inflates* the reported information.

Conventions.  Generators are ``G_i = P_i / 2`` so that the pure-state QFIM is
``F_ij = <P_i P_j> - <P_i><P_j>``, matching :mod:`covq.qfim`.  Every Fisher
matrix here is per shot.

Naming.  These gates are labelled M1-M3 rather than continuing the C-series.
The manuscript's contract declares gates C1-C16; the local prototype
implements C1-C12 and the contract file that fixes the numbering has not been
received.  Renumbering into the C-series is deferred rather than guessed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .paulis import PauliSet, apply_pauli
from .qfim import qfim_from_statevector


# ----------------------------------------------------------------------
# Equatorial single-qubit readout bases
# ----------------------------------------------------------------------

def equatorial_axis(alpha: float) -> np.ndarray:
    """``M(alpha) = cos(alpha) X + sin(alpha) Y``."""
    return np.array([[0.0, np.exp(-1j * alpha)],
                     [np.exp(1j * alpha), 0.0]], dtype=complex)


def equatorial_basis_change(alpha: float) -> np.ndarray:
    """Unitary sending the ``M(alpha)`` eigenbasis to the computational basis.

    ``M(alpha)`` has eigenvectors ``(|0> +- e^{i alpha}|1>)/sqrt(2)``; the rows
    of the returned matrix are their conjugates, so applying it and then
    reading ``Z`` implements the measurement.
    """
    e = np.exp(-1j * alpha)
    return np.array([[1.0, e], [1.0, -e]], dtype=complex) / math.sqrt(2.0)


def _apply_1q(psi: np.ndarray, u: np.ndarray, q: int, n: int) -> np.ndarray:
    t = psi.reshape(2 ** (n - q - 1), 2, 2**q)
    return np.einsum("ij,ajb->aib", u, t, optimize=True).reshape(-1)


def readout_amplitudes(psi: np.ndarray, alphas) -> np.ndarray:
    """Amplitudes of ``psi`` in the product equatorial basis ``M(alpha_q)``."""
    alphas = list(alphas)
    n = len(alphas)
    if psi.size != 1 << n:
        raise ValueError(f"state has {psi.size} amplitudes, expected {1 << n}")
    out = psi.astype(complex, copy=True)
    for q, a in enumerate(alphas):
        out = _apply_1q(out, equatorial_basis_change(float(a)), q, n)
    return out


# ----------------------------------------------------------------------
# Classical Fisher information of an emitted readout
# ----------------------------------------------------------------------

def outcome_model(psi: np.ndarray, ps: PauliSet, alphas) -> tuple[np.ndarray, np.ndarray]:
    """Outcome law and its exact parameter derivatives.

    With ``|psi(theta)> = exp[-(i/2) sum_i theta_i P_i]|psi>`` and rank-one
    projectors ``Pi_x``,

        d_i p(x) = 2 Re <d_i psi| Pi_x |psi> = -Im <psi| P_i Pi_x |psi>,

    which is analytic -- no finite differences enter anywhere in this module.
    """
    c = readout_amplitudes(psi, alphas)
    p = np.abs(c) ** 2
    dp = np.empty((ps.m, p.size))
    for i in range(ps.m):
        d = readout_amplitudes(apply_pauli(psi, ps, i), alphas)
        dp[i] = -np.imag(np.conjugate(d) * c)
    return p, dp


def classical_fisher_matrix(p: np.ndarray, dp: np.ndarray,
                            tol: float = 1e-12) -> np.ndarray:
    """``CFI_ij = sum_x d_i p(x) d_j p(x) / p(x)`` over the support."""
    keep = p > tol
    w = dp[:, keep] / np.sqrt(p[keep])
    return w @ w.T


def cfi_of_readout(psi: np.ndarray, ps: PauliSet, alphas,
                   tol: float = 1e-12) -> np.ndarray:
    p, dp = outcome_model(psi, ps, alphas)
    return classical_fisher_matrix(p, dp, tol)


# ----------------------------------------------------------------------
# M1/M2: compiling the readout for one branch
# ----------------------------------------------------------------------

def compile_branch_readout(psi: np.ndarray, ps: PauliSet,
                           blocks, n_grid: int = 24) -> np.ndarray:
    """Choose per-qubit equatorial axes that attain the branch QFIM.

    For a cat block ``B`` with sign vector ``s`` the state is
    ``(e^{-i phi/2}|s> + e^{i phi/2}|-s>)/sqrt(2)`` with
    ``phi = sum_{i in B} s_i theta_i``, and the product readout gives block
    parity law ``(1 + P cos(phi - A))/4...`` with ``A`` the summed axis angle.
    Its Fisher information is *exactly* one for every ``A`` off the degenerate
    set ``sin(phi - A) = 0``, so only one angle per block is a real degree of
    freedom and a coarse scan lands on an exact optimum rather than near one.

    Blocks are independent -- the state and the readout both factorize -- so
    the scan is done blockwise, which keeps the cost linear in the number of
    blocks instead of exponential in the number of qubits.
    """
    n = int(np.log2(psi.size))
    alphas = np.zeros(n)
    grid = np.linspace(0.0, math.pi, n_grid, endpoint=False)
    for blk in blocks:
        blk = list(blk)
        idx = np.array(blk, dtype=int)
        best, best_val = 0.0, -np.inf
        for a in grid:
            trial = alphas.copy()
            trial[blk[0]] = a
            sub = cfi_of_readout(psi, ps, trial)[np.ix_(idx, idx)]
            val = float(np.trace(sub))
            if val > best_val:
                best, best_val = float(a), val
        alphas[blk[0]] = best
    return alphas


@dataclass
class ReadoutCertificate:
    """Per-branch attainability record for a labelled schedule."""

    weight: float
    blocks: list
    alphas: np.ndarray
    F_branch: np.ndarray
    CFI_branch: np.ndarray
    fixed_x_CFI: np.ndarray

    @property
    def gap(self) -> float:
        """``||F - CFI||_2``.  Zero means the readout attains the branch QFIM."""
        return float(np.linalg.norm(self.F_branch - self.CFI_branch, 2))

    @property
    def dominated(self) -> bool:
        """``CFI <= F`` in the Loewner order, as the quantum bound requires."""
        w = np.linalg.eigvalsh(self.F_branch - self.CFI_branch)
        return bool(w.min() >= -1e-9)

    @property
    def fixed_x_information(self) -> float:
        return float(np.trace(self.fixed_x_CFI))


def compile_schedule_readout(program, ps: PauliSet,
                             theta: np.ndarray | None = None) -> dict:
    """M1 + M2 for a whole labelled schedule.

    ``theta`` is the operating point; the readout is compiled *there*, which is
    the content of Section 11.1(ii).  The comparison arm is the fixed all-X
    readout, compiled nowhere.
    """
    from .sim import data_statevector

    if program.kind != "labelled":
        raise ValueError("measurement compilation is defined for labelled schedules")
    m = ps.m
    theta = np.zeros(m) if theta is None else np.asarray(theta, dtype=float)

    certs: list[ReadoutCertificate] = []
    total_F = np.zeros((m, m))
    total_CFI = np.zeros((m, m))
    total_fixed = np.zeros((m, m))
    for br in program.branches:
        psi = data_statevector(br.circuit)
        psi = _rotate(psi, ps, theta)
        F = qfim_from_statevector(psi, ps)
        alphas = compile_branch_readout(psi, ps, br.blocks)
        cfi = cfi_of_readout(psi, ps, alphas)
        fixed = cfi_of_readout(psi, ps, np.zeros(int(np.log2(psi.size))))
        certs.append(ReadoutCertificate(br.weight, list(br.blocks), alphas, F, cfi, fixed))
        total_F += br.weight * F
        total_CFI += br.weight * cfi
        total_fixed += br.weight * fixed

    return {
        "theta": theta,
        "branches": certs,
        "F_schedule": total_F,
        "CFI_schedule": total_CFI,
        "CFI_schedule_fixed_x": total_fixed,
        "max_branch_gap": max((c.gap for c in certs), default=0.0),
        "schedule_gap": float(np.linalg.norm(total_F - total_CFI, 2)),
        "all_branches_attain": all(c.gap < 1e-9 for c in certs),
        "all_branches_dominated": all(c.dominated for c in certs),
        "fixed_x_total_information": float(np.trace(total_fixed)),
        "compiled_total_information": float(np.trace(total_CFI)),
    }


def _rotate(psi: np.ndarray, ps: PauliSet, theta: np.ndarray) -> np.ndarray:
    """``exp[-(i/2) sum_i theta_i P_i] |psi>`` for commuting Pauli generators.

    Exact by Trotter-free sequential exponentiation: each factor is
    ``cos(theta_i/2) I - i sin(theta_i/2) P_i``.
    """
    out = psi.astype(complex, copy=True)
    for i, t in enumerate(theta):
        if abs(t) < 1e-15:
            continue
        out = math.cos(t / 2.0) * out - 1j * math.sin(t / 2.0) * apply_pauli(out, ps, i)
    return out


# ----------------------------------------------------------------------
# M3: Equation (110) versus Equation (109)
# ----------------------------------------------------------------------

def mixed_state_qfim(rho: np.ndarray, ps: PauliSet, tol: float = 1e-12) -> np.ndarray:
    """SLD QFIM of a mixed state, Equation (110), with ``G_i = P_i / 2``.

        (F_Q)_ij = 2 sum_{a,b: l_a + l_b > 0}
                       (l_a - l_b)^2 / (l_a + l_b) Re[<a|G_i|b><b|G_j|a>].

    On pure states this reduces to ``4 Cov`` and hence to
    :func:`~covq.qfim.qfim_from_statevector`; that reduction is checked in the
    test suite rather than taken on faith.
    """
    lam, vecs = np.linalg.eigh(np.asarray(rho, dtype=complex))
    lam = np.clip(lam.real, 0.0, None)
    m = ps.m
    # G_i in the eigenbasis.
    gs = []
    for i in range(m):
        cols = np.column_stack([apply_pauli(vecs[:, k], ps, i) for k in range(vecs.shape[1])])
        gs.append(0.5 * (vecs.conj().T @ cols))
    s = lam[:, None] + lam[None, :]
    d = lam[:, None] - lam[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(s > tol, d**2 / np.where(s > tol, s, 1.0), 0.0)
    out = np.empty((m, m))
    for i in range(m):
        for j in range(i, m):
            val = 2.0 * float(np.real(np.sum(w * gs[i] * gs[j].conj().T)))
            out[i, j] = out[j, i] = val
    return out


def covariance_surrogate(rho: np.ndarray, ps: PauliSet) -> np.ndarray:
    """Equation (109) read off a state: ``4 Cov_rho(G_i, G_j)``.

    Exact for pure states, **not** the QFIM otherwise.  Provided so the gap can
    be reported instead of silently inherited from a covariance table.
    """
    rho = np.asarray(rho, dtype=complex)
    m = ps.m
    cols = [_apply_pauli_mat(rho, ps, i) for i in range(m)]
    means = np.array([float(np.real(np.trace(c))) for c in cols])
    out = np.empty((m, m))
    for i in range(m):
        for j in range(i, m):
            second = float(np.real(np.trace(_apply_pauli_mat(cols[i], ps, j))))
            out[i, j] = out[j, i] = second - means[i] * means[j]
    return out


def _apply_pauli_mat(a: np.ndarray, ps: PauliSet, i: int) -> np.ndarray:
    """``P_i @ a`` via the statevector action, column by column."""
    return np.column_stack([apply_pauli(a[:, k], ps, i) for k in range(a.shape[1])])


def depolarize(rho: np.ndarray, p: float, qubits=None) -> np.ndarray:
    """Independent single-qubit depolarizing noise of strength ``p``."""
    n = int(np.log2(rho.shape[0]))
    qubits = range(n) if qubits is None else qubits
    out = np.asarray(rho, dtype=complex).copy()
    for q in qubits:
        acc = (1.0 - p) * out
        for pauli in ("x", "y", "z"):
            k = _single_qubit_op(pauli, q, n)
            acc = acc + (p / 3.0) * (k @ out @ k.conj().T)
        out = acc
    return out


def dephase(rho: np.ndarray, p: float, qubits=None) -> np.ndarray:
    """Independent single-qubit ``Z`` dephasing of strength ``p``."""
    n = int(np.log2(rho.shape[0]))
    qubits = range(n) if qubits is None else qubits
    out = np.asarray(rho, dtype=complex).copy()
    for q in qubits:
        k = _single_qubit_op("z", q, n)
        out = (1.0 - p) * out + p * (k @ out @ k.conj().T)
    return out


_PAULI_1Q = {
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _single_qubit_op(name: str, q: int, n: int) -> np.ndarray:
    op = _PAULI_1Q[name]
    left = np.eye(1 << (n - q - 1), dtype=complex)
    right = np.eye(1 << q, dtype=complex)
    return np.kron(np.kron(left, op), right)


def noisy_report(rho: np.ndarray, ps: PauliSet) -> dict:
    """The two matrices Section 11.3(a),(b) requires a noisy backend to report."""
    fq = mixed_state_qfim(rho, ps)
    cov = covariance_surrogate(rho, ps)
    w = np.linalg.eigvalsh(cov - fq)
    tr_fq, tr_cov = float(np.trace(fq)), float(np.trace(cov))
    return {
        "sld_qfim": fq,
        "covariance_surrogate": cov,
        "surrogate_minus_qfim_min_eig": float(w.min()),
        "surrogate_dominates": bool(w.min() >= -1e-9),
        "trace_sld_qfim": tr_fq,
        "trace_covariance_surrogate": tr_cov,
        "inflation_factor": (tr_cov / tr_fq) if tr_fq > 1e-12 else float("inf"),
        "purity": float(np.real(np.trace(np.asarray(rho) @ np.asarray(rho)))),
    }

def cfi_of_readout_mixed(rho: np.ndarray, ps: PauliSet, alphas,
                         tol: float = 1e-12) -> np.ndarray:
    """Classical Fisher matrix of the product readout on a *mixed* state.

    ``d_i p(x) = tr(Pi_x d_i rho)`` with ``d_i rho = -(i/2)[P_i, rho]`` gives
    ``d_i p(x) = Im <x| P_i rho |x>``, which agrees with the pure-state
    expression in :func:`outcome_model` and is likewise analytic.

    This is the matrix Section 11.3(b) calls "the classical Fisher matrix of a
    declared measurement", and it is the only one of the three that an
    experiment can actually claim.
    """
    rho = np.asarray(rho, dtype=complex)
    n = int(np.log2(rho.shape[0]))
    v = np.eye(rho.shape[0], dtype=complex)
    b = np.column_stack([readout_amplitudes(v[:, k], alphas) for k in range(v.shape[1])])
    rho_b = b @ rho @ b.conj().T          # rho in the readout basis
    p = np.real(np.diag(rho_b)).copy()
    dp = np.empty((ps.m, p.size))
    for i in range(ps.m):
        pr = b @ _apply_pauli_mat(rho, ps, i) @ b.conj().T
        dp[i] = np.imag(np.diag(pr))
    return classical_fisher_matrix(p, dp, tol)


def noise_chain_report(rho: np.ndarray, ps: PauliSet, alphas) -> dict:
    """The full honest chain for a noisy backend, in Loewner order.

        CFI(declared readout)  <=  SLD QFIM (110)  <=  ... and 4 Cov (109),
                                                        which bounds neither.

    The covariance surrogate is reported *outside* the chain because it is not
    an information quantity off pure states; it is included only so the size of
    the error committed by reading it as one can be seen.
    """
    fq = mixed_state_qfim(rho, ps)
    cfi = cfi_of_readout_mixed(rho, ps, alphas)
    cov = covariance_surrogate(rho, ps)
    return {
        "cfi_declared_readout": cfi,
        "sld_qfim": fq,
        "covariance_surrogate": cov,
        "cfi_within_qfim": bool(np.linalg.eigvalsh(fq - cfi).min() >= -1e-9),
        "trace_cfi": float(np.trace(cfi)),
        "trace_qfim": float(np.trace(fq)),
        "trace_surrogate": float(np.trace(cov)),
    }
