"""Quantum Fisher information for commuting Pauli parameter programs.

For ``U_theta = exp[-(i/2) sum_i theta_i P_i]`` with commuting Hermitian
``P_i`` and a pure probe ``|psi>``,

    F_ij = <P_i P_j> - <P_i><P_j>      (the ordinary covariance).

The factor 4 in the pure-state QFI cancels the (1/2)^2 from the generator
normalisation, and the symmetrisation is free *because the generators commute*
(so ``P_i P_j`` is Hermitian and its expectation is real).

Normalisation declared once: every QFIM in this package is **per shot**.  A
labelled schedule that spends a fraction ``p_r`` of its shots on branch ``r``
has per-shot QFIM ``sum_r p_r F_r``; the total after ``N`` shots is ``N`` times
that.  No result mixes the two conventions.
"""

from __future__ import annotations

import numpy as np

from .paulis import PauliSet, apply_pauli


def qfim_from_statevector(psi: np.ndarray, ps: PauliSet) -> np.ndarray:
    """Exact per-shot QFIM of ``|psi>`` for the generator set ``ps``."""
    vecs = [apply_pauli(psi, ps, i) for i in range(ps.m)]
    means = np.array([float(np.real(np.vdot(psi, v))) for v in vecs])
    second = np.empty((ps.m, ps.m))
    for i in range(ps.m):
        for j in range(i, ps.m):
            val = np.vdot(vecs[i], vecs[j])
            second[i, j] = second[j, i] = float(np.real(val))
    return second - np.outer(means, means)


def pauli_means(psi: np.ndarray, ps: PauliSet) -> np.ndarray:
    return np.array([float(np.real(np.vdot(psi, apply_pauli(psi, ps, i)))) for i in range(ps.m)])


def sld_commutator_defect(psi: np.ndarray, ps: PauliSet) -> float:
    """Weak-commutativity residual ``max_ij |Im <psi| P_i P_j |psi>|``.

    The multiparameter Cramer-Rao bound is attainable for a pure-state model iff
    the mean Uhlmann curvature vanishes, i.e. ``Im <d_i psi | d_j psi> = 0``.
    Here that quantity equals ``-(1/4) Im <psi| P_i P_j |psi>``, which vanishes
    identically for commuting Hermitian generators.  This function measures the
    residual so the claim is checked rather than asserted.
    """
    vecs = [apply_pauli(psi, ps, i) for i in range(ps.m)]
    worst = 0.0
    for i in range(ps.m):
        for j in range(ps.m):
            worst = max(worst, abs(float(np.imag(np.vdot(vecs[i], vecs[j])))))
    return worst


def relative_frobenius_error(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(b, "fro")
    denom = denom if denom > 0 else 1.0
    return float(np.linalg.norm(a - b, "fro") / denom)


def operator_norm_error(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b, 2))


def estimable_subspace(F: np.ndarray, tol: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Eigen-decomposition split into estimable and gauge directions.

    Returns ``(eigenvalues, eigenvectors)`` sorted descending.  Directions with
    eigenvalue <= ``tol`` carry no information and must be quotiented out before
    any Cramer-Rao statement is made.
    """
    w, v = np.linalg.eigh(np.asarray(F, dtype=float))
    order = np.argsort(w)[::-1]
    return w[order], v[:, order]


def crb_conditioning(F: np.ndarray, tol: float = 1e-10) -> dict:
    """Report the conditioning of the CRB on the estimable subspace.

    The extreme points of the feasible set are rank-one, so a target that sits
    near a vertex is *maximally* ill-conditioned for estimation.  Any precision
    claim has to be reported alongside these numbers, not just alongside
    ``||F - F_star||``.
    """
    w, _ = estimable_subspace(F, tol)
    pos = w[w > tol]
    return {
        "rank": int(pos.size),
        "eig_max": float(w[0]) if w.size else 0.0,
        "eig_min_positive": float(pos[-1]) if pos.size else 0.0,
        "condition_number": float(w[0] / pos[-1]) if pos.size else float("inf"),
        "trace_crb_estimable": float(np.sum(1.0 / pos)) if pos.size else float("inf"),
    }


def pseudoinverse_stability(F: np.ndarray, F_star: np.ndarray,
                            tol: float = 1e-10) -> dict:
    """Conditioning-aware downstream stability, replacing submultiplicativity.

    The operational object is ``F^+`` on the identifiable quotient, not ``F``.
    If ``F_star >= gamma I`` on its support and ``||F - F_star|| <= eps < gamma``
    then

        ||F^+ - F_star^+||  <=  eps / (gamma (gamma - eps)),

    and if the two kernels differ, *no* stable inverse guarantee follows at all.
    The largest principal angle between the kernels is therefore reported
    alongside the norms: a small ``||F - F_star||`` with a rotated kernel is a
    failure, not a success.
    """
    F = np.asarray(F, dtype=float)
    F_star = np.asarray(F_star, dtype=float)
    eps = float(np.linalg.norm(F - F_star, 2))
    w_s, v_s = estimable_subspace(F_star, tol)
    w_f, v_f = estimable_subspace(F, tol)
    r_s = int((w_s > tol).sum())
    r_f = int((w_f > tol).sum())
    gamma = float(w_s[r_s - 1]) if r_s else 0.0

    ker_s = v_s[:, r_s:]
    ker_f = v_f[:, r_f:]
    if ker_s.shape[1] and ker_f.shape[1]:
        sv = np.linalg.svd(ker_s.T @ ker_f, compute_uv=False)
        max_angle = float(np.arccos(np.clip(sv.min(), -1.0, 1.0)))
    else:
        max_angle = 0.0 if ker_s.shape[1] == ker_f.shape[1] else float(np.pi / 2)

    bound = (eps / (gamma * (gamma - eps))) if (gamma > 0 and eps < gamma) else None
    return {
        "operator_norm_error": eps,
        "rank_target": r_s,
        "rank_realised": r_f,
        "ranks_agree": r_s == r_f,
        "lambda_min_positive_target": gamma,
        "condition_number_target": float(w_s[0] / gamma) if gamma > 0 else float("inf"),
        "max_principal_angle_between_kernels": max_angle,
        "pinv_error": float(np.linalg.norm(np.linalg.pinv(F, rcond=tol)
                                           - np.linalg.pinv(F_star, rcond=tol), 2)),
        "pinv_error_bound": bound,
        "bound_applies": bound is not None and r_s == r_f,
    }
