"""J1: independent regression on the mixed-state QFIM path.

The spectral implementation of Eq (110) had a real defect -- a
transpose-versus-conjugate-transpose error that was invisible on every
real-valued fixture and produced a non-PSD matrix on complex ones.  A second
implementation that shares none of its machinery is therefore the right check.

This file solves the defining Lyapunov equation directly,

    L_i rho + rho L_i = 2 d_i rho,      d_i rho = -i [G_i, rho],  G_i = P_i / 2,

with a general Sylvester solver, and forms

    F_ij = (1/2) Tr[ rho (L_i L_j + L_j L_i) ].

No eigendecomposition of rho appears on this path, so agreement with the
spectral formula is evidence rather than restatement.

Scope of the tolerance.  The Lyapunov equation has a unique solution exactly
when rho is full rank; on a rank-deficient state any operator supported on the
kernel may be added to L_i, so the SLD is *not unique* and an equality check to
a fixed tolerance is not the right question there.  The registered
``spectral_vs_sld_max_abs = 1e-10`` is therefore asserted on full-rank states,
where the solve is well posed.  Rank-deficient states are checked as the
regularised limit instead, and are required to converge at the predicted linear
rate -- which is what confirms the spectral formula computes the correct
support-restricted SLD rather than merely something nearby.
"""

from __future__ import annotations

import numpy as np
import pytest

from covq.measurement import mixed_state_qfim, sld_qfim
from covq.paulis import random_commuting_paulis, z_generators
from covq.qfim import qfim_from_statevector

SPECTRAL_VS_SLD_MAX_ABS = 1.0e-10
PSD_FLOOR = -1.0e-10


def _random_full_rank(dim, rng):
    a = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    rho = a @ a.conj().T
    return rho / np.trace(rho).real


def test_j1_sld_solve_matches_the_spectral_formula_on_full_rank_states():
    """The registered acceptance criterion, on the states where it is well posed."""
    rng = np.random.default_rng(1234)
    worst = 0.0
    checked = 0
    for m, n in ((1, 1), (2, 2), (3, 3), (2, 3), (3, 4)):
        ps = z_generators(m) if n == m else random_commuting_paulis(m, rng, n=n)
        for _ in range(6):
            rho = _random_full_rank(1 << n, rng)
            assert np.linalg.eigvalsh(rho).min() > 1e-6, "fixture must be full rank"
            worst = max(worst, float(np.abs(sld_qfim(rho, ps)
                                            - mixed_state_qfim(rho, ps)).max()))
            checked += 1
    assert checked == 30
    assert worst <= SPECTRAL_VS_SLD_MAX_ABS, worst


def test_j1_both_paths_are_positive_semidefinite():
    rng = np.random.default_rng(4321)
    for m, n in ((2, 2), (3, 3)):
        ps = z_generators(m) if n == m else random_commuting_paulis(m, rng, n=n)
        for _ in range(8):
            rho = _random_full_rank(1 << n, rng)
            for F in (mixed_state_qfim(rho, ps), sld_qfim(rho, ps)):
                assert np.linalg.eigvalsh(F).min() >= PSD_FLOOR


def test_j1_rank_deficient_states_converge_at_the_predicted_linear_rate():
    """Singular by construction, so checked as a limit rather than an equality.

    ``rho + eps I`` restores uniqueness and the residual against the spectral
    formula must fall linearly in ``eps``.  Linear convergence is the
    substantive claim: it identifies the spectral formula as the ``eps -> 0``
    limit of the well-posed problem, which is what makes it the correct
    support-restricted SLD.
    """
    rng = np.random.default_rng(99)
    for m in (2, 3):
        ps = z_generators(m)
        dim = 1 << m
        v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
        v /= np.linalg.norm(v)
        pure = np.outer(v, v.conj())
        b = rng.normal(size=(dim, 2)) + 1j * rng.normal(size=(dim, 2))
        rank2 = b @ b.conj().T
        rank2 /= np.trace(rank2).real

        for rho, ref in ((pure, qfim_from_statevector(v, ps)),
                         (rank2, mixed_state_qfim(rank2, ps))):
            assert np.linalg.matrix_rank(rho, tol=1e-10) < dim
            coarse = float(np.abs(sld_qfim(rho, ps, ridge=1e-5) - ref).max())
            fine = float(np.abs(sld_qfim(rho, ps, ridge=1e-7) - ref).max())
            assert fine < coarse
            # Linear in eps: two decades of eps buy two decades of residual.
            assert 30.0 < coarse / fine < 300.0, (coarse, fine)
            assert fine < 1e-5


def test_j1_pure_state_limit_recovers_the_covariance_formula():
    """Closing the loop: SLD solve -> spectral Eq (110) -> pure-state covariance."""
    rng = np.random.default_rng(7)
    for m in (2, 3):
        ps = z_generators(m)
        v = rng.normal(size=1 << m) + 1j * rng.normal(size=1 << m)
        v /= np.linalg.norm(v)
        rho = np.outer(v, v.conj())
        ref = qfim_from_statevector(v, ps)
        assert mixed_state_qfim(rho, ps) == pytest.approx(ref, abs=1e-10)
        extrapolated = 2.0 * sld_qfim(rho, ps, ridge=1e-6) - sld_qfim(rho, ps, ridge=2e-6)
        assert extrapolated == pytest.approx(ref, abs=1e-8)
