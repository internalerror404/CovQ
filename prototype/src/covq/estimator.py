"""The estimator layer of Eq (112): from a compiled readout to a usable number.

Attaining the classical Fisher information at a point is not the same as an
estimator that reaches it at finite shots.  Eq (112) asks for

    (probe schedule, parameter circuit, branch-conditioned POVM, estimator),

and the first three are supplied by :mod:`covq.measurement`.  This module
supplies the fourth and checks it against the bound it is supposed to saturate.

The likelihood is explicit, which is the whole point of compiling to signed cat
blocks.  By the block-local probability law, a shot of branch ``r`` observed
through the matched readout yields, for each block ``B`` of that branch, a
parity ``P in {+1,-1}`` with

    Pr[P | theta] = (1 + P v_B cos(s_B . theta_B - A_B)) / 2,

independently across blocks (the state and the readout both factorise), and the
remaining ``k-1`` outcome bits per block are ancillary uniform randomness that
the estimator discards.  So the sufficient statistic of the whole campaign is
one parity count per (branch, block), and the log-likelihood is a sum of
Bernoulli terms in the block phases.

Its Fisher matrix is ``sum_r p_r sum_B s_B s_B^T = sum_r p_r F_r = F_Pi`` at
unit visibility, so a consistent efficient estimator must satisfy
``N Cov(theta_hat) -> F_Pi^+`` on the identifiable quotient.  That is the
statement tested here, by Monte Carlo, rather than assumed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize


@dataclass
class BlockChannel:
    """One (branch, block) parity channel."""

    branch: int
    block: tuple[int, ...]
    sign: np.ndarray          # s_B, entries +-1, aligned with ``block``
    analyzer: float           # A_B
    visibility: float = 1.0

    def functional(self, m: int) -> np.ndarray:
        v = np.zeros(m)
        for q, s in zip(self.block, self.sign):
            v[q] = float(s)
        return v


@dataclass
class ScheduleModel:
    """A labelled schedule together with its matched readout."""

    m: int
    weights: np.ndarray                       # p_r
    channels: list[BlockChannel] = field(default_factory=list)

    def design(self) -> np.ndarray:
        """Rows ``s_B`` stacked, one per channel."""
        return np.array([c.functional(self.m) for c in self.channels])

    def fisher(self) -> np.ndarray:
        """Per-shot Fisher matrix of the parity model."""
        F = np.zeros((self.m, self.m))
        for c in self.channels:
            u = c.functional(self.m)
            F += self.weights[c.branch] * (c.visibility ** 2) * np.outer(u, u)
        return F

    def phases(self, theta: np.ndarray) -> np.ndarray:
        return self.design() @ np.asarray(theta, dtype=float)

    def p_plus(self, theta: np.ndarray) -> np.ndarray:
        det = self.phases(theta) - np.array([c.analyzer for c in self.channels])
        vis = np.array([c.visibility for c in self.channels])
        return 0.5 * (1.0 + vis * np.cos(det))


def model_from_decomposition(dec, alphas_by_branch=None,
                             visibility: float = 1.0) -> ScheduleModel:
    """Build the estimator model from a :class:`~covq.width.WidthDecomposition`.

    With ``alphas_by_branch`` omitted the analyzers are placed at matched
    quadrature for the reference point ``theta = 0``, i.e. ``A_B = -pi/2``.
    """
    channels = []
    for r, (_, blocks, s) in enumerate(dec.branches):
        for b, blk in enumerate(blocks):
            blk = tuple(int(q) for q in blk)
            a = (-math.pi / 2.0 if alphas_by_branch is None
                 else float(alphas_by_branch[r][b]))
            channels.append(BlockChannel(r, blk, np.array([int(s[q]) for q in blk]),
                                         a, visibility))
    return ScheduleModel(dec.m, np.array([w for w, _, _ in dec.branches]), channels)


def sample_counts(model: ScheduleModel, theta: np.ndarray, n_shots: int,
                  rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Allocate shots by branch weight and draw one parity count per channel."""
    per_branch = np.floor(model.weights * n_shots).astype(int)
    per_branch[np.argmax(model.weights)] += n_shots - per_branch.sum()
    n = np.array([per_branch[c.branch] for c in model.channels])
    return rng.binomial(n, model.p_plus(theta)), n


def neg_log_likelihood(theta: np.ndarray, model: ScheduleModel,
                       k_plus: np.ndarray, n: np.ndarray) -> float:
    p = np.clip(model.p_plus(theta), 1e-15, 1.0 - 1e-15)
    return float(-(k_plus * np.log(p) + (n - k_plus) * np.log1p(-p)).sum())


def mle(model: ScheduleModel, k_plus: np.ndarray, n: np.ndarray,
        theta0: np.ndarray | None = None) -> np.ndarray:
    """Branch-conditioned maximum likelihood over the block-parity counts."""
    x0 = np.zeros(model.m) if theta0 is None else np.asarray(theta0, dtype=float)
    res = minimize(neg_log_likelihood, x0, args=(model, k_plus, n), method="BFGS")
    return res.x


def efficiency_report(model: ScheduleModel, theta_true: np.ndarray, n_shots: int,
                      n_reps: int, rng: np.random.Generator,
                      tol: float = 1e-10) -> dict:
    """Monte Carlo: does the MLE saturate the Cramer-Rao bound of the readout?

    Compares the empirical covariance of the estimator against ``F^+ / N`` on
    the identifiable quotient.  Reported as a spectral ratio rather than a
    single scalar so that an estimator that is efficient in one direction and
    wasteful in another cannot hide behind a trace.
    """
    F = model.fisher()
    w, V = np.linalg.eigh(F)
    keep = w > tol
    proj = V[:, keep]

    est = np.array([mle(model, *sample_counts(model, theta_true, n_shots, rng),
                        theta0=theta_true) for _ in range(n_reps)])
    q = est @ proj                       # estimates on the identifiable quotient
    d = int(keep.sum())
    bias = q.mean(axis=0) - proj.T @ theta_true
    emp = np.atleast_2d(np.cov(q, rowvar=False))
    crb = np.diag(1.0 / w[keep]) / n_shots

    # Symmetric whitening.  ``crb^{-1/2} emp crb^{-1/2}`` is similar to
    # ``crb^{-1} emp`` but is genuinely symmetric, so its spectrum can be taken
    # with eigvalsh; applying eigvalsh to the unsymmetric product silently reads
    # one triangle and returns the wrong numbers.
    half = np.diag(1.0 / np.sqrt(np.diag(crb)))
    ratio = np.linalg.eigvalsh(half @ emp @ half)

    # Monte Carlo reference band.  Sample-covariance eigenvalues of an exactly
    # efficient estimator already spread to the Marchenko-Pastur edges at this
    # replicate count, so a ratio inside the band is agreement, not slack.
    edge = math.sqrt(d / n_reps)
    stderr = float(np.sqrt(np.diag(emp)).max() / math.sqrt(n_reps))
    return {
        "rank": d,
        "n_shots": n_shots,
        "n_reps": n_reps,
        "max_abs_bias": float(np.abs(bias).max()),
        "bias_standard_error": stderr,
        "bias_within_4_stderr": bool(np.abs(bias).max() <= 4.0 * stderr),
        "efficiency_eigenvalues": ratio,
        "worst_efficiency_ratio": float(ratio.max()),
        "best_efficiency_ratio": float(ratio.min()),
        "mp_band": [(1.0 - edge) ** 2, (1.0 + edge) ** 2],
        "efficiency_within_mp_band": bool(ratio.min() >= (1.0 - edge) ** 2 * 0.9
                                          and ratio.max() <= (1.0 + edge) ** 2 * 1.1),
        "empirical_covariance": emp,
        "cramer_rao_bound": crb,
    }


# ----------------------------------------------------------------------
# Adaptive recentering: where the operating point actually comes from
# ----------------------------------------------------------------------

def _stage_counts(model: ScheduleModel, theta: np.ndarray, analyzers: np.ndarray,
                  n_shots: int, rng: np.random.Generator):
    shifted = ScheduleModel(model.m, model.weights,
                            [BlockChannel(c.branch, c.block, c.sign, float(a), c.visibility)
                             for c, a in zip(model.channels, analyzers)])
    return shifted, *sample_counts(shifted, theta, n_shots, rng)


def adaptive_recentering(model: ScheduleModel, theta_true: np.ndarray, n_shots: int,
                         rng: np.random.Generator, pilot_fraction: float = 0.1) -> dict:
    """Two-stage protocol that needs no oracle knowledge of the operating point.

    M2 shows the matched readout must be phased to the operating point, which
    raises the obvious objection: the operating point is what we are trying to
    estimate.  The resolution is standard two-stage adaptive estimation, and it
    costs almost nothing here because the fringe is a pure sinusoid.

    Stage 1 spends ``pilot_fraction`` of the budget split between analyzers
    ``A = 0`` and ``A = pi/2``.  Those two settings determine ``cos phi_B`` and
    ``sin phi_B``, hence ``phi_B`` itself including its sign -- one setting
    alone cannot, since ``cos`` is even.  Stage 2 re-phases every block to
    ``A_B = phi_B_hat - pi/2`` and spends the rest at unit regularity margin.

    The final estimate is a single maximum likelihood over *all* counts from
    both stages with their respective analyzers, so no pilot data is discarded.
    """
    n_pilot = max(2, int(round(pilot_fraction * n_shots)))
    half = max(1, n_pilot // 2)
    n_main = n_shots - 2 * half

    zeros = np.zeros(len(model.channels))
    quarter = np.full(len(model.channels), math.pi / 2.0)
    m_a, k_a, n_a = _stage_counts(model, theta_true, zeros, half, rng)
    m_b, k_b, n_b = _stage_counts(model, theta_true, quarter, half, rng)

    # cos and sin of each block phase from the two pilot settings.
    vis = np.array([c.visibility for c in model.channels])
    cos_hat = np.clip((2.0 * k_a / np.maximum(n_a, 1) - 1.0) / np.maximum(vis, 1e-12), -1, 1)
    sin_hat = np.clip((2.0 * k_b / np.maximum(n_b, 1) - 1.0) / np.maximum(vis, 1e-12), -1, 1)
    phi_hat = np.arctan2(sin_hat, cos_hat)

    analyzers = phi_hat - math.pi / 2.0
    m_c, k_c, n_c = _stage_counts(model, theta_true, analyzers, n_main, rng)

    def nll(theta):
        return (neg_log_likelihood(theta, m_a, k_a, n_a)
                + neg_log_likelihood(theta, m_b, k_b, n_b)
                + neg_log_likelihood(theta, m_c, k_c, n_c))

    # Seed the likelihood from the pilot itself, never from an oracle.  The
    # pilot already gives every block phase, and the phases are linear in
    # theta, so least squares on the design matrix is a consistent starting
    # estimate.  Without it the likelihood is periodic and multimodal and BFGS
    # from an arbitrary point lands in the wrong mode often enough to destroy
    # the efficiency -- which is a real property of this problem, not a
    # numerical nuisance, and is precisely why the pilot stage has to exist.
    theta_seed, *_ = np.linalg.lstsq(model.design(), phi_hat, rcond=None)
    res = minimize(nll, theta_seed, method="BFGS")
    achieved = np.abs(np.sin(model.phases(theta_true) - analyzers))
    return {
        "estimate": res.x,
        "pilot_seed": theta_seed,
        "pilot_shots": 2 * half,
        "main_shots": n_main,
        "phi_hat": phi_hat,
        "analyzers": analyzers,
        "achieved_regularity_margin": float(achieved.min()),
        "mean_regularity_margin": float(achieved.mean()),
    }


def adaptive_efficiency_report(model: ScheduleModel, theta_true: np.ndarray, n_shots: int,
                               n_reps: int, rng: np.random.Generator,
                               pilot_fraction: float = 0.1, tol: float = 1e-10) -> dict:
    """Does the two-stage protocol still reach the oracle Cramer-Rao bound?

    The reference is the *oracle* bound ``F^+/N`` -- the bound available to an
    experimenter who already knew the operating point.  Paying a pilot fraction
    to discover it must show up as an efficiency ratio above one, and the
    overhead should track the pilot fraction rather than exceed it.
    """
    F = model.fisher()
    w, V = np.linalg.eigh(F)
    keep = w > tol
    proj = V[:, keep]
    d = int(keep.sum())

    runs = [adaptive_recentering(model, theta_true, n_shots, rng, pilot_fraction)
            for _ in range(n_reps)]
    q = np.array([r["estimate"] for r in runs]) @ proj
    bias = q.mean(axis=0) - proj.T @ theta_true
    emp = np.atleast_2d(np.cov(q, rowvar=False))
    crb = np.diag(1.0 / w[keep]) / n_shots
    half = np.diag(1.0 / np.sqrt(np.diag(crb)))
    ratio = np.linalg.eigvalsh(half @ emp @ half)
    edge = math.sqrt(d / n_reps)
    return {
        "rank": d, "n_shots": n_shots, "n_reps": n_reps,
        "pilot_fraction": pilot_fraction,
        "max_abs_bias": float(np.abs(bias).max()),
        "bias_standard_error": float(np.sqrt(np.diag(emp)).max() / math.sqrt(n_reps)),
        "efficiency_eigenvalues": ratio,
        "worst_efficiency_ratio": float(ratio.max()),
        "median_regularity_margin": float(np.median([r["achieved_regularity_margin"]
                                                     for r in runs])),
        "mp_band": [(1.0 - edge) ** 2, (1.0 + edge) ** 2],
    }
