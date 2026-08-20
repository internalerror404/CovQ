"""Benchmark target families.

Every instance carries a ``ground_truth`` field that is either

* ``{"feasible": True, "hidden_decomposition": ...}`` -- built from a known
  convex combination of signed cats, so feasibility is a construction fact; or
* ``{"feasible": False, "certificate": ...}`` -- built to violate a specific
  hypermetric inequality with a *stated margin*, so infeasibility is a proof; or
* ``{"feasible": None}`` -- status genuinely unknown before the solver runs.

The third kind is deliberate.  A gate that only ever sees instances whose answer
was baked in is a gate that cannot fail.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

PROTOCOL_SEEDS = [2026, 3407, 9181, 17041, 27183]


@dataclass
class Instance:
    name: str
    family: str
    m: int
    F: np.ndarray
    ground_truth: dict = field(default_factory=lambda: {"feasible": None})
    meta: dict = field(default_factory=dict)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def psd_shrink(F: np.ndarray, atol: float = 1e-10) -> tuple[np.ndarray, float]:
    """Mix ``F`` toward ``I`` by the smallest ``alpha`` that restores PSD.

    Returns ``(F_shrunk, alpha)``; ``alpha = 0`` means no shrinkage was needed.
    The shrinkage factor is always reported, never silently applied.
    """
    F = np.asarray(F, dtype=float).copy()
    np.fill_diagonal(F, 1.0)
    w = np.linalg.eigvalsh(F)
    if w.min() >= -atol:
        return F, 0.0
    lo = float(w.min())
    alpha = -lo / (1.0 - lo) + 1e-9
    G = (1.0 - alpha) * F + alpha * np.eye(F.shape[0])
    np.fill_diagonal(G, 1.0)
    return G, float(alpha)


def normalise_to_unit_diagonal(F: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Declared normalisation ``F -> D^{-1/2} F D^{-1/2}``, ``D = diag(F)``.

    Returns the normalised matrix and the scale vector ``sqrt(diag(F))`` so the
    transformation can be undone and so no downstream number is reported in a
    convention that has not been stated.
    """
    F = np.asarray(F, dtype=float)
    d = np.sqrt(np.clip(np.diag(F), 1e-300, None))
    return F / np.outer(d, d), d


def from_cat_mixture(signs: np.ndarray, weights: np.ndarray) -> np.ndarray:
    S = np.asarray(signs, dtype=float)
    w = np.asarray(weights, dtype=float)
    w = w / w.sum()
    F = np.einsum("r,ri,rj->ij", w, S, S)
    np.fill_diagonal(F, 1.0)
    return F


# ----------------------------------------------------------------------
# Feasible by construction
# ----------------------------------------------------------------------

def random_signed_cat_mixture(m: int, support: int, rng: np.random.Generator) -> Instance:
    S = rng.choice(np.array([-1, 1]), size=(support, m))
    S[:, 0] = np.abs(S[:, 0])  # quotient s ~ -s
    w = rng.dirichlet(np.ones(support))
    F = from_cat_mixture(S, w)
    return Instance(f"cat_mix_m{m}_q{support}", "feasible_hidden_decomposition", m, F,
                    {"feasible": True,
                     "hidden_decomposition": {"signs": S.tolist(), "weights": w.tolist()}},
                    {"support": int(support)})


def signed_community(m: int, n_comm: int, rng: np.random.Generator,
                     support: int = 8) -> Instance:
    """Cat mixture whose sign vectors respect a community partition."""
    comm = rng.integers(0, n_comm, size=m)
    S = np.ones((support, m), dtype=int)
    for r in range(support):
        flip = rng.choice(np.array([-1, 1]), size=n_comm)
        S[r] = flip[comm]
        if S[r, 0] < 0:
            S[r] *= -1
    w = rng.dirichlet(np.ones(support))
    F = from_cat_mixture(S, w)
    return Instance(f"community_m{m}_c{n_comm}", "feasible_hidden_decomposition", m, F,
                    {"feasible": True,
                     "hidden_decomposition": {"signs": S.tolist(), "weights": w.tolist()}},
                    {"communities": comm.tolist()})


def block_diagonal(m: int, block: int, rng: np.random.Generator) -> Instance:
    """Independent cat mixtures on disjoint blocks -- feasible, width <= block."""
    blocks = [tuple(range(i, min(i + block, m))) for i in range(0, m, block)]
    F = np.eye(m)
    per_block = {}
    for blk in blocks:
        k = len(blk)
        if k == 1:
            continue
        S = rng.choice(np.array([-1, 1]), size=(min(4, 1 << (k - 1)), k))
        S[:, 0] = np.abs(S[:, 0])
        w = rng.dirichlet(np.ones(S.shape[0]))
        Fb = from_cat_mixture(S, w)
        for a, i in enumerate(blk):
            for b, j in enumerate(blk):
                F[i, j] = Fb[a, b]
        per_block[str(blk)] = {"signs": S.tolist(), "weights": w.tolist()}
    np.fill_diagonal(F, 1.0)
    return Instance(f"blockdiag_m{m}_b{block}", "feasible_hidden_decomposition", m, F,
                    {"feasible": True, "hidden_decomposition": {"per_block": per_block}},
                    {"blocks": [list(b) for b in blocks], "max_block": block})


def matching_target(m: int, rng: np.random.Generator, load: float = 0.9) -> Instance:
    """Signed weight on a single random matching: width 2 by construction."""
    if abs(load) > 1.0:
        raise ValueError("a single matching edge needs |F_ij| <= 1 to stay PSD")
    perm = rng.permutation(m)
    F = np.eye(m)
    edges = []
    for a in range(0, m - 1, 2):
        i, j = int(perm[a]), int(perm[a + 1])
        val = float(load) * (1.0 if rng.random() < 0.5 else -1.0)
        F[i, j] = F[j, i] = val
        edges.append((i, j))
    return Instance(f"matching_m{m}_load{load:g}", "width2_probe", m, F,
                    {"feasible": True, "expected_min_width": 2 if edges else 1},
                    {"load": float(load), "edges": edges})


def path_target(m: int, c: float) -> Instance:
    """Nearest-neighbour chain ``F_{i,i+1} = c``.

    Sharp width-2 probe: interior degree sums are ``2c``, so the target is in
    ``Q^prog_{m,2}`` iff ``|c| <= 1/2``, while PSD survives up to
    ``1 / (2 cos(pi/(m+1)))``, which is strictly larger.  The window between the
    two is where a target is a perfectly good QFIM that *provably* needs a
    width-3 branch.
    """
    F = np.eye(m)
    for i in range(m - 1):
        F[i, i + 1] = F[i + 1, i] = c
    psd_limit = 1.0 / (2.0 * np.cos(np.pi / (m + 1)))
    return Instance(f"path_m{m}_c{c:g}", "width2_probe", m, F,
                    {"feasible": None,
                     "expected_min_width": 2 if abs(c) <= 0.5 + 1e-12 else None},
                    {"c": float(c), "width2_threshold": 0.5,
                     "psd_threshold": float(psd_limit),
                     "is_psd": bool(abs(c) < psd_limit)})


def star_target(m: int, c: float) -> Instance:
    """Star ``F_{0,i} = c``: width 2 iff ``(m-1)|c| <= 1``, PSD iff ``|c| <= 1/sqrt(m-1)``."""
    F = np.eye(m)
    for i in range(1, m):
        F[0, i] = F[i, 0] = c
    return Instance(f"star_m{m}_c{c:g}", "width2_probe", m, F,
                    {"feasible": None,
                     "expected_min_width": 2 if (m - 1) * abs(c) <= 1 + 1e-12 else None},
                    {"c": float(c), "width2_threshold": 1.0 / (m - 1),
                     "psd_threshold": float(1.0 / np.sqrt(m - 1))})


# ----------------------------------------------------------------------
# Status unknown before solving
# ----------------------------------------------------------------------

def toeplitz_like(m: int, rho: float) -> Instance:
    idx = np.arange(m)
    F = rho ** np.abs(idx[:, None] - idx[None, :])
    np.fill_diagonal(F, 1.0)
    return Instance(f"toeplitz_m{m}_rho{rho:g}", "structured_unknown", m, F,
                    {"feasible": None}, {"rho": float(rho)})


def banded(m: int, bandwidth: int, rho: float) -> Instance:
    idx = np.arange(m)
    d = np.abs(idx[:, None] - idx[None, :])
    F = np.where(d <= bandwidth, rho ** d, 0.0)
    np.fill_diagonal(F, 1.0)
    F, alpha = psd_shrink(F)
    return Instance(f"banded_m{m}_bw{bandwidth}", "structured_unknown", m, F,
                    {"feasible": None}, {"bandwidth": bandwidth, "rho": float(rho),
                                         "psd_shrinkage": alpha})


def hierarchical(m: int, rng: np.random.Generator, levels: int = 2) -> Instance:
    F = np.eye(m)
    idx = np.arange(m)
    for lvl in range(levels):
        size = max(2, m >> (lvl + 1))
        val = 0.6 / (lvl + 1)
        grp = idx // size
        F = F + val * (grp[:, None] == grp[None, :]).astype(float)
    np.fill_diagonal(F, 0.0)
    scale = np.abs(F).max()
    if scale > 0:
        F = F / scale * 0.8
    np.fill_diagonal(F, 1.0)
    F, alpha = psd_shrink(F)
    return Instance(f"hierarchical_m{m}_l{levels}", "structured_unknown", m, F,
                    {"feasible": None}, {"levels": levels, "psd_shrinkage": alpha})


def low_effective_rank(m: int, rank: int, rng: np.random.Generator) -> Instance:
    G = rng.standard_normal((rank, m))
    G = G / np.linalg.norm(G, axis=0, keepdims=True)
    F = G.T @ G
    np.fill_diagonal(F, 1.0)
    F, alpha = psd_shrink(F)
    return Instance(f"lowrank_m{m}_r{rank}", "structured_unknown", m, F,
                    {"feasible": None}, {"rank": rank, "psd_shrinkage": alpha})


# ----------------------------------------------------------------------
# Infeasible with a proof
# ----------------------------------------------------------------------

def hypermetric_violator(m: int, k: int, delta: float,
                         rng: np.random.Generator | None = None) -> Instance:
    """PSD, unit-diagonal, and provably outside ``Q_m`` with margin ``delta``.

    Take ``b`` the indicator of an odd ``k``-subset ``T`` and set
    ``F_ij = -1/k - delta/(k(k-1))`` inside ``T``.  Then

        b^T F b = 1 - delta   <   1,

    while the ``T`` block has eigenvalues ``(1 - delta)/k`` and ``1 - c > 0``,
    so ``F`` is positive semidefinite for ``delta`` in ``(0, 1]``.  For
    ``m = 3, k = 3, delta = 0.4`` this is ``F_ij = -0.4`` with spectrum
    ``{0.2, 1.4, 1.4}``: positive definite, unit diagonal, and infeasible.
    """
    if k % 2 == 0 or k < 3:
        raise ValueError("k must be odd and at least 3")
    if not 0.0 < delta <= 1.0:
        raise ValueError("delta must lie in (0, 1]")
    rng = np.random.default_rng(0) if rng is None else rng
    T = sorted(rng.choice(m, size=k, replace=False).tolist()) if m > k else list(range(k))
    c = -1.0 / k - delta / (k * (k - 1))
    F = np.eye(m)
    for a in T:
        for b in T:
            if a != b:
                F[a, b] = c
    bvec = np.zeros(m, dtype=int)
    bvec[T] = 1
    return Instance(f"violator_m{m}_k{k}_d{delta:g}", "infeasible_certified", m, F,
                    {"feasible": False,
                     "certificate": {"kind": "hypermetric", "b": bvec.tolist()},
                     "certified_margin": float(delta)},
                    {"subset": T, "offdiag_value": float(c)})


def near_facet_family(m: int, k: int = 3,
                      deltas=(0.01, 0.05, 0.2, 0.6),
                      rng: np.random.Generator | None = None) -> list[Instance]:
    """Targets approaching a facet from outside, with exactly known margins."""
    return [hypermetric_violator(m, k, d, rng) for d in deltas]


def on_facet(m: int, k: int = 3) -> Instance:
    """The boundary case ``b^T F b = 1`` exactly: feasible, but with no slack."""
    T = list(range(k))
    c = -1.0 / k
    F = np.eye(m)
    for a in T:
        for b in T:
            if a != b:
                F[a, b] = c
    return Instance(f"on_facet_m{m}_k{k}", "boundary", m, F, {"feasible": None},
                    {"subset": T, "offdiag_value": float(c),
                     "note": "hypermetric inequality is tight"})


# ----------------------------------------------------------------------
# Suites
# ----------------------------------------------------------------------

def task0_suite(m_values=(3, 4, 5, 6), seed: int = PROTOCOL_SEEDS[0]) -> list[Instance]:
    """The small-``m`` suite used by the exact prototype and the gates."""
    rng = np.random.default_rng(seed)
    out: list[Instance] = []
    for m in m_values:
        out.append(random_signed_cat_mixture(m, min(4, 1 << (m - 1)), rng))
        out.append(random_signed_cat_mixture(m, min(2, 1 << (m - 1)), rng))
        if m >= 4:
            out.append(block_diagonal(m, 2, rng))
            out.append(signed_community(m, 2, rng))
        out.append(matching_target(m, rng, load=0.9))
        out.append(path_target(m, 0.3))
        if m >= 4:
            psd_lim = 1.0 / (2.0 * np.cos(np.pi / (m + 1)))
            out.append(path_target(m, round(0.5 * (0.5 + psd_lim), 3)))
            out.append(star_target(m, 0.9 / np.sqrt(m - 1)))
        out.append(toeplitz_like(m, 0.5))
        out.append(low_effective_rank(m, max(1, m // 2), rng))
        out.append(hypermetric_violator(m, 3, 0.4, rng))
        out.append(on_facet(m, 3))
        if m >= 5:
            out.extend(near_facet_family(m, 3, (0.02, 0.2), rng))
            out.append(hypermetric_violator(m, 5, 0.3, rng))
    return out
