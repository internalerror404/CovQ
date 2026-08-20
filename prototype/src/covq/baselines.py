"""Baselines.

The point of this module is to make the comparison *hard* rather than
flattering.  Two of the protocol's listed baselines are weak by construction and
are labelled as such here:

``generic_full_state_preparation_then_transpile``
    Preparing a dense ``2^m`` statevector is not what a competent baseline would
    do.  Theorem A already implies that any feasible target is realised by a
    pure state with at most ``2(C(m,2)+1)`` non-zero amplitudes, so the honest
    single-state baseline is the *Caratheodory-sparse* one in
    :func:`covq.programs.single_pure_state_program`.  The dense variant is kept
    only so the gap between the two can be reported, never as the headline
    comparator.

``global_GHZ`` / ``product_probe`` / ``Hadamard_sign_schedule``
    These are fixed probes: they realise exactly one QFIM each (``J``, ``I``,
    and ``I`` respectively).  They are reference points for resource cost at a
    known target, not competitors on arbitrary targets, and reporting their
    QFIM error on a general target proves nothing.

The baseline this codebase cannot yet run is QUEST (arXiv:2605.02367).  Until it
is run, ``STATUS.md`` records the QUEST comparison as ABSENT and no claim of
advantage over it is admissible.
"""

from __future__ import annotations

import math

import numpy as np

from .circuits import Circuit
from .polytope import exact_decompose, vertex_pair_matrix, offdiag, pairs
from .programs import Program, labelled_schedule_from_signs, prepare_sparse_real, signed_cat_circuit


def global_ghz(m: int) -> Program:
    """Single GHZ probe: realises ``F = J`` (the all-ones matrix), rank one."""
    s = np.ones(m, dtype=int)
    prog = labelled_schedule_from_signs(s[None, :], np.array([1.0]))
    prog.metadata["baseline"] = "global_GHZ"
    prog.metadata["realises"] = "J (all ones)"
    return prog


def product_probe(m: int) -> Program:
    """``|+>^m``: realises ``F = I``, the identity -- zero entanglement."""
    c = Circuit(m)
    for i in range(m):
        c.h(i)
    from .programs import Branch
    prog = Program("labelled", m, [Branch(1.0, c, [(i,) for i in range(m)])])
    prog.metadata["baseline"] = "product_probe"
    prog.metadata["realises"] = "I"
    return prog


def hadamard_sign_schedule(m: int) -> Program | None:
    """Uniform schedule over Hadamard rows: realises ``F = I`` with ``m`` settings.

    Only defined when a Hadamard matrix of order ``m`` exists (here: powers of
    two).  Returns ``None`` otherwise rather than silently substituting
    something else.
    """
    if m & (m - 1) != 0:
        return None
    H = np.array([[1]])
    while H.shape[0] < m:
        H = np.block([[H, H], [H, -H]])
    S = H.astype(int)
    S = S * S[:, :1]  # quotient s ~ -s
    prog = labelled_schedule_from_signs(S, np.full(m, 1.0 / m))
    prog.metadata["baseline"] = "Hadamard_sign_schedule"
    prog.metadata["realises"] = "I"
    return prog


def random_signed_cat_schedule(F: np.ndarray, n_atoms: int,
                               rng: np.random.Generator) -> tuple[Program, float]:
    """Random settings, best convex weights: the 'no structure exploited' arm.

    Returns the program and the residual it achieves, so the value of
    *choosing* settings can be separated from the value of having many.
    """
    from scipy.optimize import linprog

    m = F.shape[0]
    S = rng.choice(np.array([-1, 1]), size=(n_atoms, m))
    S[:, 0] = np.abs(S[:, 0])
    A = vertex_pair_matrix(m, S.astype(np.int8))
    npair = A.shape[0]
    c = np.concatenate([np.zeros(n_atoms), np.ones(2 * npair)])
    A_eq = np.zeros((npair + 1, n_atoms + 2 * npair))
    A_eq[:npair, :n_atoms] = A
    A_eq[:npair, n_atoms:n_atoms + npair] = np.eye(npair)
    A_eq[:npair, n_atoms + npair:] = -np.eye(npair)
    A_eq[npair, :n_atoms] = 1.0
    b_eq = np.concatenate([offdiag(F), [1.0]])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0.0, None), method="highs")
    if not res.success:
        return labelled_schedule_from_signs(S[:1].astype(np.int8), np.array([1.0])), float("inf")
    lam = np.asarray(res.x[:n_atoms], dtype=float)
    keep = lam > 1e-12
    prog = labelled_schedule_from_signs(S[keep].astype(np.int8), lam[keep])
    prog.metadata["baseline"] = "random_signed_cat_schedule"
    recon = np.einsum("r,ri,rj->ij", lam[keep], S[keep].astype(float), S[keep].astype(float))
    np.fill_diagonal(recon, 1.0)
    return prog, float(np.linalg.norm(recon - F, "fro"))


def dense_state_preparation(signs: np.ndarray, weights: np.ndarray,
                            max_m: int = 8) -> Program:
    """Prepare the same state with a *dense* ``2^m``-amplitude routine.

    Deliberately weak; reported only to quantify how much the Caratheodory
    sparsity is worth.  See the module docstring.
    """
    m = signs.shape[1]
    if m > max_m:
        raise ValueError(f"dense state preparation capped at m={max_m}")
    probs = np.zeros(1 << m)
    for w, s in zip(weights, signs):
        b = int(np.dot(((1 - np.asarray(s, dtype=int)) // 2), 1 << np.arange(m)))
        probs[b] += float(w) / 2.0
        probs[b ^ ((1 << m) - 1)] += float(w) / 2.0
    eps = 1e-18
    probs = probs + eps
    probs = probs / probs.sum()
    amps = {k: math.sqrt(v) for k, v in enumerate(probs)}
    circ = Circuit(m)
    prepare_sparse_real(circ, range(m), amps, atol=0.0)
    return Program("single_state", m, [], circ,
                   metadata={"baseline": "generic_full_state_preparation",
                             "support": int(len(amps)),
                             "deliberately_weak": True})


def fixed_state_cat_backend(F: np.ndarray) -> Program:
    """Compile one *prescribed* state: the first atom of an exact decomposition.

    This is the 'no freedom over the equivalence class' control.  It answers the
    charter's kill criterion 5 directly: if optimising over QFIM-equivalent
    states never beats compiling one fixed state, the compiler abstraction has
    no content.
    """
    dec = exact_decompose(F)
    if not dec.feasible or dec.support == 0:
        raise ValueError("target is infeasible; no fixed-state control exists")
    s = dec.signs[int(np.argmax(dec.weights))]
    m = F.shape[0]
    from .programs import Branch
    prog = Program("labelled", m,
                   [Branch(1.0, signed_cat_circuit(s, n_qubits=m), [tuple(range(m))])])
    prog.metadata["baseline"] = "fixed_state_cat_backend"
    prog.metadata["realises"] = "s s^T for the heaviest atom only"
    return prog


BASELINES_IMPLEMENTED = [
    "global_GHZ", "product_probe", "Hadamard_sign_schedule",
    "random_signed_cat_schedule", "generic_full_state_preparation",
    "fixed_state_cat_backend", "caratheodory_single_pure_state",
]
BASELINES_ABSENT = [
    "QUEST_expectation_targeting",
    "hardware_efficient_VQA_moment_matching",
    "ADAPT_style_Pauli_rotation_targeting",
]
