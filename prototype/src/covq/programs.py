"""Program semantics: signed-cat branches, labelled schedules, coherent flags.

Three backends are distinguished throughout, and the distinction is *measured*
rather than assumed:

labelled schedule
    Sample branch ``r`` with probability ``p_r``, run ``U_r``, and keep ``r`` in
    the classical record.  Fisher information is additive over independent runs
    with known labels, so the per-shot QFIM is exactly ``sum_r p_r F_r``.  This
    is optimal experimental design over a discrete design space; the branches
    are settings and ``p`` is the design measure.

coherent accessible flag
    Prepare ``sum_r sqrt(p_r) |r>|psi_r>`` and let the generators act only on
    data.  Because the generators are ``P (x) I``, the program QFIM depends only
    on ``rho_data``, giving

        F_Psi = sum_r p_r F_r + Cov_{r~p}(mu_r),

    with equality to the labelled schedule iff the branch means *coincide*
    (not merely vanish).  Inside the declared unit-diagonal sector every branch
    is forced zero-mean, so the two backends are indistinguishable there -- a
    degeneracy the gates below are built to expose rather than hide.

unlabelled mixed state
    Forbidden as a substitute for either.  ``F(rho) <= sum_r p_r F_r`` by
    convexity of the QFI, and the gap is real.  :func:`unlabelled_mixed_qfi`
    exists only so the gap can be *reported*.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .circuits import Circuit, lower_to_cx, resources, route
from .paulis import PauliSet, z_generators
from .qfim import qfim_from_statevector, pauli_means
from .sim import data_statevector, finest_product_partition, marginal_partition_is_product, simulate


# ----------------------------------------------------------------------
# Primitives
# ----------------------------------------------------------------------

def bits_of_signs(s: np.ndarray) -> np.ndarray:
    """``|s>`` has qubit ``i`` in state ``(1 - s_i)/2``."""
    return ((1 - np.asarray(s, dtype=int)) // 2).astype(int)


def signed_cat_circuit(s: np.ndarray, qubits=None, n_qubits: int | None = None) -> Circuit:
    """Prepare ``|C_s> = (|s> + |-s>)/sqrt(2)`` on the given qubits.

    Cost: 1 H, ``k-1`` CX, and at most ``k`` X gates for a block of size ``k``.
    """
    s = np.asarray(s, dtype=int)
    k = s.size
    qubits = list(range(k)) if qubits is None else list(qubits)
    n = k if n_qubits is None else n_qubits
    c = Circuit(n)
    if k == 1:
        c.h(qubits[0])  # |+>: zero mean, unit variance
        return c
    c.h(qubits[0])
    for i in range(1, k):
        c.cx(qubits[0], qubits[i])
    for i, b in enumerate(bits_of_signs(s)):
        if b:
            c.x(qubits[i])
    return c


def branch_circuit(blocks, s: np.ndarray, m: int) -> Circuit:
    """A width-bounded branch: one signed cat per block, tensored."""
    c = Circuit(m)
    for blk in blocks:
        blk = list(blk)
        sub = signed_cat_circuit(np.asarray(s)[blk], qubits=blk, n_qubits=m)
        c.gates.extend(sub.gates)
    return c


def prepare_sparse_real(circ: Circuit, qubits, amps: dict[int, float],
                        atol: float = 1e-14) -> None:
    """Append a preparation of ``sum_b amps[b] |b>`` on ``qubits`` (real, >= 0).

    Prefix-tree method: one multi-controlled RY per non-empty prefix, so at most
    ``d * k`` multi-controlled rotations for ``d`` non-zero amplitudes on ``k``
    qubits.  This is a *correct* generic sparse preparation, not a
    state-of-the-art one -- see ``STATUS.md``: any gate-count comparison against
    a baseline that uses it must be read as an upper bound on the baseline.
    """
    qubits = list(qubits)
    k = len(qubits)
    amps = {int(b): float(a) for b, a in amps.items() if abs(a) > atol}
    if any(a < 0 for a in amps.values()):
        raise ValueError("prepare_sparse_real expects non-negative amplitudes")
    total = sum(a * a for a in amps.values())
    if abs(total - 1.0) > 1e-8:
        raise ValueError(f"amplitudes are not normalised (norm^2 = {total})")

    def weight(prefix_bits: list[tuple[int, int]]) -> float:
        w = 0.0
        for b, a in amps.items():
            if all(((b >> pos) & 1) == val for pos, val in prefix_bits):
                w += a * a
        return w

    for level in range(k):
        pos = k - 1 - level  # bit index being decided
        prefixes: set[tuple[tuple[int, int], ...]] = set()
        for b in amps:
            prefixes.add(tuple((k - 1 - j, (b >> (k - 1 - j)) & 1) for j in range(level)))
        for prefix in sorted(prefixes):
            w_all = weight(list(prefix))
            if w_all <= atol:
                continue
            w0 = weight(list(prefix) + [(pos, 0)])
            ratio = min(1.0, max(0.0, w0 / w_all))
            theta = 2.0 * math.acos(math.sqrt(ratio))
            if abs(theta) <= 1e-15:
                continue
            ctrl_qubits = [qubits[p] for p, _ in prefix]
            zeros = [qubits[p] for p, v in prefix if v == 0]
            for q in zeros:
                circ.x(q)
            circ.mcry(ctrl_qubits, qubits[pos], theta)
            for q in zeros:
                circ.x(q)


# ----------------------------------------------------------------------
# Programs
# ----------------------------------------------------------------------

@dataclass
class Branch:
    weight: float
    circuit: Circuit
    blocks: list[tuple[int, ...]]

    @property
    def declared_width(self) -> int:
        return max((len(b) for b in self.blocks), default=1)


@dataclass
class Program:
    """A compiled program plus everything needed to audit it."""

    kind: str  # "labelled" | "coherent_flag" | "single_state"
    m: int
    branches: list[Branch] = field(default_factory=list)
    circuit: Circuit | None = None  # for single-setting backends
    n_flag: int = 0
    metadata: dict = field(default_factory=dict)

    @property
    def n_settings(self) -> int:
        return len(self.branches) if self.kind == "labelled" else 1

    def realised_qfim(self, ps: PauliSet) -> np.ndarray:
        """Per-shot QFIM, computed by simulating the emitted circuits."""
        if self.kind == "labelled":
            total = np.zeros((ps.m, ps.m))
            for br in self.branches:
                psi = data_statevector(br.circuit)
                total += br.weight * qfim_from_statevector(psi, ps)
            return total
        psi = simulate(self.circuit)
        return _qfim_on_data(psi, ps, self.circuit.n_qubits)

    def branch_means(self, ps: PauliSet) -> np.ndarray:
        out = []
        for br in self.branches:
            psi = data_statevector(br.circuit)
            out.append(pauli_means(psi, ps))
        return np.array(out) if out else np.zeros((0, ps.m))

    def emitted(self, edges=None) -> dict:
        """Lower (and optionally route) every circuit and total the resources."""
        parts = []
        circuits = [br.circuit for br in self.branches] if self.kind == "labelled" else [self.circuit]
        n_swaps = 0
        for c in circuits:
            low = lower_to_cx(c)
            if edges is not None:
                low, sw, _ = route(low, edges)
                n_swaps += sw
            parts.append(resources(low))
        agg = {
            "n_settings": len(parts),
            "two_qubit_gate_count": sum(p["two_qubit_gate_count"] for p in parts),
            "cx_count": sum(p["cx_count"] for p in parts),
            "max_two_qubit_depth": max(p["two_qubit_depth"] for p in parts),
            "sum_two_qubit_depth": sum(p["two_qubit_depth"] for p in parts),
            "max_depth": max(p["depth"] for p in parts),
            "swap_count": n_swaps,
            "n_qubits": max(p["n_qubits"] for p in parts),
            "n_ancilla": max(p["n_ancilla"] for p in parts),
        }
        agg["per_setting"] = parts
        return agg


def _qfim_on_data(psi: np.ndarray, ps: PauliSet, n_total: int) -> np.ndarray:
    """QFIM of a global pure state for generators acting on the low ``m`` wires.

    Valid because ``P (x) I`` expectations depend only on ``rho_data``.  Written
    as a sum over the ancilla index so the cost is ``2^n_anc`` statevector
    operations rather than a dense ``2^m x 2^m`` density matrix.
    """
    from .paulis import apply_pauli

    nd = ps.n
    na = n_total - nd
    if na == 0:
        return qfim_from_statevector(psi, ps)
    block = psi.reshape(1 << na, 1 << nd)
    means = np.zeros(ps.m)
    second = np.zeros((ps.m, ps.m))
    for a in range(1 << na):
        u = block[a]
        if not np.any(u):
            continue
        vecs = [apply_pauli(u, ps, i) for i in range(ps.m)]
        for i in range(ps.m):
            means[i] += float(np.real(np.vdot(u, vecs[i])))
            for j in range(i, ps.m):
                val = float(np.real(np.vdot(vecs[i], vecs[j])))
                second[i, j] += val
                if j != i:
                    second[j, i] += val
    return second - np.outer(means, means)


def labelled_schedule_from_signs(signs: np.ndarray, weights: np.ndarray) -> Program:
    """One global signed cat per atom of a ``Q_m`` decomposition."""
    m = signs.shape[1]
    branches = []
    for w, s in zip(weights, signs):
        branches.append(Branch(float(w), signed_cat_circuit(s, n_qubits=m),
                               [tuple(range(m))]))
    return Program("labelled", m, branches, metadata={"source": "global_cat"})


def labelled_schedule_from_width(dec) -> Program:
    """Labelled schedule from a :class:`~covq.width.WidthDecomposition`."""
    m = dec.m
    branches = []
    for w, blocks, s in dec.branches:
        branches.append(Branch(float(w), branch_circuit(blocks, s, m),
                               [tuple(b) for b in blocks]))
    return Program("labelled", m, branches,
                   metadata={"source": f"width_{dec.width}", "declared_width": dec.width})


def single_pure_state_program(signs: np.ndarray, weights: np.ndarray) -> Program:
    """One pure state, one setting: ``|psi> = sum_r sqrt(w_r) |C_{s_r}>``.

    This is the Caratheodory-sparse realisation implied by Theorem A.  It is the
    *strong* baseline: a feasible target never needs more than ``C(m,2)+1``
    atoms, hence at most ``2(C(m,2)+1)`` non-zero amplitudes and a single
    setting.  Comparing a CovQ schedule against dense full-state preparation
    instead of against this is a strawman.
    """
    m = signs.shape[1]
    amps: dict[int, float] = {}
    for w, s in zip(weights, signs):
        b = bits_of_signs(s)
        idx = int(np.dot(b, 1 << np.arange(m)))
        for k in (idx, idx ^ ((1 << m) - 1)):
            amps[k] = amps.get(k, 0.0) + float(w) / 2.0
    amps = {k: math.sqrt(v) for k, v in amps.items() if v > 0}
    circ = Circuit(m)
    prepare_sparse_real(circ, range(m), amps)
    return Program("single_state", m, [], circ,
                   metadata={"source": "caratheodory_single_state",
                             "support": len(amps),
                             "gate_counts_are_upper_bound": True})


def coherent_flag_program(signs: np.ndarray, weights: np.ndarray) -> Program:
    """Coherently flagged program with an accessible ``ceil(log2 q)`` flag.

    Construction: prepare ``sum_r sqrt(p_r)|r>`` on the flag, write ``b_r`` into
    the data register multiplexed on the flag, then open the cat with one H and
    ``m`` CX on a shared branch qubit.  The branch qubit is part of the retained
    flag register, which is exactly why it does not spoil the QFIM: the program
    QFIM is the covariance of ``rho_data``, and ``rho_data`` is unchanged.
    """
    m = signs.shape[1]
    q = signs.shape[0]
    nf = max(1, int(math.ceil(math.log2(q)))) if q > 1 else 1
    n = m + nf + 1
    flag = list(range(m, m + nf))
    aux = m + nf
    circ = Circuit(n, n_ancilla=nf + 1)
    amps = {r: math.sqrt(float(w)) for r, w in enumerate(weights) if w > 0}
    norm = math.sqrt(sum(a * a for a in amps.values()))
    amps = {k: v / norm for k, v in amps.items()}
    prepare_sparse_real(circ, flag, amps)
    for r, s in enumerate(signs):
        if weights[r] <= 0:
            continue
        b = bits_of_signs(s)
        zeros = [flag[p] for p in range(nf) if not ((r >> p) & 1)]
        for qz in zeros:
            circ.x(qz)
        for i in range(m):
            if b[i]:
                circ.mcx(flag, i)
        for qz in zeros:
            circ.x(qz)
    circ.h(aux)
    for i in range(m):
        circ.cx(aux, i)
    return Program("coherent_flag", m, [], circ, n_flag=nf + 1,
                   metadata={"source": "flagged_cat_superposition", "n_atoms": int(q)})


def _controlled_h(circ: Circuit, controls, t: int) -> None:
    """Multi-controlled H, using ``H = RY(pi/2) . Z`` (Z acts first)."""
    circ.h(t)
    circ.mcx(list(controls), t)
    circ.h(t)
    circ.mcry(list(controls), t, math.pi / 2.0)


def coherent_flag_program_literal(signs: np.ndarray, weights: np.ndarray) -> Program:
    """The charter's literal flagged normal form: ``sum_r sqrt(p_r)|r>|C_{s_r}>``.

    Every branch's cat preparation is controlled on the flag pattern.  Kept
    alongside :func:`coherent_flag_program` so the two can be *compared*: the
    relaxed construction realises the same QFIM because the program QFIM depends
    only on ``rho_data``, and it is markedly cheaper.  If that holds up across
    the benchmark families it is a statement about the normal form, not an
    implementation detail.
    """
    m = signs.shape[1]
    q = signs.shape[0]
    nf = max(1, int(math.ceil(math.log2(q)))) if q > 1 else 1
    n = m + nf
    flag = list(range(m, m + nf))
    circ = Circuit(n, n_ancilla=nf)
    amps = {r: math.sqrt(float(w)) for r, w in enumerate(weights) if w > 0}
    norm = math.sqrt(sum(a * a for a in amps.values()))
    amps = {k: v / norm for k, v in amps.items()}
    prepare_sparse_real(circ, flag, amps)
    for r, s in enumerate(signs):
        if weights[r] <= 0:
            continue
        zeros = [flag[p] for p in range(nf) if not ((r >> p) & 1)]
        for qz in zeros:
            circ.x(qz)
        _controlled_h(circ, flag, 0)
        for i in range(1, m):
            circ.mcx(flag + [0], i)
        for i, b in enumerate(bits_of_signs(s)):
            if b:
                circ.mcx(flag, i)
        for qz in zeros:
            circ.x(qz)
    return Program("coherent_flag", m, [], circ, n_flag=nf,
                   metadata={"source": "flagged_controlled_cat_literal", "n_atoms": int(q)})


def branch_states(program: Program) -> list[tuple[float, np.ndarray]]:
    """The ``(weight, statevector)`` list of a labelled program's branches."""
    return [(br.weight, data_statevector(br.circuit)) for br in program.branches]


def unlabelled_mixed_qfi(states: list[tuple[float, np.ndarray]],
                         ps: PauliSet) -> np.ndarray:
    """SLD QFI of the *unlabelled* mixture -- computed only to report the gap.

    Never a substitute for a program QFIM.  For ``rho = sum_k q_k |k><k|`` and
    the generator convention ``exp[-(i/2) sum theta_i P_i]``,

        F_ij = (1/2) sum_{k,l} (q_k - q_l)^2 / (q_k + q_l)
                     * Re[<k|P_i|l><l|P_j|k>],

    the prefactor being ``1/2`` rather than ``2`` because the generator is
    ``P_i / 2``.  On a pure state this reduces to ``Cov(P_i, P_j)``, which is the
    consistency check :func:`covq.gates.c12_convexity_gap` runs.
    """
    m = ps.n
    if m > 8:
        raise ValueError("unlabelled_mixed_qfi builds a dense density matrix; capped at m=8")
    dim = 1 << m
    rho = np.zeros((dim, dim), dtype=complex)
    for w, psi in states:
        rho += float(w) * np.outer(psi, psi.conj())
    vals, vecs = np.linalg.eigh(rho)
    from .paulis import apply_pauli
    cols = np.eye(dim, dtype=complex)
    P = [np.column_stack([apply_pauli(cols[:, c], ps, i) for c in range(dim)])
         for i in range(ps.m)]
    Pt = [vecs.conj().T @ p @ vecs for p in P]
    denom = vals[:, None] + vals[None, :]
    coef = np.where(denom > 1e-12,
                    (vals[:, None] - vals[None, :]) ** 2 / np.where(denom > 1e-12, denom, 1.0),
                    0.0)
    F = np.zeros((ps.m, ps.m))
    for i in range(ps.m):
        for j in range(i, ps.m):
            # Re[<k|P_i|l><l|P_j|k>] = Re[ A_i[k,l] * A_j[l,k] ]
            val = 0.5 * float(np.real(np.sum(coef * (Pt[i] * Pt[j].T))))
            F[i, j] = F[j, i] = val
    return F


# ----------------------------------------------------------------------
# Width auditing of emitted branches
# ----------------------------------------------------------------------

def audit_branch_width(program: Program) -> dict:
    """Verify declared branch widths against the *emitted* statevectors.

    Reports both directions of failure: a branch that does not factorise across
    its declared blocks (under-declared width), and a branch whose finest
    genuine product partition is strictly finer than declared (over-declared).
    """
    rows = []
    for br in program.branches:
        psi = data_statevector(br.circuit)
        ok, defect = marginal_partition_is_product(psi, br.blocks)
        finest = finest_product_partition(psi, br.blocks)
        rows.append({
            "weight": br.weight,
            "declared_width": br.declared_width,
            "declared_blocks": [list(b) for b in br.blocks],
            "factorises_as_declared": bool(ok),
            "purity_defect": defect,
            "finest_width": max((len(b) for b in finest), default=1),
            "finest_blocks": [list(b) for b in finest],
        })
    return {
        "branches": rows,
        "all_declared_widths_hold": all(r["factorises_as_declared"] for r in rows),
        "max_declared_width": max((r["declared_width"] for r in rows), default=1),
        "max_finest_width": max((r["finest_width"] for r in rows), default=1),
    }


def physical_program(program: Program, norm) -> Program:
    """Append ``U^dagger`` to every branch so the probe realises the *physical*
    generator frame: ``F(U^dag psi; P) = F(psi; U P U^dag)``."""
    inv = norm.circuit.inverse()
    if program.kind == "labelled":
        branches = [Branch(br.weight, br.circuit.copy().compose(inv), br.blocks)
                    for br in program.branches]
        out = Program("labelled", program.m, branches, metadata=dict(program.metadata))
    else:
        circ = program.circuit.copy()
        circ.compose(inv, qubit_map=range(inv.n_qubits))
        out = Program(program.kind, program.m, [], circ, program.n_flag,
                      dict(program.metadata))
    out.metadata["frame"] = "physical"
    return out


def physical_width_after_normalisation(program: Program, norm) -> dict:
    """Width of the *physical* branch states once the Clifford is appended.

    Branch width is defined in the Z frame.  For a Clifford-conjugated generator
    class the physical probe is ``U^dagger|psi_Z>`` with ``U`` generically
    globally entangling, so the physical state is entangled across all qubits
    even when the Z-frame branch is a product of pairs.  This function measures
    that gap instead of letting the paper's language paper over it.
    """
    rows = []
    for br in program.branches:
        full = br.circuit.copy()
        full.compose(norm.circuit.inverse())
        psi = data_statevector(full)
        finest = finest_product_partition(psi, [tuple(range(program.m))])
        rows.append({
            "z_frame_width": br.declared_width,
            "physical_width": max((len(b) for b in finest), default=1),
        })
    return {
        "branches": rows,
        "max_z_frame_width": max((r["z_frame_width"] for r in rows), default=1),
        "max_physical_width": max((r["physical_width"] for r in rows), default=1),
        "frame_dependent": any(r["physical_width"] > r["z_frame_width"] for r in rows),
    }


# ----------------------------------------------------------------------
# The coherent-flag loophole, made measurable
# ----------------------------------------------------------------------

def coherent_sign_purification(signs: np.ndarray, weights: np.ndarray) -> Program:
    """``|Psi_p> = sum_s sqrt(p(s)) |s>_flag |s>_data`` for a zero-mean ``p``.

    Every conditional data branch is a *product eigenstate*, so the conditional
    data width is 1 for **every feasible target**, while

        F_Psi = Cov_p(s) = F_star .

    That is the loophole: branch width is not a resource measure for a general
    coherent-flag backend.  Nothing has become cheap -- the cost has moved into
    flag-data entanglement (Schmidt rank ``2q`` across the cut), an ``m``-qubit
    flag register, and a joint flag-data readout.  Dephasing the flag collapses
    the QFIM to zero, because ``rho_data`` is diagonal in the Z basis and so
    commutes with every generator.

    Reported by :func:`flag_resource_report`; never counted as width 1.
    """
    m = signs.shape[1]
    n = 2 * m
    flag = list(range(m, 2 * m))
    amps: dict[int, float] = {}
    for w, s in zip(weights, signs):
        b = bits_of_signs(s)
        idx = int(np.dot(b, 1 << np.arange(m)))
        for k in (idx, idx ^ ((1 << m) - 1)):
            amps[k] = amps.get(k, 0.0) + float(w) / 2.0
    amps = {k: math.sqrt(v) for k, v in amps.items() if v > 0}
    circ = Circuit(n, n_ancilla=m)
    prepare_sparse_real(circ, flag, amps)
    for i in range(m):
        circ.cx(flag[i], i)
    return Program("coherent_flag", m, [], circ, n_flag=m,
                   metadata={"source": "sign_purification_loophole",
                             "conditional_data_width": 1,
                             "schmidt_rank_flag_cut": len(amps),
                             "width_is_not_a_valid_resource_here": True})


def flag_resource_report(program: Program, ps: PauliSet) -> dict:
    """Resources a coherent-flag program actually spends.

    Conditional data width alone is not one of them.  Reported instead:
    flag dimension, Schmidt rank across the flag|data cut, the QFIM with the
    flag retained coherently, and the QFIM after the flag is dephased -- the
    last being what a classical branch label would deliver.
    """
    if program.circuit is None:
        raise ValueError("flag_resource_report expects a single-setting program")
    n = program.circuit.n_qubits
    nd = program.m
    na = n - nd
    psi = simulate(program.circuit)
    block = psi.reshape(1 << na, 1 << nd)
    sv = np.linalg.svd(block, compute_uv=False)
    schmidt = int((sv > 1e-9).sum())

    F_coh = _qfim_on_data(psi, ps, n)

    # Dephasing the flag in the computational basis of the flag register turns
    # the program into a labelled schedule over the conditional data branches.
    F_deph = np.zeros((ps.m, ps.m))
    from .paulis import apply_pauli
    for a in range(1 << na):
        u = block[a]
        w = float(np.real(np.vdot(u, u)))
        if w <= 1e-15:
            continue
        u = u / math.sqrt(w)
        F_deph += w * qfim_from_statevector(u, ps)

    finest = finest_product_partition(psi, [tuple(range(n))])

    # Conditional data width: worst case over *populated* flag branches, not
    # over flag index 0, which frequently carries no weight at all.
    cond_width = 1
    n_pop = 0
    for a in range(1 << na):
        nrm = float(np.linalg.norm(block[a]))
        if nrm <= 1e-12:
            continue
        n_pop += 1
        parts = finest_product_partition(block[a] / nrm, [tuple(range(nd))])
        cond_width = max(cond_width, max((len(b) for b in parts), default=1))

    return {
        "n_flag_qubits": na,
        "flag_dimension": 1 << na,
        "schmidt_rank_flag_cut": schmidt,
        "conditional_data_width": cond_width,
        "n_populated_flag_branches": n_pop,
        "global_product_width": max((len(b) for b in finest), default=1),
        "qfim_coherent_flag": F_coh.tolist(),
        "qfim_dephased_flag": F_deph.tolist(),
        "coherence_value_trace": float(np.trace(F_coh) - np.trace(F_deph)),
        "coherence_is_load_bearing": bool(np.trace(F_coh) - np.trace(F_deph) > 1e-9),
    }
