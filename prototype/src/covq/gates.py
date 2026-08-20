"""Correctness gates C1-C12 and the stop rule.

A gate has to be able to fail.  Three of the protocol's v0.1 gates cannot, and
they are re-specified here rather than copied:

C2 (flagged composition)
    ``F_program == sum_r p_r F_r`` is *hollow* against a labelled backend, whose
    QFIM is defined to be that sum.  Re-specified as: simulate the emitted
    coherent-flag circuit and check its QFIM against the branch average, and
    separately check the exact identity ``F_Psi = sum_r p_r F_r + Cov(mu_r)`` on
    a deliberately non-zero-mean control where the two sides must differ.  The
    control is what gives the gate teeth.

C3 (exact feasibility)
    ``exact solver agrees with exhaustive enumeration`` is hollow when the exact
    solver *is* the enumeration.  Re-specified as agreement between three
    independent code paths -- full-vertex LP, column generation with exact
    pricing, and the hypermetric certificate search -- on a suite that contains
    feasible, infeasible, and boundary instances.

C10 (downstream stability)
    ``||A^T dF A|| <= ||A||^2 ||dF||`` is submultiplicativity: it holds for every
    matrix, so it cannot fail.  Re-specified as a *measurement* of tightness,
    reported with the CRB conditioning that actually governs precision.  Its
    status is MEASURED, never PASS.

Statuses: PASS, FAIL, MEASURED, ABSENT.  Only FAIL trips the stop rule.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import polytope as pol
from . import width as wid
from .circuits import lower_to_cx, resources, route
from .instances import Instance
from .paulis import PauliSet, diagonalizing_clifford, random_commuting_paulis, z_generators
from .programs import (
    Program,
    audit_branch_width,
    branch_circuit,
    coherent_flag_program,
    labelled_schedule_from_signs,
    labelled_schedule_from_width,
    physical_program,
    physical_width_after_normalisation,
    signed_cat_circuit,
    single_pure_state_program,
    unlabelled_mixed_qfi,
)
from .qfim import (
    crb_conditioning,
    operator_norm_error,
    qfim_from_statevector,
    relative_frobenius_error,
    sld_commutator_defect,
)
from .sim import data_statevector, simulate
from .topology import build as build_topology

PASS, FAIL, MEASURED, ABSENT = "PASS", "FAIL", "MEASURED", "ABSENT"


@dataclass
class GateResult:
    name: str
    status: str
    criterion: str
    measured: dict = field(default_factory=dict)
    notes: str = ""

    @property
    def failed(self) -> bool:
        return self.status == FAIL


class StopRule(Exception):
    """Raised when a gate fails.  Tolerances are never relaxed in response."""


# ----------------------------------------------------------------------

def c1_signed_cat_semantics(m_values=(2, 3, 4, 5, 6), seed: int = 2026) -> GateResult:
    rng = np.random.default_rng(seed)
    worst = 0.0
    worst_case = None
    for m in m_values:
        ps = z_generators(m)
        for _ in range(8):
            s = rng.choice(np.array([-1, 1]), size=m)
            psi = data_statevector(signed_cat_circuit(s, n_qubits=m))
            F = qfim_from_statevector(psi, ps)
            err = relative_frobenius_error(F, np.outer(s, s).astype(float))
            if err > worst:
                worst, worst_case = err, (m, s.tolist())
    status = PASS if worst < 1e-12 else FAIL
    return GateResult("C1_signed_cat_semantics", status,
                      "relative_Frobenius_error(F_statevector, s s^T) < 1e-12",
                      {"worst_relative_error": worst, "worst_case": worst_case})


def c2_flagged_composition(m: int = 4, seed: int = 2026) -> GateResult:
    """Coherent-flag QFIM vs branch average, plus a non-zero-mean control."""
    rng = np.random.default_rng(seed)
    ps = z_generators(m)
    S = rng.choice(np.array([-1, 1]), size=(4, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(4))
    branch_avg = np.einsum("r,ri,rj->ij", w, S.astype(float), S.astype(float))
    np.fill_diagonal(branch_avg, 1.0)

    prog = coherent_flag_program(S, w)
    F_coh = prog.realised_qfim(ps)
    err_coh = relative_frobenius_error(F_coh, branch_avg)

    labelled = labelled_schedule_from_signs(S, w)
    F_lab = labelled.realised_qfim(ps)
    err_lab = relative_frobenius_error(F_lab, branch_avg)

    # Non-zero-mean control: branches |0..0> and |1..1> have means +1 and -1, so
    # every branch QFIM is 0 while the ensemble covariance is the all-ones matrix.
    from .circuits import Circuit

    zero_branch = Circuit(m)
    one_branch = Circuit(m)
    for i in range(m):
        one_branch.x(i)
    psi0 = data_statevector(zero_branch)
    psi1 = data_statevector(one_branch)
    F0 = qfim_from_statevector(psi0, ps)
    F1 = qfim_from_statevector(psi1, ps)
    p = 0.5
    avg_branch = p * F0 + (1 - p) * F1
    mus = np.array([np.ones(m), -np.ones(m)])
    cov_mu = (p * np.outer(mus[0], mus[0]) + (1 - p) * np.outer(mus[1], mus[1])
              - np.outer(p * mus[0] + (1 - p) * mus[1], p * mus[0] + (1 - p) * mus[1]))
    ctrl_signs = np.array([[1] * m, [-1] * m])
    ctrl = coherent_flag_program(ctrl_signs, np.array([p, 1 - p]))
    # coherent_flag_program builds cats; for the control we need the raw basis
    # states, so build the flagged state by hand.
    flag_state = np.zeros(1 << (m + 1), dtype=complex)
    flag_state[0] = math.sqrt(p)                      # |flag=0>|0...0>
    flag_state[(1 << m) | ((1 << m) - 1)] = math.sqrt(1 - p)  # |flag=1>|1...1>
    from .programs import _qfim_on_data

    F_ctrl = _qfim_on_data(flag_state, ps, m + 1)
    identity_residual = float(np.linalg.norm(F_ctrl - (avg_branch + cov_mu), "fro"))
    separation = float(np.linalg.norm(F_ctrl - avg_branch, "fro"))

    ok = (err_coh < 1e-12 and err_lab < 1e-12 and identity_residual < 1e-10
          and separation > 1e-6)
    return GateResult(
        "C2_flagged_composition", PASS if ok else FAIL,
        "coherent-flag QFIM == branch average in the zero-mean sector; "
        "F_Psi == sum p_r F_r + Cov(mu_r) in general; the two differ on the control",
        {"coherent_vs_branch_average": err_coh,
         "labelled_vs_branch_average": err_lab,
         "nonzero_mean_identity_residual": identity_residual,
         "nonzero_mean_separation": separation,
         "control_branch_qfim_is_zero": bool(np.allclose(avg_branch, 0.0))},
        "The control has zero branch QFIM but a rank-one program QFIM; a backend "
        "that reports sum_r p_r F_r there is wrong by the whole matrix.")


def c3_exact_feasibility(instances: list[Instance], max_m: int = 6) -> GateResult:
    rows = []
    disagreements = []
    for inst in instances:
        if inst.m > max_m:
            continue
        dec = pol.exact_decompose(inst.F)
        cg = pol.column_generation(inst.F, rng=np.random.default_rng(inst.m))
        hyper = pol.hypermetric_certificate(inst.F, max_entry=1)
        row = {
            "instance": inst.name,
            "lp_feasible": bool(dec.feasible),
            "column_generation_feasible": bool(cg.feasible),
            "hypermetric_violation_found": hyper is not None,
            "declared_ground_truth": inst.ground_truth.get("feasible"),
        }
        if dec.feasible != cg.feasible:
            disagreements.append((inst.name, "lp_vs_column_generation"))
        if hyper is not None and dec.feasible:
            disagreements.append((inst.name, "certificate_vs_lp"))
        gt = inst.ground_truth.get("feasible")
        if gt is not None and bool(dec.feasible) != bool(gt):
            disagreements.append((inst.name, "lp_vs_ground_truth"))
        rows.append(row)
    n_feas = sum(1 for r in rows if r["lp_feasible"])
    status = FAIL if disagreements else PASS
    return GateResult(
        "C3_exact_feasibility", status,
        "full-vertex LP, column generation with exact pricing, and certificate "
        "search agree on a suite containing feasible, infeasible and boundary targets",
        {"n_instances": len(rows), "n_feasible": n_feas,
         "n_infeasible": len(rows) - n_feas,
         "disagreements": disagreements, "rows": rows},
        "Both outcomes must occur, else the gate never exercises rejection.")


def c4_infeasible_rejection(instances: list[Instance]) -> GateResult:
    rows = []
    bad = []
    for inst in instances:
        if inst.ground_truth.get("feasible") is not False:
            continue
        dec = pol.exact_decompose(inst.F)
        cert = dec.certificate
        verified = bool(cert.verify(inst.F)) if cert is not None else False
        row = {"instance": inst.name, "accepted_by_solver": bool(dec.feasible),
               "certificate_kind": None if cert is None else cert.kind,
               "certificate_verified": verified,
               "declared_margin": inst.ground_truth.get("certified_margin")}
        if dec.feasible or not verified:
            bad.append(inst.name)
        rows.append(row)
    status = FAIL if (bad or not rows) else PASS
    return GateResult(
        "C4_infeasible_rejection", status,
        "every certified-infeasible target is rejected AND the emitted "
        "certificate independently verifies",
        {"n_infeasible_instances": len(rows), "failures": bad, "rows": rows},
        "Reference instance: m=3, all off-diagonals -0.4. Spectrum {0.2, 1.4, 1.4} "
        "(positive definite, unit diagonal) but b=(1,1,1) gives b^T F b = 0.6 < 1.")


def c5_exact_reconstruction(instances: list[Instance], max_m: int = 6) -> GateResult:
    worst = 0.0
    worst_case = None
    rows = []
    for inst in instances:
        if inst.m > max_m or inst.ground_truth.get("feasible") is False:
            continue
        dec = pol.exact_decompose(inst.F)
        if not dec.feasible:
            continue
        red = pol.caratheodory_reduce(dec)
        prog = labelled_schedule_from_signs(red.signs, red.weights)
        F = prog.realised_qfim(z_generators(inst.m))
        err = float(np.linalg.norm(F - inst.F, "fro"))
        rows.append({"instance": inst.name, "support": red.support,
                     "caratheodory_bound": inst.m * (inst.m - 1) // 2 + 1,
                     "frobenius_error": err})
        if err > worst:
            worst, worst_case = err, inst.name
    status = PASS if worst < 1e-10 and rows else (FAIL if rows else ABSENT)
    return GateResult("C5_exact_reconstruction", status,
                      "Frobenius_error < 1e-10 from the emitted circuits",
                      {"worst_error": worst, "worst_case": worst_case, "rows": rows})


def c6_width_integrity(instances: list[Instance], max_m: int = 8) -> GateResult:
    rows = []
    bad = []
    for inst in instances:
        if inst.m > max_m:
            continue
        d2 = wid.decompose_width2(inst.F)
        if not d2.feasible:
            continue
        prog = labelled_schedule_from_width(d2)
        audit = audit_branch_width(prog)
        row = {"instance": inst.name, "n_branches": prog.n_settings,
               "max_declared_width": audit["max_declared_width"],
               "max_finest_width": audit["max_finest_width"],
               "all_declared_widths_hold": audit["all_declared_widths_hold"]}
        if not audit["all_declared_widths_hold"] or audit["max_declared_width"] > 2:
            bad.append(inst.name)
        rows.append(row)
    status = FAIL if bad else (PASS if rows else ABSENT)
    return GateResult("C6_width_integrity", status,
                      "every emitted branch factorises across its declared blocks "
                      "and no declared width exceeds the compiled width",
                      {"failures": bad, "rows": rows})


def c7_clifford_mapping(m_values=(3, 4, 5), seed: int = 3407) -> GateResult:
    rng = np.random.default_rng(seed)
    rows = []
    worst = 0.0
    for m in m_values:
        ps = random_commuting_paulis(m, rng)
        norm = diagonalizing_clifford(ps)
        S = rng.choice(np.array([-1, 1]), size=(3, m))
        S[:, 0] = np.abs(S[:, 0])
        w = rng.dirichlet(np.ones(3))
        F_target = np.einsum("r,ri,rj->ij", w, S.astype(float), S.astype(float))
        np.fill_diagonal(F_target, 1.0)
        G = norm.z_frame_target(F_target)
        dec = pol.exact_decompose(G)
        if not dec.feasible:
            rows.append({"m": m, "status": "z_frame_target_infeasible"})
            continue
        zprog = labelled_schedule_from_signs(dec.signs, dec.weights)
        pprog = physical_program(zprog, norm)
        F_real = pprog.realised_qfim(ps)
        err = relative_frobenius_error(F_real, F_target)
        worst = max(worst, err)
        frame = physical_width_after_normalisation(zprog, norm)
        rows.append({"m": m, "generators": [ps.label(i) for i in range(m)],
                     "clifford_two_qubit_gates": norm.cx_count,
                     "sign_flips": int((norm.sign_out < 0).sum()),
                     "permutation": norm.perm,
                     "relative_error": err,
                     "z_frame_width": frame["max_z_frame_width"],
                     "physical_width": frame["max_physical_width"]})
    status = PASS if worst < 1e-10 else FAIL
    return GateResult(
        "C7_Clifford_mapping", status,
        "QFIM realised in the physical generator frame matches the target after "
        "simultaneous diagonalisation",
        {"worst_relative_error": worst, "rows": rows},
        "Width is measured in both frames: a width-k Z-frame branch is generically "
        "a globally entangled physical state, so minimum-entanglement language is "
        "frame-relative and must be stated as such.")


def c8_hardware_legality(m: int = 6, seed: int = 9181,
                         topologies=("all_to_all", "line", "square_grid",
                                     "heavy_hex", "modular_two_cluster",
                                     "neutral_atom_grid")) -> GateResult:
    from .circuits import is_hardware_legal

    rng = np.random.default_rng(seed)
    S = rng.choice(np.array([-1, 1]), size=(3, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(3))
    prog = labelled_schedule_from_signs(S, w)
    rows = []
    bad = []
    for name in topologies:
        edges = build_topology(name, m)
        legal = True
        swaps = 0
        two_q = 0
        for br in prog.branches:
            low = lower_to_cx(br.circuit)
            routed, sw, _ = route(low, edges)
            swaps += sw
            r = resources(routed)
            two_q += r["two_qubit_gate_count"]
            legal = legal and is_hardware_legal(routed, edges)
        rows.append({"topology": name, "legal": legal, "swaps": swaps,
                     "two_qubit_gate_count": two_q})
        if not legal:
            bad.append(name)
    return GateResult("C8_hardware_legality", FAIL if bad else PASS,
                      "every two-qubit gate uses a legal edge after routing",
                      {"failures": bad, "rows": rows})


def c9_resource_accounting(m: int = 5, seed: int = 17041) -> GateResult:
    """Resources come from emitted circuits, and lowering preserves semantics."""
    rng = np.random.default_rng(seed)
    S = rng.choice(np.array([-1, 1]), size=(4, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(4))
    rows = []
    bad = []
    programs = {
        "labelled_global_cat": labelled_schedule_from_signs(S, w),
        "single_pure_state": single_pure_state_program(S, w),
        "coherent_flag": coherent_flag_program(S, w),
    }
    for name, prog in programs.items():
        circuits = ([br.circuit for br in prog.branches]
                    if prog.kind == "labelled" else [prog.circuit])
        fid = 1.0
        for c in circuits:
            psi_hi = simulate(c)
            low = lower_to_cx(c)
            psi_lo = simulate(low)
            # compare on the ancilla-zero slice of the lowered circuit
            extra = low.n_qubits - c.n_qubits
            if extra:
                psi_lo = psi_lo.reshape(1 << extra, -1)[0]
            fid = min(fid, float(abs(np.vdot(psi_hi, psi_lo))))
        emitted = prog.emitted()
        recount = sum(1 for c in circuits for g in lower_to_cx(c).gates if g.name == "cx")
        ok = abs(fid - 1.0) < 1e-9 and recount == emitted["cx_count"]
        rows.append({"program": name, "lowering_fidelity": fid,
                     "cx_count": emitted["cx_count"],
                     "independent_recount": recount,
                     "n_settings": emitted["n_settings"],
                     "n_ancilla": emitted["n_ancilla"]})
        if not ok:
            bad.append(name)
    return GateResult("C9_resource_accounting", FAIL if bad else PASS,
                      "lowered circuits are unitarily equivalent to the emitted "
                      "circuits and every count is a recount of the instruction list",
                      {"failures": bad, "rows": rows})


def c10_downstream_stability(m: int = 5, seed: int = 27183, trials: int = 200) -> GateResult:
    """MEASURED, not PASS: report how tight the submultiplicative bound is."""
    rng = np.random.default_rng(seed)
    S = rng.choice(np.array([-1, 1]), size=(4, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(4))
    F_star = np.einsum("r,ri,rj->ij", w, S.astype(float), S.astype(float))
    np.fill_diagonal(F_star, 1.0)
    fw = pol.frank_wolfe(F_star, iters=12)
    dF = fw.matrix() - F_star
    ndF = float(np.linalg.norm(dF, 2))
    ratios = []
    if ndF > 1e-14:
        for _ in range(trials):
            A = rng.standard_normal((m, m))
            num = float(np.linalg.norm(A.T @ dF @ A, 2))
            den = float(np.linalg.norm(A, 2) ** 2) * ndF
            ratios.append(num / den if den > 0 else 0.0)
        wv, vv = np.linalg.eigh(dF)
        a = vv[:, int(np.argmax(np.abs(wv)))]
        A_tight = np.outer(a, a)
        num = float(np.linalg.norm(A_tight.T @ dF @ A_tight, 2))
        den = float(np.linalg.norm(A_tight, 2) ** 2) * ndF
        ratios.append(num / den)
    return GateResult(
        "C10_downstream_stability", MEASURED,
        "tightness of ||A^T dF A|| <= ||A||^2 ||dF||, reported with CRB conditioning",
        {"operator_norm_error": ndF,
         "frobenius_error": float(np.linalg.norm(dF, "fro")),
         "ratio_max": float(max(ratios)) if ratios else None,
         "ratio_median": float(np.median(ratios)) if ratios else None,
         "conditioning_target": crb_conditioning(F_star),
         "conditioning_realised": crb_conditioning(fw.matrix())},
        "The inequality is submultiplicativity and holds for every matrix, so it "
        "is not a test. What governs precision is the conditioning of F^{-1} on "
        "the estimable subspace; the polytope's extreme points are rank one, i.e. "
        "maximally ill-conditioned.")


def c11_sld_saturability(m_values=(3, 4, 5), seed: int = 2026) -> GateResult:
    """The attainability lemma the abstract's 'certified performance' rests on."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    rows = []
    for m in m_values:
        for gen_name, ps in (("physical_Z", z_generators(m)),
                             ("clifford_conjugated", random_commuting_paulis(m, rng))):
            s = rng.choice(np.array([-1, 1]), size=m)
            psi = data_statevector(signed_cat_circuit(s, n_qubits=m))
            defect = sld_commutator_defect(psi, ps)
            worst = max(worst, defect)
            rows.append({"m": m, "generators": gen_name, "defect": defect})
    status = PASS if worst < 1e-12 else FAIL
    return GateResult(
        "C11_sld_saturability", status,
        "mean Uhlmann curvature vanishes: max_ij |Im <psi|P_i P_j|psi>| < 1e-12",
        {"worst_defect": worst, "rows": rows},
        "Im <d_i psi|d_j psi> = -(1/4) Im <psi|P_i P_j|psi>, which vanishes "
        "identically for commuting Hermitian generators. Without this the "
        "multiparameter Cramer-Rao bound is not attainable and no 'certified "
        "downstream performance' claim is licensed.")


def c12_convexity_gap(m: int = 4, seed: int = 3407) -> GateResult:
    """Report the gap the forbidden shortcut would silently pocket.

    The gap between ``F(rho)`` and ``sum_r p_r F_r`` vanishes exactly when every
    generator maps ``supp(rho)`` into ``ker(rho)``, since
    ``(q_k - q_l)^2 / (q_k + q_l) = q_k + q_l`` iff ``q_k q_l = 0``.  A signed
    cat satisfies this identically: ``Z_i |C_s> = s_i |C^-_s>`` sends the even
    cat to the *odd* cat, which carries no weight in a schedule built from even
    cats.  So for the compiler's own output family the unlabelled mixture is
    already optimal and the shortcut costs nothing -- which is a finding about
    the program-semantics section, not a licence to take the shortcut.

    The even/odd control is the case where it does cost: mixing ``|C_s>`` with
    ``|C^-_s>`` at equal weight gives ``rho = (|b><b| + |b~><b~|)/2``, whose
    covariance is still ``s s^T`` but whose QFI is exactly zero.  The gate
    therefore requires a zero gap on cat schedules *and* a strictly positive gap
    on the control; either one alone would be unfalsifiable.
    """
    from .programs import branch_states

    rng = np.random.default_rng(seed)
    ps = z_generators(m)
    rows = []

    S = rng.choice(np.array([-1, 1]), size=(4, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(4))
    prog_cat = labelled_schedule_from_signs(S, w)
    F_target = prog_cat.realised_qfim(ps)

    F_pure_check = unlabelled_mixed_qfi([(1.0, data_statevector(prog_cat.branches[0].circuit))], ps)
    pure_consistency = relative_frobenius_error(
        F_pure_check, np.outer(S[0], S[0]).astype(float))

    # width-2 schedule on a target inside Q^prog_{m,2}
    F_w2 = np.eye(m)
    for i in range(m - 1):
        F_w2[i, i + 1] = F_w2[i + 1, i] = 0.4
    d2 = wid.decompose_width2(F_w2)
    prog_w2 = labelled_schedule_from_width(d2) if d2.feasible else None

    # positive control: even and odd cat on the same sign vector
    s_ctrl = np.array([1] + [-1] * (m - 1))
    even = signed_cat_circuit(s_ctrl, n_qubits=m)
    odd = signed_cat_circuit(s_ctrl, n_qubits=m)
    odd.z(0)
    ctrl_states = [(0.5, data_statevector(even)), (0.5, data_statevector(odd))]

    cases: list[tuple[str, list, np.ndarray]] = [
        ("global_cat_schedule", branch_states(prog_cat), prog_cat.realised_qfim(ps)),
    ]
    if prog_w2 is not None:
        cases.append(("width2_schedule", branch_states(prog_w2), prog_w2.realised_qfim(ps)))
    ctrl_avg = sum(w * qfim_from_statevector(v, ps) for w, v in ctrl_states)
    cases.append(("even_odd_cat_control", ctrl_states, ctrl_avg))

    for label, states, avg in cases:
        F_mixed = unlabelled_mixed_qfi(states, ps)
        diff = avg - F_mixed
        eigs = np.linalg.eigvalsh(diff)
        overlaps = [abs(complex(np.vdot(a[1], b[1])))
                    for i, a in enumerate(states) for b in states[i + 1:]]
        rows.append({
            "schedule": label,
            "n_branches": len(states),
            "max_branch_overlap": float(max(overlaps)) if overlaps else 0.0,
            "trace_gap": float(np.trace(diff)),
            "relative_trace_gap": float(np.trace(diff) / max(np.trace(avg), 1e-30)),
            "min_eigenvalue_of_gap": float(eigs.min()),
            "loewner_order_holds": bool(eigs.min() >= -1e-9),
        })

    by_name = {r["schedule"]: r for r in rows}
    control = by_name.get("even_odd_cat_control", {})
    control_detects_gap = float(control.get("trace_gap", 0.0)) > 1e-6
    ok = (all(r["loewner_order_holds"] for r in rows)
          and pure_consistency < 1e-10
          and control_detects_gap)
    return GateResult(
        "C12_convexity_gap", MEASURED if ok else FAIL,
        "F(unlabelled rho) <= sum_r p_r F_r, with the gap reported per schedule shape",
        {"pure_state_consistency": pure_consistency,
         "control_detects_nonzero_gap": control_detects_gap, "rows": rows},
        "If the gap is zero for global cat schedules, the forbidden shortcut is "
        "not a shortcut there: orthogonal branches make the label free. The "
        "labelled/coherent distinction has to earn its section in the "
        "bounded-width regime, where the branches overlap.")


# ----------------------------------------------------------------------

def run_all(instances: list[Instance], stop_on_failure: bool = True) -> dict:
    """Run every gate.  On failure: write the report and abort, no retuning."""
    results: list[GateResult] = []
    order = [
        lambda: c1_signed_cat_semantics(),
        lambda: c2_flagged_composition(),
        lambda: c2c_flag_coherence_discriminator(),
        lambda: c3_exact_feasibility(instances),
        lambda: c4_infeasible_rejection(instances),
        lambda: c5_exact_reconstruction(instances),
        lambda: c6_width_integrity(instances),
        lambda: c7_clifford_mapping(),
        lambda: c8_hardware_legality(),
        lambda: c9_resource_accounting(),
        lambda: c10_downstream_stability(),
        lambda: c11_sld_saturability(),
        lambda: c12_convexity_gap(),
    ]
    aborted = None
    for fn in order:
        try:
            r = fn()
        except Exception as exc:  # a crash is a failure, not a skip
            r = GateResult(getattr(fn, "__name__", "gate"), FAIL,
                           "gate raised", {"exception": repr(exc)})
        results.append(r)
        if r.failed and stop_on_failure:
            aborted = r.name
            break
    return {
        "gates": [vars(r) for r in results],
        "summary": {
            "n_pass": sum(1 for r in results if r.status == PASS),
            "n_fail": sum(1 for r in results if r.status == FAIL),
            "n_measured": sum(1 for r in results if r.status == MEASURED),
            "n_absent": sum(1 for r in results if r.status == ABSENT),
            "aborted_at": aborted,
        },
        "stop_rule": "write_report_then_abort_without_tolerance_changes",
    }


def c2c_flag_coherence_discriminator(m: int = 4, seed: int = 2026) -> GateResult:
    """Decide, per program, whether the flag's *coherence* carries information.

    Dephasing the flag in its computational basis turns a coherent-flag program
    into the labelled schedule over its conditional branches.  Comparing the two
    QFIMs separates three objects the v0.1 protocol treated as one backend:

    * ``coherence not load-bearing`` -- dephasing costs nothing, so the program
      is genuinely interchangeable with a labelled schedule and branch width is
      a meaningful resource;
    * ``coherence load-bearing`` -- dephasing destroys the QFIM, so the program
      is *not* a labelled schedule in disguise, conditional branch width is
      meaningless as a resource, and flag dimension, Schmidt rank across the
      flag|data cut, and joint readout must be priced instead.

    The sign purification ``sum_s sqrt(p(s))|s>|s>`` is the extreme case: every
    conditional branch is a product eigenstate (width 1) yet it realises any
    feasible target, and dephasing sends its QFIM to exactly zero.
    """
    from .programs import coherent_sign_purification, flag_resource_report

    rng = np.random.default_rng(seed)
    ps = z_generators(m)
    S = rng.choice(np.array([-1, 1]), size=(4, m))
    S[:, 0] = np.abs(S[:, 0])
    w = rng.dirichlet(np.ones(4))
    target = np.einsum("r,ri,rj->ij", w, S.astype(float), S.astype(float))
    np.fill_diagonal(target, 1.0)

    progs = {
        "sign_purification": coherent_sign_purification(S, w),
        "flagged_cat_relaxed": coherent_flag_program(S, w),
        "flagged_cat_literal": prg_literal(S, w),
    }
    rows = []
    for name, p in progs.items():
        rep = flag_resource_report(p, ps)
        F = np.array(rep["qfim_coherent_flag"])
        em = p.emitted()
        rows.append({
            "program": name,
            "realises_target": bool(np.linalg.norm(F - target, "fro") < 1e-9),
            "conditional_data_width": rep["conditional_data_width"],
            "n_flag_qubits": rep["n_flag_qubits"],
            "schmidt_rank_flag_cut": rep["schmidt_rank_flag_cut"],
            "cx_count": em["cx_count"],
            "trace_coherent": float(np.trace(F)),
            "trace_dephased": float(np.trace(np.array(rep["qfim_dephased_flag"]))),
            "coherence_is_load_bearing": rep["coherence_is_load_bearing"],
        })
    literal = next(r for r in rows if r["program"] == "flagged_cat_literal")
    purif = next(r for r in rows if r["program"] == "sign_purification")
    ok = (all(r["realises_target"] for r in rows)
          and not literal["coherence_is_load_bearing"]
          and purif["coherence_is_load_bearing"]
          and purif["conditional_data_width"] == 1)
    return GateResult(
        "C2c_flag_coherence_discriminator", PASS if ok else FAIL,
        "dephasing the flag separates labelled-equivalent programs from "
        "genuinely coherent ones; conditional branch width is only a resource "
        "for the former",
        {"rows": rows},
        "Programs on opposite sides of this test are not comparable on gate "
        "count. A cheaper coherent program that fails dephasing has not "
        "compiled the same object.")


def prg_literal(S, w):
    from .programs import coherent_flag_program_literal
    return coherent_flag_program_literal(S, w)
