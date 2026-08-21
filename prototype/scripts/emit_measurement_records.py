"""Machine-generate the canonical measurement and noise records.

Every field is computed here; nothing is transcribed by hand.  Schema follows the
v0.3 handoff's canonical-record requirements plus the readout-artifact fields
required by the Section 11 review.  Run from ``prototype/``:

    PYTHONPATH=src python3 scripts/emit_measurement_records.py
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from covq.instances import matching_target
from covq.measurement import (cfi_of_readout, cfi_of_readout_mixed, compile_branch_readout,
                              compile_schedule_readout, covariance_surrogate, dephase,
                              depolarize, mixed_state_qfim, readout_contract, _rotate)
from covq.paulis import z_generators
from covq.programs import labelled_schedule_from_width, signed_cat_circuit
from covq.qfim import qfim_from_statevector
from covq.sim import data_statevector
from covq.width import decompose_width2

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "covq.measurement.record/0.1"
TOL = {"attainment_abs": 1e-9, "closed_form_abs": 1e-9, "pinching_abs": 1e-10}


def _hash(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=_json).encode()).hexdigest()


def _json(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    raise TypeError(type(o))


def _git(*args, default="unknown"):
    try:
        return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return default


SOURCE_PATHS = ("prototype/src", "prototype/tests", "prototype/scripts",
                "STATUS.md", "REGISTRATION.md", "docs", "charter")


def _source_dirty() -> bool:
    """Uncommitted modifications to the **source** tree.

    Deliberately excludes ``results/``.  A record's provenance answers "which
    source produced this?", and every emitter necessarily modifies ``results/``
    while running -- so a whole-tree check makes the first emitter of a run
    report clean and every later one report dirty, which is an artefact of
    ordering rather than a fact about the code.  Scoping the check to source
    paths makes the flag mean what it is read as meaning.
    """
    return bool(_git("status", "--porcelain", "--", *SOURCE_PATHS))


def _provenance() -> dict:
    return {
        "schema_version": SCHEMA,
        "source_commit": _git("rev-parse", "HEAD"),
        "dirty": _source_dirty(),
        "dirty_scope": list(SOURCE_PATHS),
        "tolerances": TOL,
        "qfim_evaluator": "covq.qfim.qfim_from_statevector / covq.measurement.mixed_state_qfim",
        "cfi_evaluator": "covq.measurement.cfi_of_readout(_mixed)",
    }


def _write(rel: str, payload: dict) -> Path:
    path = ROOT / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {**_provenance(), **payload}
    payload["record_hash"] = _hash({k: v for k, v in payload.items() if k != "record_hash"})
    path.write_text(json.dumps(payload, indent=2, default=_json) + "\n")
    return path


# ----------------------------------------------------------------------

def ideal_block_readout() -> dict:
    """C10a: block-local attainability, per block width."""
    cases = []
    for k in (2, 3, 4):
        ps = z_generators(k)
        for phi in (0.0, 0.37, 1.3):
            psi = _rotate(data_statevector(signed_cat_circuit(np.ones(k, int), n_qubits=k)),
                          ps, np.array([phi] + [0.0] * (k - 1)))
            alphas = compile_branch_readout(psi, ps, [tuple(range(k))])
            fq = qfim_from_statevector(psi, ps)
            cfi = cfi_of_readout(psi, ps, alphas)
            (rec,) = readout_contract(psi, [tuple(range(k))])
            cases.append({
                "block_width": k, "operating_point": [phi] + [0.0] * (k - 1),
                "block_sign_vector": [1] * k, "pivot": rec["pivot"],
                "local_analyzer_angles": alphas, "aggregate_analyzer_phase": rec["alpha"],
                "readout_regularity_margin": rec["regularity_margin"],
                "fringe_visibility": rec["visibility"],
                "parity_postprocessing_map": "product of outcome signs over the block",
                "readout_2q_gate_count": 0, "readout_ancilla_count": 0,
                "qfim": fq, "declared_readout_cfi": cfi,
                "attainment_gap": float(np.linalg.norm(fq - cfi, 2)),
                "readout_hash": _hash(list(np.round(alphas, 12))),
            })
    worst = max(c["attainment_gap"] for c in cases)
    return {"gate": "C10a", "alias": "M1", "cases": cases, "worst_attainment_gap": worst,
            "status": "DONE" if worst < TOL["attainment_abs"] else "FAILS"}


def schedule_attainability() -> dict:
    """C10b: retained-label schedule CFI additivity."""
    rng = np.random.default_rng(17)
    cases = []
    for m in (3, 4, 5, 6):
        inst = matching_target(m, np.random.default_rng(m * 13), load=0.8)
        dec = decompose_width2(inst.F)
        if not dec.feasible:
            continue
        prog = labelled_schedule_from_width(dec)
        ps = z_generators(m)
        for name, theta in (("origin", np.zeros(m)), ("generic", rng.uniform(-2.0, 2.0, m))):
            rep = compile_schedule_readout(prog, ps, theta)
            cases.append({
                "m": m, "operating_point_name": name, "operating_point": theta,
                "target_hash": _hash(np.round(inst.F, 12)), "n_branches": len(rep["branches"]),
                "qfim_schedule": rep["F_schedule"], "declared_readout_cfi": rep["CFI_schedule"],
                "schedule_gap": rep["schedule_gap"],
                "all_branches_loewner_dominated": rep["all_branches_dominated"],
            })
    worst = max(c["schedule_gap"] for c in cases)
    return {"gate": "C10b", "alias": "M1_schedule", "cases": cases,
            "worst_schedule_gap": worst,
            "status": "DONE" if worst < TOL["attainment_abs"] else "FAILS"}


def readout_singularities() -> dict:
    """C10c: singular hyperplanes of an unmatched readout, and their repair."""
    dec = decompose_width2(matching_target(4, np.random.default_rng(11), load=0.8).F)
    prog = labelled_schedule_from_width(dec)
    ps = z_generators(4)
    signs = [[int(v) for v in s] for _, _, s in dec.branches]
    cases = []
    for name, theta in (("origin", np.zeros(4)), ("uniform_ray", np.full(4, 0.3)),
                        ("generic", np.array([0.7, -0.35, 1.1, 0.2]))):
        rep = compile_schedule_readout(prog, ps, theta)
        cases.append({
            "operating_point_name": name, "operating_point": theta,
            "matched_readout_gap": rep["schedule_gap"],
            "matched_readout_trace": rep["compiled_total_information"],
            "fixed_all_x_trace": rep["fixed_x_total_information"],
            "fixed_all_x_gap": float(np.linalg.norm(rep["F_schedule"]
                                                    - rep["CFI_schedule_fixed_x"], 2)),
            "qfim_trace": float(np.trace(rep["F_schedule"])),
        })
    return {"gate": "C10c", "alias": "M2", "branch_sign_vectors": signs,
            "singular_set": "phi_B(theta) - A_B in pi Z",
            "regularity_claim": "zero first-order score; nonregular, not distinguishability-free",
            "cases": cases,
            "status": "DONE" if all(c["matched_readout_gap"] < TOL["attainment_abs"]
                                    for c in cases) else "FAILS"}


def qfi_covariance_separation() -> dict:
    """C15a/b/c: mixed QFIM, the three-matrix separation, and pinching."""
    rows = []
    ps2 = z_generators(2)
    psi = _rotate(data_statevector(signed_cat_circuit(np.array([1, 1]), n_qubits=2)),
                  ps2, np.array([0.4, 0.0]))
    alphas = compile_branch_readout(psi, ps2, [(0, 1)])
    rho0 = np.outer(psi, psi.conj())
    ideal = qfim_from_statevector(psi, ps2)
    for channel, fn in (("dephasing", dephase), ("depolarizing", depolarize)):
        for p in (0.0, 0.1, 0.2, 0.35, 0.5):
            rho = fn(rho0, p)
            rows.append({
                "noise_channel": channel, "noise_parameters": {"p": p},
                "state_hash": _hash(np.round(rho.view(float), 10)),
                "generator_hash": _hash([ps2.label(i) for i in range(ps2.m)]),
                "readout_hash": _hash(list(np.round(alphas, 12))),
                "operating_point": [0.4, 0.0], "analyzer_angles": alphas,
                "regularity_margin": readout_contract(psi, [(0, 1)])[0]["regularity_margin"],
                "ideal_qfim": ideal,
                "noisy_generator_covariance": covariance_surrogate(rho, ps2),
                "mixed_state_qfim": mixed_state_qfim(rho, ps2),
                "declared_readout_cfi": cfi_of_readout_mixed(rho, ps2, alphas),
            })
    cov_flat = [r for r in rows if r["noise_channel"] == "dephasing"]
    cov_invariant = max(float(np.abs(r["noisy_generator_covariance"]
                                     - cov_flat[0]["noisy_generator_covariance"]).max())
                        for r in cov_flat)
    gap = max(abs(float(np.trace(r["declared_readout_cfi"]) - np.trace(r["mixed_state_qfim"])))
              for r in rows if r["noise_channel"] == "dephasing")
    return {"gate": ["C15a", "C15b", "C15c"], "alias": ["M3_qfi", "M3_sep", "M3_dephase"],
            "proposition": "Cov_{D(rho)}(P) = Cov_rho(P) and F_Q(D(rho)) = 0",
            "dephasing_covariance_max_drift": cov_invariant,
            "dephasing_cfi_minus_qfim_max": gap,
            "note": "the CFI column is the declared matched-quadrature readout, which "
                    "attains the mixed-state QFIM at every visibility",
            "rows": rows,
            "status": "DONE" if cov_invariant < 1e-12 and gap < TOL["attainment_abs"] else "FAILS"}


def estimator_efficiency() -> dict:
    """Eq (112): does the branch-conditioned MLE saturate the readout's CRB?"""
    from covq.estimator import efficiency_report, model_from_decomposition

    rng = np.random.default_rng(9)
    dec = decompose_width2(matching_target(3, np.random.default_rng(39), load=0.8).F)
    model = model_from_decomposition(dec)
    cases = []
    for theta_name, theta in (("origin", np.zeros(3)),
                              ("generic", np.array([0.21, -0.14, 0.33]))):
        for n_shots in (500, 5_000, 50_000):
            rep = efficiency_report(model, theta, n_shots, 4_000, rng)
            cases.append({
                "operating_point_name": theta_name, "operating_point": theta,
                "n_shots": rep["n_shots"], "n_reps": rep["n_reps"], "rank": rep["rank"],
                "max_abs_bias": rep["max_abs_bias"],
                "bias_standard_error": rep["bias_standard_error"],
                "bias_within_4_stderr": rep["bias_within_4_stderr"],
                "efficiency_eigenvalues": rep["efficiency_eigenvalues"],
                "marchenko_pastur_band": rep["mp_band"],
                "efficiency_within_mp_band": rep["efficiency_within_mp_band"],
            })
    ok = all(c["efficiency_within_mp_band"] and c["bias_within_4_stderr"] for c in cases)
    return {"gate": "C10d", "alias": "M4",
            "claim": "branch-conditioned MLE over block-parity counts is consistent and "
                     "attains F_Pi^+ / N on the identifiable quotient",
            "estimator": "maximum likelihood, BFGS, sufficient statistic = one parity "
                         "count per (branch, block)",
            "cases": cases, "status": "DONE" if ok else "FAILS"}


def adaptive_recentering_record() -> dict:
    """Where the operating point comes from, without an oracle."""
    from covq.estimator import adaptive_efficiency_report, model_from_decomposition

    rng = np.random.default_rng(3)
    dec = decompose_width2(matching_target(3, np.random.default_rng(39), load=0.8).F)
    model = model_from_decomposition(dec)
    theta = np.array([0.6, -0.4, 0.9])
    cases = []
    for fraction in (0.05, 0.1, 0.2, 0.4):
        for n_shots in (2_000, 20_000):
            rep = adaptive_efficiency_report(model, theta, n_shots, 800, rng,
                                             pilot_fraction=fraction)
            cases.append({
                "pilot_fraction": fraction, "n_shots": rep["n_shots"],
                "n_reps": rep["n_reps"], "operating_point": theta,
                "max_abs_bias": rep["max_abs_bias"],
                "bias_standard_error": rep["bias_standard_error"],
                "worst_efficiency_ratio": rep["worst_efficiency_ratio"],
                "naive_discard_pilot_penalty": 1.0 / (1.0 - fraction),
                "median_regularity_margin": rep["median_regularity_margin"],
                "marchenko_pastur_band": rep["mp_band"],
            })
    ok = all(c["worst_efficiency_ratio"] < c["naive_discard_pilot_penalty"] + 0.35
             for c in cases)
    return {"gate": "C10e", "alias": "M5",
            "protocol": "two-stage: pilot split between A=0 and A=pi/2 fixes each block "
                        "phase including sign; main stage at matched quadrature; single "
                        "MLE over all counts, seeded by least squares on the design matrix",
            "reference_bound": "oracle F^+/N",
            "cases": cases, "status": "DONE" if ok else "FAILS"}


if __name__ == "__main__":
    for rel, builder in (
        ("results/measurements/ideal_block_readout.json", ideal_block_readout),
        ("results/measurements/schedule_attainability.json", schedule_attainability),
        ("results/measurements/readout_singularities.json", readout_singularities),
        ("results/noise/qfi_covariance_separation.json", qfi_covariance_separation),
        ("results/measurements/estimator_efficiency.json", estimator_efficiency),
        ("results/measurements/adaptive_recentering.json", adaptive_recentering_record),
    ):
        payload = builder()
        path = _write(rel, payload)
        print(f"{payload['status']:<6} {rel}")
