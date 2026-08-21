"""Eq (111) campaign records.  New namespace; the six earlier records are immutable."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from emit_measurement_records import _hash, _write  # noqa: E402

from covq.measurement import covariance_surrogate, mixed_state_qfim
from covq.noise import (BlockLocalNoise, branch_template, delta_edge,
                        noise_aware_floor_compile, pair_channel, pair_template,
                        price_branch_bruteforce, price_branch_noisy, product_template,
                        setting_cost_crossover)
from covq.paulis import z_generators

ROOT = Path(__file__).resolve().parents[2]


def _edges(m):
    return [(i, j) for i in range(m) for j in range(i + 1, m)]


def _noise(m, edges, deph=0.02, depol=0.0, qe=0.05, idle=0.01, conf=None):
    return BlockLocalNoise(dephasing={q: deph for q in range(m)},
                           depolarizing={q: depol for q in range(m)},
                           edge_depolarizing={e: qe for e in edges},
                           idle_dephasing=idle,
                           readout_confusion=conf or {})


def _noise_fields(m, edges, nz, theta):
    return {
        "noise_channel_and_location": {
            "local_dephasing": "per data qubit, before readout",
            "local_depolarizing": "per data qubit",
            "edge_depolarizing": "two-qubit, after each pair preparation",
            "idle_dephasing": "unmatched qubits, during the single pair-prep layer",
            "readout_confusion": "per qubit, asymmetric (e0, e1)",
        },
        "noise_parameters": {"dephasing": nz.dephasing, "depolarizing": nz.depolarizing,
                             "idle_dephasing": nz.idle_dephasing,
                             "readout_confusion": {str(k): v
                                                   for k, v in nz.readout_confusion.items()}},
        "edge_error_map": {str(e): nz.e_depol(e) for e in edges},
        "operating_point": theta,
        "symmetric_readout": nz.is_symmetric_readout,
    }


def noisy_branch_templates() -> dict:
    """N1-N4: noiseless limit, analytic quadrature value, pinching, edge additivity."""
    rng = np.random.default_rng(31)
    cases = []
    for m in (3, 4, 5):
        edges = _edges(m)
        nz = BlockLocalNoise(
            dephasing={q: float(rng.uniform(0, 0.1)) for q in range(m)},
            depolarizing={q: float(rng.uniform(0, 0.05)) for q in range(m)},
            edge_depolarizing={e: float(rng.uniform(0, 0.12)) for e in edges},
            idle_dephasing=float(rng.uniform(0, 0.05)))
        theta = rng.uniform(-1, 1, m)
        matching = [(0, 1)] if m == 3 else [(0, 1), (2, 3)]
        worst = 0.0
        for signs in ([1] * len(matching), [-1] * len(matching)):
            got = branch_template(matching, signs, m, theta, nz)
            want = product_template(m, theta, nz, idle=True) + sum(
                delta_edge(e, s, m, theta, nz) for e, s in zip(matching, signs))
            worst = max(worst, float(np.abs(got - want).max()))
        cases.append({"m": m, **_noise_fields(m, edges, nz, theta),
                      "matching": matching, "edge_additivity_residual": worst})
    # N1 and N2 analytic anchors.
    ideal = branch_template([(0, 1), (2, 3)], [1, -1], 4,
                            np.array([0.3, -0.2, 0.5, 0.1]), BlockLocalNoise())
    n2 = []
    for p in (0.0, 0.05, 0.15, 0.3):
        blk = pair_template((0.4, 0.0), (0, 1), 1, BlockLocalNoise(dephasing={0: p, 1: p}))
        v = (1.0 - 2.0 * p) ** 2
        n2.append({"p": p, "visibility": v, "declared_readout_cfi": blk,
                   "analytic_v_squared": v * v,
                   "residual": float(np.abs(blk - v * v * np.ones((2, 2))).max())})
    dead = pair_channel((0.4, 0.0), (0, 1), 1, BlockLocalNoise(dephasing={0: 0.5, 1: 0.5}))
    clean = pair_channel((0.4, 0.0), (0, 1), 1, BlockLocalNoise())
    ps = z_generators(2)
    ok = (max(c["edge_additivity_residual"] for c in cases) < 1e-12
          and max(r["residual"] for r in n2) < 1e-9)
    return {"gates": ["N1", "N2", "N3", "N4"], "cases": cases,
            "n1_noiseless_template_residual": float(np.abs(
                ideal - np.array([[1, 1, 0, 0], [1, 1, 0, 0],
                                  [0, 0, 1, -1], [0, 0, -1, 1]], float)).max()),
            "n2_dephased_quadrature": n2,
            "n3_pinching": {
                "covariance_drift": float(np.abs(covariance_surrogate(dead, ps)
                                                 - covariance_surrogate(clean, ps)).max()),
                "mixed_state_qfim_max": float(np.abs(mixed_state_qfim(dead, ps)).max()),
                "declared_readout_cfi_max": float(np.abs(pair_template(
                    (0.4, 0.0), (0, 1), 1, BlockLocalNoise(dephasing={0: 0.5, 1: 0.5}))).max()),
            },
            "status": "DONE" if ok else "FAILS"}


def noisy_pricing_oracle() -> dict:
    """N5: does block-local noise preserve matching-based pricing?"""
    rng = np.random.default_rng(7)
    cases, worst = [], 0.0
    for m in (3, 4, 5):
        edges = _edges(m)
        for _ in range(4):
            theta = rng.uniform(-1, 1, m)
            nz = BlockLocalNoise(
                dephasing={q: float(rng.uniform(0, 0.12)) for q in range(m)},
                depolarizing={q: float(rng.uniform(0, 0.05)) for q in range(m)},
                edge_depolarizing={e: float(rng.uniform(0, 0.15)) for e in edges},
                idle_dephasing=float(rng.uniform(0, 0.06)))
            a = rng.standard_normal((m, m))
            Q = a @ a.T / m
            costs = {e: float(rng.uniform(0.0, 0.35)) for e in edges}
            oracle = price_branch_noisy(Q, m, edges, theta, nz, costs)
            brute = price_branch_bruteforce(Q, m, edges, theta, nz, costs)
            gap = abs(oracle["value"] - brute["value"])
            worst = max(worst, gap)
            cases.append({"m": m, **_noise_fields(m, edges, nz, theta),
                          "oracle_value": oracle["value"], "brute_value": brute["value"],
                          "gap": gap, "winning_branch_kind": oracle["branch_kind"],
                          "oracle_matching": oracle["matching"],
                          "brute_matching": brute["matching"]})
    return {"gate": "N5",
            "theorem": "block-local noise reweights edges but preserves the maximum-weight "
                       "matching oracle; the pairless branch is a separate column because "
                       "it activates no pair and so carries no idle dephasing",
            "scope_excluded": ["crosstalk", "correlated branch noise", "route collisions",
                               "coherent errors coupling distinct blocks"],
            "cases": cases, "n_cases": len(cases), "worst_gap": worst,
            "status": "DONE" if worst < 1e-9 else "FAILS"}


def noise_aware_floor_compiler() -> dict:
    """N6: convex exposure problem, primal/dual, and the entanglement threshold."""
    m = 4
    edges = _edges(m)
    G = np.full((m, m), 0.6)
    np.fill_diagonal(G, 1.2)
    theta = np.zeros(m)
    sweep = []
    for qe in (0.0, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30):
        nz = _noise(m, edges, qe=qe)
        r = noise_aware_floor_compile(G, m, edges, theta, nz,
                                      costs={e: 0.1 for e in edges}, c0=1.0)
        sweep.append({"edge_depolarizing": qe, **_noise_fields(m, edges, nz, theta),
                      "primal_cost": r["cost"], "dual_bound": r["dual_bound"],
                      "optimality_gap": r["relative_gap"],
                      "downstream_floor_margin": r["floor_slack_min_eig"],
                      "settings": r["n_settings_used"],
                      "entangled_settings": r["n_entangled_settings_used"],
                      "status": r["status"]})
    on = [s["edge_depolarizing"] for s in sweep if s["entangled_settings"] > 0]
    off = [s["edge_depolarizing"] for s in sweep if s["entangled_settings"] == 0]
    ok = all(s["optimality_gap"] < 1e-9 and s["status"] == "solved" for s in sweep)
    return {"gate": "N6", "target": G, "sweep": sweep,
            "entanglement_threshold_between": [max(on) if on else None,
                                               min(off) if off else None],
            "worst_optimality_gap": max(s["optimality_gap"] for s in sweep),
            "status": "DONE" if ok else "FAILS"}


def setting_cost_phase_diagram() -> dict:
    """Turn the unpriced-settings caveat into an amortization boundary."""
    src = json.loads((ROOT / "results/baselines_k3_quest.json").read_text())
    rows = []
    for c in src["cases"]:
        if "skipped" in c or c["quest"]["stop_reason"] != "converged":
            continue
        entry = {"instance": c["instance"],
                 "q_covq": c["covq"]["settings"],
                 "expected_cx_covq": c["covq"]["expected_cx_per_shot"],
                 "q_quest": c["quest"]["settings"],
                 "cx_quest": c["quest"]["cx_per_shot"], "crossover": {}}
        for c_setup in (10, 100, 1000, 10000):
            r = setting_cost_crossover(entry["expected_cx_covq"], entry["q_covq"],
                                       entry["cx_quest"], entry["q_quest"], c_setup)
            entry["crossover"][str(c_setup)] = r["crossover_shots"]
        rows.append(entry)
    worst = max(r["crossover"]["10000"] for r in rows)
    return {"gate": "K3-amortization",
            "cost_unit": "one two-qubit gate; c_setup quoted in the same unit",
            "switching_policy": "batched execution, N_switch = q - 1",
            "formula": "N* = (q_C - q_Q) c_setup / (c_Q - cbar_C)",
            "caveat": "a fixed charge per distinct setting is a cardinality cost; pricing "
                      "it turns the compiler into a mixed-integer conic program and voids "
                      "the oracle-polynomial claim, so it is reported here as a boundary "
                      "rather than folded into the convex objective",
            "rows": rows, "worst_crossover_at_c_setup_1e4": worst,
            "status": "DONE"}


if __name__ == "__main__":
    for rel, builder in (
        ("results/noise/noisy_branch_templates.json", noisy_branch_templates),
        ("results/noise/noisy_pricing_oracle.json", noisy_pricing_oracle),
        ("results/noise/noise_aware_floor_compiler.json", noise_aware_floor_compiler),
        ("results/noise/setting_cost_phase_diagram.json", setting_cost_phase_diagram),
    ):
        payload = builder()
        _write(rel, payload)
        print(f"{payload['status']:<6} {rel}")
