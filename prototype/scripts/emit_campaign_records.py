"""N7, N9, N10: deployable readout, fixed setting cost, and the equal-accounting campaign."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from emit_measurement_records import _hash, _write  # noqa: E402

from covq.measurement import mixed_state_qfim
from covq.circuits import lower_to_cx, resources
from covq.noise import (FROZEN_PILOT_POLICY, BlockLocalNoise, branch_template,
                        deployable_exposure, fixed_setting_cost_compile, floor_margin,
                        noise_aware_floor_compile, pair_channel, pilot_recentered_block,
                        product_template, simulate_circuit_noisy, simulate_noisy_state)
from covq.paulis import z_generators
from covq.polytope import exact_decompose
from covq.programs import single_pure_state_program
from covq.quest import canonicalise_rotations, moment_constraints, quest_circuit, quest_published
from covq.sim import simulate
from covq.width import decompose_width2


def _nz(m, edges, qe):
    return BlockLocalNoise(dephasing={q: 0.02 for q in range(m)},
                           edge_depolarizing={e: qe for e in edges},
                           idle_dephasing=0.01)


def adaptive_readout_performance() -> dict:
    """N7: the deployable pilot arm against the oracle-angle upper bound."""
    rng = np.random.default_rng(5)
    total = 20_000
    arms = []
    for tag, noise in (("clean", BlockLocalNoise()),
                       ("light", BlockLocalNoise(dephasing={0: 0.03, 1: 0.03},
                                                 edge_depolarizing={(0, 1): 0.04})),
                       ("heavy", BlockLocalNoise(dephasing={0: 0.12, 1: 0.12},
                                                 edge_depolarizing={(0, 1): 0.15}))):
        rho = pair_channel((0.4, 0.0), (0, 1), 1, noise)
        rows = []
        for fraction in (0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40):
            n_pilot = max(2, int(fraction * total))
            rep = pilot_recentered_block(rho, (0, 1), noise, n_pilot, total - n_pilot,
                                         rng, n_reps=300)
            rows.append({"pilot_fraction": fraction, "total_shots": total,
                         "deployable_over_oracle": rep["ratio"]})
        best = max(rows, key=lambda r: r["deployable_over_oracle"])
        arms.append({"noise_level": tag, "sweep": rows,
                     "optimal_pilot_fraction": best["pilot_fraction"],
                     "best_ratio": best["deployable_over_oracle"]})
    ok = all(r["deployable_over_oracle"] <= 1.0 + 1e-9 for a in arms for r in a["sweep"])
    return {"gate": "N7",
            "claim": "oracle analyzer angles are an upper bound; the deployable arm pays "
                     "pilot-estimation error (vanishing in pilot size) plus off-quadrature "
                     "pilot shots (proportional to pilot fraction)",
            "arms": arms, "status": "DONE" if ok else "FAILS"}


def fixed_setting_cost_solver() -> dict:
    """N9: greedy support selection against exhaustive enumeration.  It fails."""
    m = 3
    edges = [(0, 1), (0, 2), (1, 2)]
    G = np.full((m, m), 0.5)
    np.fill_diagonal(G, 1.1)
    noise = _nz(m, edges, 0.04)
    rows = []
    for lam in (0.0, 0.1, 0.2, 0.3, 1.0):
        g = fixed_setting_cost_compile(G, m, edges, np.zeros(m), noise, lam,
                                       costs={e: 0.1 for e in edges}, method="greedy")
        e = fixed_setting_cost_compile(G, m, edges, np.zeros(m), noise, lam,
                                       costs={e: 0.1 for e in edges},
                                       method="exhaustive", max_support=3)
        rows.append({"lambda_q": lam, "greedy_cost": g["total_cost"],
                     "greedy_settings": g["n_settings"],
                     "exhaustive_cost": e["total_cost"],
                     "exhaustive_settings": e["n_settings"],
                     "greedy_excess_percent": 100 * (g["total_cost"] - e["total_cost"])
                     / e["total_cost"]})
    worst = max(r["greedy_excess_percent"] for r in rows)
    return {"gate": "N9",
            "verdict": "FAILS as an exactness claim, by design of the check: forward "
                       "selection commits to the best single column, which need not "
                       "belong to the best pair.  It is exact only once the cardinality "
                       "penalty makes one setting genuinely optimal -- outside the "
                       "multi-setting regime the schedule exists to exploit.",
            "note": "cardinality cost makes this mixed-integer conic; no oracle-polynomial "
                    "claim is made for it, and it is kept out of the convex compiler",
            "rows": rows, "worst_greedy_excess_percent": worst,
            "status": "FAILS"}


def quest_operational_comparison() -> dict:
    """N10: same target, same noise map, same accounting, four arms.

    QUEST is given its **noisy QFI** as an optimistic upper bound -- no readout
    is compiled for it -- while CovQ reports the CFI of its actual emitted
    readout.  That is deliberately unkind to CovQ, and it is the only honest
    form: handing QUEST an unpriced, unspecified optimal POVM and calling the
    result operational is exactly the asymmetry the handoff warns against.
    Two-qubit depolarization is charged per emitted two-qubit gate on both sides.
    """
    m = 4
    edges = [(i, j) for i in range(m) for j in range(i + 1, m)]
    theta = np.zeros(m)
    G = np.full((m, m), 0.6)
    np.fill_diagonal(G, 1.2)
    ps = z_generators(m)

    F_target = np.full((m, m), 0.4)
    np.fill_diagonal(F_target, 1.0)
    qr = quest_published(moment_constraints(ps, F_target), m, variant='tE',
                         max_depth=80, tol=1e-13)
    _qres = resources(lower_to_cx(quest_circuit(qr, m)))
    quest_cx = _qres["cx_count"]
    quest_depth = _qres["two_qubit_depth"]

    dec = exact_decompose(F_target)
    sparse_prog = single_pure_state_program(dec.signs, dec.weights)
    sparse_cx = resources(lower_to_cx(sparse_prog.circuit))["cx_count"]

    rows = []
    for qe in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20):
        noise = _nz(m, edges, qe)
        covq = noise_aware_floor_compile(G, m, edges, theta, noise,
                                         costs={e: 0.0 for e in edges}, c0=1.0)
        prod = 1.0 / floor_margin(product_template(m, theta, noise, idle=False), G)

        def shots_from(rho):
            g = floor_margin(mixed_state_qfim(rho / np.trace(rho).real, ps), G)
            return (1.0 / g) if g > 1e-15 else None

        # Route QUEST through the *same* lowering and per-emitted-gate noise
        # path as every other circuit arm.  Charging noise per two-qubit
        # rotation instead of per emitted CX under-billed QUEST by exactly 2x,
        # since each two-qubit Pauli rotation lowers to two CX.
        q_shots = shots_from(simulate_circuit_noisy(quest_circuit(qr, m), noise))
        s_shots = shots_from(simulate_circuit_noisy(sparse_prog.circuit, noise))
        # Two independent channels for the deployable template.  "exact"
        # enumerates the binomial pilot support and does no sampling at all;
        # "mc" draws pilot counts.  Their agreement bounds both the sampler's
        # noise and any bias in how the expectation was formed -- pilot-seed
        # spread alone would only have measured the former.  The deterministic
        # channel is the reported value.
        deploy = deployable_exposure(G, m, edges, theta, noise, covq["branches"],
                                     costs={e: 0.0 for e in edges}, c0=1.0,
                                     method="exact")
        deploy_mc = deployable_exposure(G, m, edges, theta, noise, covq["branches"],
                                        costs={e: 0.0 for e in edges}, c0=1.0,
                                        method="mc")
        rows.append({
            "edge_depolarizing": qe,
            "covq_deployable_fixed_pilot_exposure": deploy.get("cost"),
            "covq_oracle_angle_exposure": covq["cost"],
            "deployable_over_oracle": (deploy["cost"] / covq["cost"]
                                       if deploy.get("cost") else None),
            "deployable_floor_slack_min_eig": deploy.get("floor_slack_min_eig"),
            "deployable_status": deploy["status"],
            "covq_deployable_channel": "exact (binomial enumeration, no sampling)",
            "covq_deployable_monte_carlo_crosscheck": deploy_mc.get("cost"),
            "covq_deployable_channel_disagreement_relative": (
                abs(deploy["cost"] - deploy_mc["cost"]) / deploy["cost"]
                if deploy.get("cost") else None),
            "covq_shots_emitted_readout_cfi": covq["cost"],
            "covq_entangled_settings": covq["n_entangled_settings_used"],
            "covq_expected_two_qubit_gates": float(sum(
                n * len(M) for n, (M, _) in zip(covq["exposures"], covq["branches"]))),
            "product_only_shots": prod,
            "quest_shots_noisy_qfi_upper_bound": q_shots,
            "quest_two_qubit_gates_total": (qr.two_qubit_rotations * q_shots
                                            if q_shots else None),
            "sparse_caratheodory_shots_noisy_qfi_upper_bound": s_shots,
        })
    disagree = [r["covq_deployable_channel_disagreement_relative"] for r in rows
                if r.get("covq_deployable_channel_disagreement_relative") is not None]
    margins = [r["quest_shots_noisy_qfi_upper_bound"]
               / r["covq_deployable_fixed_pilot_exposure"] for r in rows
               if r.get("quest_shots_noisy_qfi_upper_bound")]
    crossing = [abs(mm - 1.0) for mm in margins if mm > 1.0]
    return {"gate": "N10",
            "deployable_uncertainty": {
                "channels": ["exact binomial enumeration", "Monte-Carlo pilot draws"],
                "max_relative_disagreement": max(disagree) if disagree else None,
                "exact_channel_grid_discretisation": 1.5e-6,
                "narrowest_crossing_margin": min(crossing) if crossing else None,
                "margin_over_uncertainty": (
                    (min(crossing) / max(disagree)) if crossing and disagree else None),
                "note": "pilot-seed spread alone measures only sampler noise; the "
                        "deterministic channel is what bounds bias in the expectation",
            },
            "frozen_deployable_policy": dict(FROZEN_PILOT_POLICY),
            "primary_covq_field": "covq_deployable_fixed_pilot_exposure",
            "diagnostic_covq_field": "covq_oracle_angle_exposure",
            "accounting": {
                "covq_metric": "classical Fisher matrix of the emitted readout under the "
                               "frozen f=0.02 two-quadrature pilot policy; the "
                               "matched-analyzer value is retained only as a ceiling",
                "quest_metric": "mixed-state SLD QFIM of QUEST-tE (optimistic upper "
                                "bound, no readout compiled)",
                "sparse_metric": "mixed-state SLD QFIM (optimistic upper bound)",
                "two_qubit_noise": "charged per emitted two-qubit gate on every arm",
                "shot_metric": "shots to satisfy A^T F A >= G_req, i.e. 1 / floor margin",
            },
            "sparse_caratheodory_preparation": {
                "cx_count": sparse_cx,
                "caveat": "this prototype's sparse-amplitude preparation is a dense "
                          "upper bound, not a competitive routine; the arm's degradation "
                          "reflects that circuit, not an intrinsic property of "
                          "single-state realisation",
            },
            "quest_variant": "tE (published: insert, then joint L-BFGS reoptimisation)",
            "quest_preparation": {
                "stop_reason": qr.stop_reason,
                "quest_total_pauli_rotations": qr.depth_adaptive_length,
                "quest_two_qubit_pauli_rotations": qr.two_qubit_rotations,
                "quest_emitted_cx_count": quest_cx,
                "quest_two_qubit_depth": quest_depth,
                "note": "a two-qubit Pauli rotation lowers to two CX; rotations and "
                        "emitted gates are reported separately and never conflated",
            },
            "lowering_symmetry": {
                "native_gate_set": "{1q, CX} after covq.circuits.lower_to_cx, both arms",
                "near_zero_angle_removal": "applied to both (no-op for CovQ, which has "
                                           "no parameterised preparation rotations)",
                "adjacent_rotation_combination": "applied to both (no-op on the "
                                                 "registered solutions)",
                "routing": "all-to-all assumed for both; no SWAP inserted on either arm",
                "noise_placement": "one two-qubit depolarizing event per emitted CX on "
                                   "every arm, via the same simulate_circuit_noisy path",
            },
            "rows": rows, "status": "DONE"}


if __name__ == "__main__":
    for rel, builder in (
        ("results/noise/adaptive_readout_performance.json", adaptive_readout_performance),
        ("results/noise/fixed_setting_cost_solver.json", fixed_setting_cost_solver),
        ("results/noise/quest_operational_comparison.json", quest_operational_comparison),
    ):
        payload = builder()
        _write(rel, payload)
        print(f"{payload['status']:<6} {rel}")
