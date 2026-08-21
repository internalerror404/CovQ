"""J6: map every headline manuscript quantity to its artifact, pointer and hash.

The manuscript LaTeX source is not present in this repository, so the paper
cannot be regenerated here.  What *can* be built -- and is what a deterministic
paper-generation script consumes -- is the binding between each named quantity
and the record that produced it.  Any table cell whose value disagrees with the
value resolved here is, by construction, hand-transcribed.

Each entry resolves a JSON pointer against a canonical record and carries that
record's content hash and source commit, so a cell can be traced to a specific
run of a specific tree.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from emit_measurement_records import _write  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]

# (label, record, JSON pointer)
BINDINGS = [
    ("readout_worst_attainment_gap",
     "results/measurements/ideal_block_readout.json", "/worst_attainment_gap"),
    ("schedule_worst_cfi_gap",
     "results/measurements/schedule_attainability.json", "/worst_schedule_gap"),
    ("fixed_readout_information_at_origin",
     "results/measurements/readout_singularities.json", "/cases/0/fixed_all_x_trace"),
    ("matched_readout_information_at_origin",
     "results/measurements/readout_singularities.json", "/cases/0/matched_readout_trace"),
    ("dephasing_covariance_max_drift",
     "results/noise/qfi_covariance_separation.json", "/dephasing_covariance_max_drift"),
    ("dephasing_cfi_minus_qfim_max",
     "results/noise/qfi_covariance_separation.json", "/dephasing_cfi_minus_qfim_max"),
    ("estimator_status",
     "results/measurements/estimator_efficiency.json", "/status"),
    ("adaptive_status",
     "results/measurements/adaptive_recentering.json", "/status"),
    ("k3_converged", "results/baselines_k3_quest.json", "/n_converged"),
    ("k3_total", "results/baselines_k3_quest.json", "/n_total"),
    ("k3_min_cx_ratio", "results/baselines_k3_quest.json",
     "/min_cx_ratio_over_converged"),
    ("k3_max_cx_ratio", "results/baselines_k3_quest.json",
     "/max_cx_ratio_over_converged"),
    ("k3_depth_always_one", "results/baselines_k3_quest.json",
     "/covq_two_qubit_depth_is_always_one"),
    ("noisy_pricing_worst_gap", "results/noise/noisy_pricing_oracle.json", "/worst_gap"),
    ("noisy_pricing_cases", "results/noise/noisy_pricing_oracle.json", "/n_cases"),
    ("floor_worst_optimality_gap",
     "results/noise/noise_aware_floor_compiler.json", "/worst_optimality_gap"),
    ("entanglement_threshold",
     "results/noise/noise_aware_floor_compiler.json", "/entanglement_threshold_between"),
    ("worst_setting_crossover_1e4",
     "results/noise/setting_cost_phase_diagram.json", "/worst_crossover_at_c_setup_1e4"),
    ("greedy_worst_excess_percent",
     "results/noise/fixed_setting_cost_solver.json", "/worst_greedy_excess_percent"),
    ("greedy_status", "results/noise/fixed_setting_cost_solver.json", "/status"),
    ("n7_clean_best_ratio",
     "results/noise/adaptive_readout_performance.json", "/arms/0/best_ratio"),
    ("n7_heavy_best_ratio",
     "results/noise/adaptive_readout_performance.json", "/arms/2/best_ratio"),
    ("n10_noiseless_covq",
     "results/noise/quest_operational_comparison.json",
     "/rows/0/covq_deployable_fixed_pilot_exposure"),
    ("n10_noiseless_quest",
     "results/noise/quest_operational_comparison.json",
     "/rows/0/quest_shots_noisy_qfi_upper_bound"),
    ("n10_primary_field",
     "results/noise/quest_operational_comparison.json", "/primary_covq_field"),
    ("quest_two_qubit_rotations",
     "results/noise/quest_operational_comparison.json",
     "/quest_preparation/two_qubit_rotations"),
    ("sparse_preparation_cx",
     "results/noise/quest_operational_comparison.json",
     "/sparse_caratheodory_preparation/cx_count"),
]


def resolve(doc, pointer: str):
    node = doc
    for raw in pointer.lstrip("/").split("/"):
        key = raw.replace("~1", "/").replace("~0", "~")
        node = node[int(key)] if isinstance(node, list) else node[key]
    return node


def build() -> dict:
    entries, missing = [], []
    for label, rel, pointer in BINDINGS:
        path = ROOT / rel
        if not path.exists():
            missing.append({"label": label, "artifact": rel, "reason": "artifact absent"})
            continue
        doc = json.loads(path.read_text())
        try:
            value = resolve(doc, pointer)
        except (KeyError, IndexError, ValueError) as exc:
            missing.append({"label": label, "artifact": rel, "pointer": pointer,
                            "reason": f"{type(exc).__name__}: {exc}"})
            continue
        entries.append({"label": label, "artifact": rel, "json_pointer": pointer,
                        "value": value, "record_hash": doc.get("record_hash"),
                        "record_source_commit": doc.get("source_commit"),
                        "record_dirty": doc.get("dirty")})
    commits = {e["record_source_commit"] for e in entries}
    return {
        "purpose": "binding between manuscript quantities and the records that produced "
                   "them; a cell disagreeing with a resolved value is hand-transcribed",
        "manuscript_source_present": False,
        "manuscript_note": "LaTeX source is not in this repository, so the paper cannot "
                           "be regenerated here; this manifest is the input a "
                           "deterministic generator would consume",
        "entries": entries, "unresolved": missing,
        "n_bound": len(entries), "n_unresolved": len(missing),
        "all_records_clean": all(e["record_dirty"] is False for e in entries),
        "single_source_commit": len(commits) == 1,
        "source_commits": sorted(c for c in commits if c),
        "status": "DONE" if not missing else "FAILS",
    }


if __name__ == "__main__":
    payload = build()
    _write("results/table_manifest.json", payload)
    print(f"{payload['status']}  bound={payload['n_bound']} "
          f"unresolved={payload['n_unresolved']} "
          f"clean={payload['all_records_clean']} "
          f"single_commit={payload['single_source_commit']}")
    for u in payload["unresolved"]:
        print("   MISSING:", u)
