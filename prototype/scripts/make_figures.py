#!/usr/bin/env python3
"""Deterministically regenerate CovQ journal Figures 3--4 from canonical JSON.

This script is intentionally read-only with respect to the result records.  It
writes figures plus a machine-readable manifest that maps each plotted series to
its source JSON pointer(s).  By default it requires release-clean records:
``dirty == false`` and one common ``source_commit`` across all inputs.

Expected inputs under REPO_ROOT:
  results/noise/adaptive_readout_performance.json
  results/noise/quest_operational_comparison.json
  results/noise/noise_aware_floor_compiler.json
  results/noise/fixed_setting_cost_solver.json
  results/baselines_k3_quest.json

Expected release field in N10:
  covq_deployable_fixed_pilot_exposure
The old development field can be used only with --allow-development-fallback.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INPUTS = {
    "n7": Path("results/noise/adaptive_readout_performance.json"),
    "n10": Path("results/noise/quest_operational_comparison.json"),
    "n6": Path("results/noise/noise_aware_floor_compiler.json"),
    "n9": Path("results/noise/fixed_setting_cost_solver.json"),
    "k3": Path("results/baselines_k3_quest.json"),
}


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise SystemExit(f"missing required input: {path}") from exc
    except json.JSONDecodeError as exc:
        raise SystemExit(f"invalid JSON in {path}: {exc}") from exc


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _finite_or_nan(value: Any) -> float:
    if value is None:
        return float("nan")
    out = float(value)
    return out if math.isfinite(out) else float("nan")


def _release_checks(records: dict[str, dict[str, Any]], require_clean: bool) -> str | None:
    commits = {str(record.get("source_commit")) for record in records.values()}
    if "None" in commits:
        raise SystemExit("at least one record lacks source_commit")
    if len(commits) != 1:
        raise SystemExit(f"records do not share one source_commit: {sorted(commits)}")
    commit = next(iter(commits))
    if require_clean:
        dirty = [name for name, record in records.items() if record.get("dirty") is not False]
        if dirty:
            raise SystemExit(f"release-clean mode rejected dirty/non-false records: {dirty}")
    return commit


def _save(fig: plt.Figure, out_base: Path) -> list[Path]:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "Title": out_base.stem,
        "Author": "CovQ release pipeline",
        "Creator": "CovQ_make_figures_v0.4.py",
        "Producer": "Matplotlib",
        "CreationDate": None,
        "ModDate": None,
    }
    pdf = out_base.with_suffix(".pdf")
    png = out_base.with_suffix(".png")
    fig.savefig(pdf, bbox_inches="tight", metadata=metadata)
    fig.savefig(png, bbox_inches="tight", dpi=220, metadata={"Software": "CovQ release pipeline"})
    plt.close(fig)
    return [pdf, png]


def _n10_covq_value(row: dict[str, Any], allow_fallback: bool) -> float:
    primary = "covq_deployable_fixed_pilot_exposure"
    if primary in row:
        return float(row[primary])
    if allow_fallback and "covq_shots_emitted_readout_cfi" in row:
        return float(row["covq_shots_emitted_readout_cfi"])
    raise SystemExit(
        "N10 lacks covq_deployable_fixed_pilot_exposure. "
        "Finish J2/J3, or use --allow-development-fallback only for a nonrelease dry run."
    )


def make_operational(
    records: dict[str, dict[str, Any]],
    out_dir: Path,
    allow_fallback: bool,
) -> tuple[list[Path], dict[str, Any]]:
    n7, n10, n6, k3 = records["n7"], records["n10"], records["n6"], records["k3"]

    fig, axes = plt.subplots(2, 3, figsize=(19.2, 10.6), constrained_layout=True)

    # (a) Pilot-recentered readout.
    ax = axes[0, 0]
    n7_ptrs: list[str] = []
    for arm_idx, arm in enumerate(n7["arms"]):
        x = [float(row["pilot_fraction"]) for row in arm["sweep"]]
        y = [float(row["deployable_over_oracle"]) for row in arm["sweep"]]
        ax.plot(x, y, marker="o", linewidth=1.8, label=str(arm["noise_level"]))
        n7_ptrs.extend(f"/arms/{arm_idx}/sweep/{i}" for i in range(len(x)))
    ax.axvline(0.02, linestyle=":", linewidth=1.4, color="0.45")
    ax.set_xscale("log")
    ax.set_ylim(0.80, 1.01)
    ax.set_xlabel("pilot fraction")
    ax.set_ylabel("deployable / oracle information")
    ax.set_title("(a) Pilot-recentered readout")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)

    # (b) Equal-accounting noisy contract campaign.
    ax = axes[0, 1]
    rows = n10["rows"]
    qe = np.array([float(r["edge_depolarizing"]) for r in rows])
    covq = np.array([_n10_covq_value(r, allow_fallback) for r in rows])
    product = np.array([_finite_or_nan(r.get("product_only_shots")) for r in rows])
    quest = np.array([_finite_or_nan(r.get("quest_shots_noisy_qfi_upper_bound")) for r in rows])
    sparse = np.array([_finite_or_nan(r.get("sparse_caratheodory_shots_noisy_qfi_upper_bound")) for r in rows])
    ax.plot(qe, covq, marker="o", label="CovQ deployable CFI")
    ax.plot(qe, product, marker="s", label="product")
    ax.plot(qe, quest, marker="^", label="QUEST noisy QFI bound")
    ax.plot(qe, sparse, marker="d", label="generic sparse-state QFI bound")
    ax.set_yscale("log")
    ax.set_xlabel("two-qubit depolarizing probability")
    ax.set_ylabel("exposure required by contract")
    ax.set_title("(b) Equal-accounting noise campaign")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)

    # (c) Noise-driven compiler transition.
    ax = axes[1, 0]
    sweep = n6["sweep"]
    x = np.array([float(r["edge_depolarizing"]) for r in sweep])
    cost = np.array([float(r["primal_cost"]) for r in sweep])
    settings = np.array([int(r["entangled_settings"]) for r in sweep])
    line1 = ax.plot(x, cost, marker="o", color="C0", label="certified cost")[0]
    ax.set_xlabel("edge depolarizing probability")
    ax.set_ylabel("certified cost", color=line1.get_color())
    ax.tick_params(axis="y", labelcolor=line1.get_color())
    ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    line2 = ax2.step(x, settings, where="post", color="C1",
                     label="entangled settings")[0]
    ax2.set_ylabel("entangled settings", color=line2.get_color())
    ax2.tick_params(axis="y", labelcolor=line2.get_color())
    ax.set_title("(c) Noise-driven compiler transition")

    # (e) Exposure ratio.  Essential: on the log-scale absolute plot the
    # crossover is invisible, and the crossover is the result.
    ax = axes[0, 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = quest / covq
    ax.plot(qe, ratio, marker="o", color="C3")
    ax.axhline(1.0, linewidth=1.2, color="0.35")
    ax.set_xlabel("two-qubit depolarizing probability")
    ax.set_ylabel("QUEST QFI ceiling / CovQ deployable CFI")
    ax.set_title("(e) Exposure ratio")
    ax.grid(True, alpha=0.25)

    # (d) Exact-target realization cost.
    ax = axes[1, 1]
    kept: list[tuple[str, float, str, int]] = []
    for idx, case in enumerate(k3["cases"]):
        if "skipped" in case:
            continue
        quest_info = case.get("quest", {})
        if quest_info.get("stop_reason") != "converged":
            continue
        name = str(case["instance"])
        family = name.split("_", 1)[0]
        ratio = float(case["cx_ratio_quest_over_covq"])
        kept.append((name, ratio, family, idx))
    labels = [name.replace("matching_", "match_").replace("toeplitz_", "toepl_")
              .replace("banded_", "band_") for name, _, _, _ in kept]
    ratios = [ratio for _, ratio, _, _ in kept]
    families = [family for _, _, family, _ in kept]
    family_order = {family: i for i, family in enumerate(dict.fromkeys(families))}
    colors = [f"C{family_order[f] % 10}" for f in families]
    ax.bar(np.arange(len(ratios)), ratios, color=colors)
    ax.axhline(1.0, linewidth=1.0, color="0.35")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=48, ha="right", fontsize=8)
    ax.set_ylabel("QUEST CX / CovQ expected CX")
    ax.set_title("(d) Exact-target realization cost")
    ax.grid(True, axis="y", alpha=0.25)

    axes[1, 2].axis("off")
    outputs = _save(fig, out_dir / "operational_results")
    manifest = {
        "figure": "operational_results",
        "panels": {
            "a": {"artifact": INPUTS["n7"].as_posix(), "json_pointers": n7_ptrs},
            "e": {"artifact": INPUTS["n10"].as_posix(),
                  "series": "quest_shots_noisy_qfi_upper_bound / "
                            "covq_deployable_fixed_pilot_exposure",
                  "reading": "points below one favour the optimistic QUEST QFI ceiling; "
                             "points above one establish that CovQ's emitted-readout CFI "
                             "beats that ceiling",
                  "json_pointers": [f"/rows/{i}" for i in range(len(rows))]},
            "b": {
                "artifact": INPUTS["n10"].as_posix(),
                "json_pointers": [f"/rows/{i}" for i in range(len(rows))],
                "covq_primary_field": (
                    "covq_deployable_fixed_pilot_exposure"
                    if all("covq_deployable_fixed_pilot_exposure" in r for r in rows)
                    else "covq_shots_emitted_readout_cfi (development fallback)"
                ),
            },
            "c": {"artifact": INPUTS["n6"].as_posix(), "json_pointers": [f"/sweep/{i}" for i in range(len(sweep))]},
            "d": {"artifact": INPUTS["k3"].as_posix(), "json_pointers": [f"/cases/{idx}" for *_, idx in kept]},
        },
    }
    return outputs, manifest


def make_fixed_setting(records: dict[str, dict[str, Any]], out_dir: Path) -> tuple[list[Path], dict[str, Any]]:
    n9 = records["n9"]
    rows = n9["rows"]
    x = np.array([float(r["lambda_q"]) for r in rows])
    greedy = np.array([float(r["greedy_cost"]) for r in rows])
    exhaustive = np.array([float(r["exhaustive_cost"]) for r in rows])

    fig, ax = plt.subplots(figsize=(9.2, 5.4), constrained_layout=True)
    ax.plot(x, exhaustive, marker="o", linewidth=2, label="exhaustive optimum")
    ax.plot(x, greedy, marker="s", linewidth=2, label="greedy forward")
    ax.fill_between(x, exhaustive, greedy, where=greedy >= exhaustive, alpha=0.16)
    ax.set_xlabel(r"fixed setting penalty $\lambda_q$")
    ax.set_ylabel("total objective")
    ax.set_title("Greedy fixed-setting selection fails in the multi-setting regime")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)

    outputs = _save(fig, out_dir / "fixed_setting_failure")
    manifest = {
        "figure": "fixed_setting_failure",
        "artifact": INPUTS["n9"].as_posix(),
        "json_pointers": [f"/rows/{i}" for i in range(len(rows))],
    }
    return outputs, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--out-dir", type=Path, default=Path("paper/figures"))
    parser.add_argument("--manifest", type=Path, default=Path("paper/generated/figure_manifest.json"))
    parser.add_argument("--allow-development-fallback", action="store_true")
    parser.add_argument("--no-require-clean", action="store_true")
    args = parser.parse_args()

    root = args.repo_root.resolve()
    records = {name: _read_json(root / rel) for name, rel in INPUTS.items()}
    commit = _release_checks(records, require_clean=not args.no_require_clean)

    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    manifest_path = args.manifest if args.manifest.is_absolute() else root / args.manifest

    op_outputs, op_manifest = make_operational(records, out_dir, args.allow_development_fallback)
    fs_outputs, fs_manifest = make_fixed_setting(records, out_dir)

    all_outputs = op_outputs + fs_outputs
    manifest = {
        "schema_version": "covq.figure_manifest/0.4",
        "source_commit": commit,
        "generator": Path(__file__).name,
        "inputs": {
            name: {
                "path": rel.as_posix(),
                "sha256": _sha256(root / rel),
                "record_hash": records[name].get("record_hash"),
            }
            for name, rel in INPUTS.items()
        },
        "figures": [op_manifest, fs_manifest],
        "outputs": [
            {"path": str(path.relative_to(root) if path.is_relative_to(root) else path), "sha256": _sha256(path)}
            for path in all_outputs
        ],
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"source_commit={commit}")
    for path in all_outputs:
        print(f"wrote {path} sha256={_sha256(path)}")
    print(f"wrote {manifest_path} sha256={_sha256(manifest_path)}")


if __name__ == "__main__":
    main()
