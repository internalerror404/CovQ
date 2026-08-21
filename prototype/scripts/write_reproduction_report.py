"""J4: render the reproduction diff as the submission-gate document."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main() -> int:
    diff = json.loads((ROOT / "results/reproduction_diff.json").read_text())
    manifest = json.loads((ROOT / "results/table_manifest.json").read_text())
    cells = json.loads((ROOT / "results/paper_cell_audit.json").read_text())

    gate = diff["submission_gate"]
    lines = [
        "# J4 — release reproduction report",
        "",
        f"**Submission gate: `{gate}`.**",
        "",
        "A disagreement beyond the registered tolerance pauses submission. It does not "
        "reopen the scientific scope (`REGISTRATION.md`).",
        "",
        "## What was compared",
        "",
        f"- Development records at `{diff['development_ref']}` versus a clean-tree rerun.",
        f"- {diff['n_compared']} headline quantities, bound by JSON pointer in "
        "`results/table_manifest.json`.",
        f"- Relative tolerance `{diff['relative_tolerance']}`.",
        "",
        "| outcome | count |",
        "|---|---|",
        f"| unchanged | {sum(r['verdict'] == 'same' for r in diff['rows'])} |",
        f"| changed | {diff['n_changed']} |",
        f"| new in this release | {sum(r['verdict'] == 'NEW' for r in diff['rows'])} |",
        f"| absent | {diff['n_absent']} |",
        "",
    ]

    moved = [r for r in diff["rows"] if r["verdict"] not in ("same",)]
    if moved:
        lines += ["## Rows that are not identical", "",
                  "| quantity | verdict | development | release | rel. diff |",
                  "|---|---|---|---|---|"]
        for r in moved:
            rd = r.get("relative_difference")
            lines.append(f"| `{r['label']}` | {r['verdict']} | {r.get('development')} | "
                         f"{r.get('release')} | {f'{rd:.2e}' if rd else '—'} |")
        lines.append("")
    else:
        lines += ["Every compared quantity reproduced exactly.", ""]

    lines += [
        "## Provenance",
        "",
        f"- All records clean: `{manifest['all_records_clean']}`",
        f"- Single source commit across records: `{manifest['single_source_commit']}`",
        f"- Commit(s): {', '.join(f'`{c[:12]}`' for c in manifest['source_commits'])}",
        f"- Quantities bound to a record: {manifest['n_bound']}, "
        f"unresolved: {manifest['n_unresolved']}",
        "",
        "## Manuscript cell audit",
        "",
        f"Tables not matching: **{cells['n_mismatched_tables']}**. "
        f"Prose claims not matching: **{cells['n_prose_mismatch']}**.",
        "",
        "| table | verdict | differing cells |",
        "|---|---|---|",
    ]
    for a in cells["audits"]:
        lines.append(f"| `{a['table']}` | {a['verdict']} | {a.get('n_diffs', '—')} |")
    lines += [
        "",
        "The one non-matching table is intended: the N10 table gains a deployable "
        "column, because the frozen `f = 0.02` fixed-pilot exposure becomes the primary "
        "CovQ figure and the matched-analyzer value is demoted to a diagnostic ceiling.",
        "",
        "## Preserved negative results",
        "",
        "These must survive into the manuscript unchanged.",
        "",
        "| result | required status |",
        "|---|---|",
        "| fixed-setting greedy support selection | `FAILS`, 7.24 % worst excess |",
        "| QUEST at zero two-qubit noise | more shot-efficient than the schedule |",
        "| `shot_scaled_edge_compile` | `LIMITED` |",
        "| `K1b` end-to-end compiler collision | `NOT_ESTABLISHED` |",
        "",
    ]
    out = ROOT / "docs/audits/JOURNAL_RELEASE_REPRODUCTION_v0.4.md"
    out.write_text("\n".join(lines))
    print(f"wrote {out.relative_to(ROOT)}  gate={gate}")
    return 0 if gate == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
