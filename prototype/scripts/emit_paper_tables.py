"""J6: generate the numeric table bodies from records, and audit the manuscript.

Two jobs.

*Generate* every numeric result table as an ``\\input``-able fragment under
``paper/generated/``, computed from the canonical records.  A table built this
way cannot drift from the evidence.

*Audit* the vendored manuscript source against those fragments, cell by cell.
Any cell that disagrees is, by definition, hand-transcribed or stale, and the
audit names it.  This is the check that "no numerical result remains manually
transcribed" is actually true rather than intended.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from emit_measurement_records import _write  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "paper/generated"
TEX = ROOT / "paper/CovQ_Paper_v0.4.tex"


def load(rel):
    return json.loads((ROOT / rel).read_text())


def _f(x, nd):
    return f"{x:.{nd}f}"


def tab_quest_exact(k3):
    rows = []
    for c in k3["cases"]:
        if "skipped" in c:
            continue
        q, v = c["quest"], c["covq"]
        status = ("converged" if q["stop_reason"] == "converged"
                  else "depth cap; feasible target")
        fam, _, size = c["instance"].rpartition("_")
        label = f"{fam}-${size}$"
        rows.append(f"{label} & {v['settings']} & "
                    f"{_f(v['expected_cx_per_shot'], 3)} & {q['cx_per_shot']} & "
                    f"{q['two_qubit_depth']} & "
                    f"{_f(c['cx_ratio_quest_over_covq'], 2)} & {status}\\\\")
    return rows


def tab_pilot(n7):
    arms = {a["noise_level"]: {r["pilot_fraction"]: r["deployable_over_oracle"]
                               for r in a["sweep"]} for a in n7["arms"]}
    best = {a["noise_level"]: a["optimal_pilot_fraction"] for a in n7["arms"]}
    rows = []
    for f in sorted(arms["clean"]):
        cells = []
        for tag in ("clean", "light", "heavy"):
            v = _f(arms[tag][f], 4)
            cells.append(f"\\textbf{{{v}}}" if best[tag] == f and tag != "clean" else v)
        rows.append(f"{f:.3f} & " + " & ".join(cells) + "\\\\")
    return rows


def tab_noise_transition(n6):
    rows = []
    for s in n6["sweep"]:
        if s["edge_depolarizing"] not in (0.0, 0.05, 0.10, 0.15, 0.20):
            continue
        rows.append(f"{_f(s['edge_depolarizing'], 2)} & {_f(s['primal_cost'], 4)} & "
                    f"{s['settings']} & {s['entangled_settings']}\\\\")
    return rows


def tab_n10(n10):
    rows = []
    for r in n10["rows"]:
        sparse = r["sparse_caratheodory_shots_noisy_qfi_upper_bound"]
        quest = r["quest_shots_noisy_qfi_upper_bound"]
        primary = r["covq_deployable_fixed_pilot_exposure"]
        rows.append(
            f"{_f(r['edge_depolarizing'], 2)} & {_f(primary, 3)} & "
            f"{_f(r['covq_oracle_angle_exposure'], 3)} & "
            f"{_f(r['product_only_shots'], 3)} & "
            f"{_f(quest, 3) if quest else '{infinite}'} & "
            f"{_f(sparse, 3) if sparse else '{infinite}'}\\\\")
    return rows


def tab_fixed_setting(n9):
    rows = []
    for r in n9["rows"]:
        exc = r["greedy_excess_percent"]
        cell = "0" if abs(exc) < 1e-9 else f"{exc:.2f}\\%"
        rows.append(f"{r['lambda_q']:.1f} & {_f(r['greedy_cost'], 5)} & "
                    f"{r['greedy_settings']} & {_f(r['exhaustive_cost'], 5)} & "
                    f"{r['exhaustive_settings']} & {cell}\\\\")
    return rows


NUM = re.compile(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?")

_WORDS = {8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}


def audit(name: str, generated: list[str], tex: str) -> dict:
    """Locate the table in the manuscript by its first data row and diff cells."""
    if not generated:
        return {"table": name, "verdict": "NO_ROWS"}
    key = generated[0].split("&")[0].strip()
    idx = tex.find("\n" + key + " &")
    if idx < 0:
        return {"table": name, "verdict": "NOT_FOUND_IN_TEX", "anchor": key}
    block = tex[idx: tex.find("\\bottomrule", idx)]
    tex_rows = [ln.strip() for ln in block.strip().splitlines()
                if ln.strip().endswith("\\\\")]
    diffs = []
    for i, gen in enumerate(generated):
        if i >= len(tex_rows):
            diffs.append({"row": i, "issue": "missing in manuscript", "generated": gen})
            continue
        g_nums = NUM.findall(gen)
        t_nums = NUM.findall(tex_rows[i])
        if len(g_nums) != len(t_nums):
            diffs.append({"row": i, "issue": "cell count differs",
                          "generated": gen, "manuscript": tex_rows[i]})
            continue
        for j, (a, b) in enumerate(zip(g_nums, t_nums)):
            fa, fb = float(a), float(b)
            denom = max(abs(fa), abs(fb), 1e-9)
            if abs(fa - fb) / denom > 5e-3:
                diffs.append({"row": i, "cell": j, "generated": fa, "manuscript": fb,
                              "relative_difference": abs(fa - fb) / denom})
    return {"table": name, "n_rows": len(generated), "n_diffs": len(diffs),
            "diffs": diffs, "verdict": "MATCH" if not diffs else "MISMATCH"}


# Numbers quoted in prose -- abstract, contribution list, discussion -- are not
# in any table and so are invisible to the cell audit, yet they are exactly the
# figures a reader remembers.  Each entry is (label, literal as it appears in the
# source, record, JSON pointer, formatter).
PROSE_CLAIMS = [
    ("noisy_pricing_cases", "12", "results/noise/noisy_pricing_oracle.json",
     "/n_cases", lambda v: str(int(v))),
    ("greedy_excess", "7.24", "results/noise/fixed_setting_cost_solver.json",
     "/worst_greedy_excess_percent", lambda v: f"{v:.2f}"),
    ("k3_min_ratio", "2.5", "results/baselines_k3_quest.json",
     "/min_cx_ratio_over_converged", lambda v: f"{v:.1f}"),
    ("k3_max_ratio", "135.6", "results/baselines_k3_quest.json",
     "/max_cx_ratio_over_converged", lambda v: f"{v:.1f}"),
    ("k3_converged", "ten converged", "results/baselines_k3_quest.json",
     "/n_converged", lambda v: f"{_WORDS.get(int(v), str(int(v)))} converged"),
    ("quest_two_qubit", "44 two-qubit rotations",
     "results/noise/quest_operational_comparison.json",
     "/quest_preparation/two_qubit_rotations", lambda v: f"{int(v)} two-qubit rotations"),
    ("sparse_cx", "232-CX", "results/noise/quest_operational_comparison.json",
     "/sparse_caratheodory_preparation/cx_count", lambda v: f"{int(v)}-CX"),
]


def audit_prose(tex: str) -> list[dict]:
    from emit_table_manifest import resolve

    out = []
    for label, literal, rel, pointer, fmt in PROSE_CLAIMS:
        path = ROOT / rel
        if not path.exists():
            out.append({"claim": label, "verdict": "ARTIFACT_ABSENT"})
            continue
        expected = fmt(resolve(json.loads(path.read_text()), pointer))
        out.append({
            "claim": label, "quoted_in_source": literal, "from_record": expected,
            "artifact": rel, "json_pointer": pointer,
            "verdict": ("MATCH" if literal == expected else "MISMATCH"),
            "present_in_source": literal in tex,
        })
    return out


# Claims withdrawn by the QUEST fidelity correction.  They must not survive in
# the active manuscript.  Archived retraction documents are exempt -- the record
# of what was withdrawn has to keep quoting the withdrawn numbers.
RETRACTED_LITERALS = [
    "135.6", "197.8", "80--198", "80-198", "ten converged",
    "44 two-qubit rotations", "five orders of magnitude",
    r"3.51\times10^5", "1.27\times10^5", "351088", "351{,}088",
]


def audit_retracted(tex: str) -> list[dict]:
    out = []
    for literal in RETRACTED_LITERALS:
        hits = tex.count(literal)
        out.append({"literal": literal, "occurrences": hits,
                    "verdict": "CLEAR" if hits == 0 else "PRESENT_MUST_REMOVE"})
    return out


def main() -> int:
    GEN.mkdir(parents=True, exist_ok=True)
    # The N10 table intentionally gains a column: the frozen f=0.02 deployable
    # exposure becomes the primary CovQ figure and the matched-analyzer value is
    # demoted to a diagnostic ceiling.  That is a required change, so it must not
    # be reported in the same voice as an accidental drift.
    EXPECTED_CHANGE = {"n10"}
    tables = {
        "quest_exact": tab_quest_exact(load("results/baselines_k3_quest.json")),
        "pilot": tab_pilot(load("results/noise/adaptive_readout_performance.json")),
        "noise_transition": tab_noise_transition(
            load("results/noise/noise_aware_floor_compiler.json")),
        "n10": tab_n10(load("results/noise/quest_operational_comparison.json")),
        "fixed_setting": tab_fixed_setting(
            load("results/noise/fixed_setting_cost_solver.json")),
    }
    for name, rows in tables.items():
        (GEN / f"tab_{name}.tex").write_text("\n".join(rows) + "\n")

    tex = TEX.read_text() if TEX.exists() else ""
    prose = audit_prose(tex) if tex else []
    retracted = audit_retracted(tex) if tex else []
    audits = []
    if tex:
        for n, r in tables.items():
            a = audit(n, r, tex)
            if n in EXPECTED_CHANGE and a["verdict"] == "MISMATCH":
                a["verdict"] = "EXPECTED_CHANGE"
                a["reason"] = ("deployable fixed-pilot column added as the primary "
                               "CovQ figure; oracle retained as a diagnostic ceiling")
            audits.append(a)
    payload = {
        "generated_fragments": sorted(str(p.relative_to(ROOT))
                                      for p in GEN.glob("tab_*.tex")),
        "manuscript_source": str(TEX.relative_to(ROOT)) if tex else None,
        "cell_relative_tolerance": 5e-3,
        "audits": audits,
        "n_mismatched_tables": sum(a["verdict"] not in ("MATCH", "EXPECTED_CHANGE")
                                   for a in audits),
        "prose_claims": prose,
        "retracted_literals": retracted,
        "n_retracted_present": sum(r["verdict"] != "CLEAR" for r in retracted),
        "retracted_scope": "active manuscript only; archived retraction documents "
                           "must keep quoting the withdrawn numbers",
        "n_prose_mismatch": sum(p["verdict"] != "MATCH" for p in prose),
        "n_prose_absent_from_source": sum(not p.get("present_in_source", True)
                                          for p in prose),
        "note": "the N10 table gains a deployable column: the frozen f=0.02 "
                "fixed-pilot exposure is the primary CovQ figure and the "
                "matched-analyzer value is retained as a diagnostic ceiling, so a "
                "cell-count difference there is expected and intended",
        "status": "DONE",
    }
    payload["status"] = ("DONE" if payload["n_mismatched_tables"] == 0
                         and payload["n_prose_mismatch"] == 0
                         and payload["n_retracted_present"] == 0 else "FAILS")
    _write("results/paper_cell_audit.json", payload)
    for a in audits:
        print(f"  {a['verdict']:<20} {a['table']:<18} diffs={a.get('n_diffs', '-')}")
    for rc in retracted:
        if rc["verdict"] != "CLEAR":
            print(f"  {rc['verdict']:<20} retracted:{rc['literal']!r} "
                  f"x{rc['occurrences']}")
    for pc in prose:
        flag = "" if pc.get("present_in_source", True) else "  (not found in source)"
        print(f"  {pc['verdict']:<20} prose:{pc['claim']:<13} "
              f"source={pc.get('quoted_in_source')!r} record={pc.get('from_record')!r}{flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
