"""J4: compare the clean release rerun against the development records.

Development records live in git at the pre-rerun commit; the rerun overwrites
the working tree.  Comparing the two answers the only question that matters for
submission: does a clean tree reproduce the numbers the manuscript quotes?

A disagreement beyond the registered tolerance **pauses submission**.  It does
not reopen the scientific scope.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from emit_table_manifest import BINDINGS, resolve  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
REL_TOL = 1e-6
ABS_FLOOR = 1e-9

# Fields introduced by the J2 deployable-pilot work.  They cannot exist in the
# development records, so their absence there is expected and must not be
# reported as a reproduction failure -- but it also must not be waved through
# silently, so they are named individually rather than pattern-matched.
EXPECTED_NEW = {
    "n10_noiseless_covq",
    "n10_primary_field",
    # Bindings added with the emitted-gate accounting fix and the bE cross-check.
    "quest_emitted_cx",
    "k3_bE_coverage",
    "k3_variants_agree",
}

# Pointers whose *field* was renamed, so the development record cannot be
# resolved even though the quantity is unchanged.  Recorded separately from
# genuinely new bindings so the distinction survives in the artifact.
EXPECTED_RENAMED = {
    "quest_two_qubit_rotations":
        "two_qubit_rotations -> quest_two_qubit_pauli_rotations, to stop rotations "
        "being read as emitted gates",
}

# Quantities the QUEST fidelity correction was *supposed* to move.  The v0.1
# baseline omitted the joint angle-reoptimisation phase and so was not QUEST;
# replacing it with the published algorithm necessarily changes every
# QUEST-derived number.  Naming them individually, with the reason, keeps the
# gate meaningful: anything else that moves is still a reproduction failure.
EXPECTED_CHANGED = {
    "k3_converged": "path_m5 no longer hits the depth cap once angles are reoptimised",
    "k3_max_cx_ratio": "greedy overestimated QUEST depth; 135.6 was an artefact",
    "k3_min_cx_ratio": "unchanged families re-measured under the published algorithm",
    "n10_noiseless_quest": "same target state, far fewer gates",
    "quest_two_qubit_rotations": "44 under terminal greedy, 3 under QUEST-tE",
    "n10_noiseless_covq":
        "deployable template switched from the sampled to the exact binomial "
        "channel; the two agree to 7.8e-6 here",
}


def git_show(ref: str, rel: str):
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "show", f"{ref}:{rel}"],
                             capture_output=True, text=True, check=True).stdout
        return json.loads(out)
    except Exception:
        return None


def compare(a, b):
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None:
        return ("same" if a == b else "CHANGED"), None
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        d = abs(a - b)
        rel = d / max(abs(a), abs(b), ABS_FLOOR)
        return ("same" if rel <= REL_TOL else "CHANGED"), rel
    if isinstance(a, list) and isinstance(b, list):
        return ("same" if a == b else "CHANGED"), None
    return ("same" if a == b else "CHANGED"), None


def main(ref: str = "HEAD") -> int:
    rows, changed, absent = [], 0, 0
    for label, rel, pointer in BINDINGS:
        dev_doc = git_show(ref, rel)
        cur_path = ROOT / rel
        cur_doc = json.loads(cur_path.read_text()) if cur_path.exists() else None
        expected_new = label in EXPECTED_NEW or label in EXPECTED_RENAMED
        if dev_doc is None or cur_doc is None:
            rows.append({"label": label, "artifact": rel,
                         "verdict": "NEW" if expected_new else "ABSENT",
                         "development": None, "release": None})
            absent += 0 if expected_new else 1
            continue
        try:
            cur = resolve(cur_doc, pointer)
        except Exception as exc:
            rows.append({"label": label, "artifact": rel, "verdict": "ABSENT",
                         "reason": f"missing in release: {exc}"})
            absent += 1
            continue
        try:
            dev = resolve(dev_doc, pointer)
        except Exception as exc:
            rows.append({"label": label, "artifact": rel,
                         "verdict": ("RENAMED" if label in EXPECTED_RENAMED
                                     else "NEW" if expected_new else "ABSENT"),
                         "development": None, "release": cur,
                         "reason": EXPECTED_RENAMED.get(
                             label, None if expected_new else str(exc))})
            absent += 0 if expected_new else 1
            continue
        verdict, rel_diff = compare(dev, cur)
        row = {"label": label, "artifact": rel, "json_pointer": pointer,
               "development": dev, "release": cur,
               "relative_difference": rel_diff, "verdict": verdict}
        if verdict == "CHANGED" and label in EXPECTED_CHANGED:
            row["verdict"] = "EXPECTED_CHANGE"
            row["reason"] = EXPECTED_CHANGED[label]
        else:
            changed += verdict == "CHANGED"
        rows.append(row)
    payload = {
        "schema_version": "covq.reproduction_diff/0.4",
        "development_ref": ref,
        "relative_tolerance": REL_TOL,
        "n_compared": len(rows), "n_changed": changed, "n_absent": absent,
        "n_expected_change": sum(r["verdict"] == "EXPECTED_CHANGE" for r in rows),
        "expected_change_reasons": EXPECTED_CHANGED,
        "expected_renamed": EXPECTED_RENAMED,
        "n_renamed": sum(r["verdict"] == "RENAMED" for r in rows),
        "rows": rows,
        "submission_gate": "PASS" if changed == 0 and absent == 0 else "PAUSE",
    }
    (ROOT / "results/reproduction_diff.json").write_text(
        json.dumps(payload, indent=2) + "\n")
    print(f"gate={payload['submission_gate']} compared={len(rows)} "
          f"changed={changed} absent={absent}")
    for r in rows:
        if r["verdict"] != "same":
            print(f"   {r['verdict']:<8} {r['label']}: "
                  f"{r.get('development')!r} -> {r.get('release')!r}")
    return 0 if payload["submission_gate"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "HEAD"))
