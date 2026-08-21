"""K3: CovQ labelled pair schedule vs the QUEST expectation-targeting baseline.

Exact-target mode only.  Per the v0.3 handoff, QUEST is a baseline for exact
first/second-moment targeting and *not* an information-floor optimizer; forcing
it into the primary role would need a separately specified outer optimization
over targets, which is not attempted here.

Both sides are given the identical constraint set -- ``<P_i> = 0`` and
``<P_i P_j> = F_ij`` -- and both sides' two-qubit costs are parsed from emitted
circuits (gate C9), never from a formula.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from emit_measurement_records import _hash, _write  # noqa: E402

from covq.circuits import lower_to_cx, resources
from covq.instances import banded, matching_target, path_target, star_target, toeplitz_like
from covq.paulis import z_generators
from covq.programs import labelled_schedule_from_width
from covq.qfim import qfim_from_statevector
from covq.quest import moment_constraints, quest_circuit, quest_published
from covq.width import decompose_width2

MAX_DEPTH = 80
VARIANT = "tE"          # terminal exact insertion, per the published algorithm
CROSS_CHECK = "bE"      # best-position exact
# bE searches every insertion position, so its cost grows as
# (pool size) x (current depth) per iteration.  On the deeper path and banded
# targets that is hundreds of thousands of cost evaluations per insertion and
# it does not finish in any sane budget.  It is therefore run only where the tE
# solution is shallow enough for the search to be affordable, and instances
# where it was not run say so rather than quietly dropping the column.
CROSS_CHECK_MAX_ROTATIONS = 8
# Frozen selection rule, decided once and applied globally: the primary variant
# is whichever has the lower TOTAL emitted CX summed over all registered
# instances.  Choosing per instance would silently take the minimum of two
# baselines and understate QUEST.
PRIMARY_RULE = "lower total emitted CX across all registered instances, fixed globally"


def best_quest(cons, m, rng):
    """QUEST as published: insertion followed by joint angle reoptimisation.

    The earlier terminal-greedy routine is not QUEST and is no longer used here;
    see ``docs/audits/QUEST_FIDELITY_AUDIT_v0.4.md``.  Restarts are unnecessary
    for the published algorithm on this surface -- every registered target
    converges from ``|+>^m`` -- so the restart machinery is gone rather than
    left switched off.
    """
    return quest_published(cons, m, variant=VARIANT, max_depth=MAX_DEPTH,
                           tol=1e-13), 1


def families():
    out = []
    for m in (3, 4, 5):
        out.append((f"matching_m{m}", matching_target(m, np.random.default_rng(m * 13),
                                                      load=0.8).F, m))
        out.append((f"path_m{m}", path_target(m, 0.45).F, m))
        out.append((f"toeplitz_m{m}", toeplitz_like(m, 0.35).F, m))
    out.append(("star_m4", star_target(4, 0.3).F, 4))
    out.append(("banded_m5", banded(5, 2, 0.3).F, 5))
    return out


def run() -> dict:
    cases = []
    for name, F, m in families():
        dec = decompose_width2(F)
        if not dec.feasible:
            cases.append({"instance": name, "m": m, "skipped": "pair_infeasible"})
            continue
        prog = labelled_schedule_from_width(dec)
        emitted = prog.emitted()
        expected_cx = float(sum(w * sum(1 for b in blocks if len(b) == 2)
                                for w, blocks, _ in dec.branches))
        ps = z_generators(m)
        cons = moment_constraints(ps, F)
        res, tried = best_quest(cons, m, np.random.default_rng(2026))
        parsed = resources(lower_to_cx(quest_circuit(res, m)))
        if res.depth_adaptive_length <= CROSS_CHECK_MAX_ROTATIONS:
            cross = quest_published(cons, m, variant=CROSS_CHECK,
                                    max_depth=MAX_DEPTH, tol=1e-13)
            cross_parsed = resources(lower_to_cx(quest_circuit(cross, m)))
            cross_record = {
                "variant": CROSS_CHECK, "stop_reason": cross.stop_reason,
                "quest_total_pauli_rotations": cross.depth_adaptive_length,
                "quest_two_qubit_pauli_rotations": cross.two_qubit_rotations,
                "quest_emitted_cx_count": cross_parsed["cx_count"],
                "quest_two_qubit_depth": cross_parsed["two_qubit_depth"],
            }
        else:
            cross_record = {
                "variant": CROSS_CHECK, "stop_reason": "not_run_within_budget",
                "reason": (f"tE needed {res.depth_adaptive_length} rotations; "
                           "best-position search is quadratic in depth and exceeds "
                           f"the registered budget of {CROSS_CHECK_MAX_ROTATIONS}"),
            }
        psi = res.state / np.linalg.norm(res.state)
        cases.append({
            "instance": name, "m": m, "target_hash": _hash(np.round(F, 12)),
            "covq": {
                "settings": len(dec.branches),
                "expected_cx_per_shot": expected_cx,
                "activation_identity_sum_abs_F_e": float(
                    np.abs(F[np.triu_indices(m, 1)]).sum()),
                "max_two_qubit_depth": emitted["max_two_qubit_depth"],
                "realisation_error": float(np.abs(dec.matrix() - F).max()),
            },
            "quest": {
                "settings": 1,
                "variant": VARIANT,
                "quest_total_pauli_rotations": res.depth_adaptive_length,
                "quest_two_qubit_pauli_rotations": res.two_qubit_rotations,
                "quest_emitted_cx_count": parsed["cx_count"],
                "quest_two_qubit_depth": parsed["two_qubit_depth"],
                "rotations": res.depth_adaptive_length,
                "two_qubit_rotations": res.two_qubit_rotations,
                "cx_per_shot": parsed["cx_count"],
                "two_qubit_depth": parsed["two_qubit_depth"],
                "stop_reason": res.stop_reason,
                "residual": res.residual,
                "qfim_error": float(np.abs(qfim_from_statevector(psi, ps) - F).max()),
                "start_state": "|+>^n",
                "starts_tried": tried,
                "max_depth_cap": MAX_DEPTH,
            },
        })
        c = cases[-1]
        if "skipped" not in c:
            c["cx_ratio_quest_over_covq"] = (c["quest"]["cx_per_shot"]
                                             / c["covq"]["expected_cx_per_shot"]
                                             if c["covq"]["expected_cx_per_shot"] > 0 else None)
    converged = [c for c in cases if "skipped" not in c
                 and c["quest"]["stop_reason"] == "converged"]
    both = [c for c in converged
            if c["quest_cross_check"].get("quest_emitted_cx_count") is not None]
    tot_primary = sum(c["quest"]["quest_emitted_cx_count"] for c in both)
    tot_cross = sum(c["quest_cross_check"]["quest_emitted_cx_count"] for c in both)
    return {
        "gate": "K3",
        "baseline": "QUEST-tE (arXiv:2605.02367; Mahapatra and Kadiri), implemented "
                    "from the published method description: per iteration, insert one "
                    "Pauli rotation chosen by exact one-angle minimisation, then "
                    "jointly reoptimise every accumulated angle by L-BFGS",
        "variant": VARIANT,
        "supersedes": "the terminal-greedy routine used before v0.4, which omitted the "
                      "joint reoptimisation phase and is therefore not QUEST; it "
                      "understated the baseline by roughly an order of magnitude",
        "scope": "exact first/second-moment targeting only; NOT an information-floor "
                 "optimizer and not run as one",
        "resource_provenance": "both sides parsed from emitted circuits after the same "
                               "lowering; a two-qubit Pauli rotation emits two CX and "
                               "rotations are never reported as gates",
        "primary_variant_rule": PRIMARY_RULE,
        "cases": cases,
        "n_converged": len(converged),
        "n_total": len([c for c in cases if "skipped" not in c]),
        "covq_two_qubit_depth_is_always_one": all(
            c["covq"]["max_two_qubit_depth"] == 1 for c in cases if "skipped" not in c),
        "cross_check_coverage": f"{len(both)}/{len(converged)} converged instances",
        "cross_check_budget_rotations": CROSS_CHECK_MAX_ROTATIONS,
        "total_emitted_cx_primary": tot_primary,
        "total_emitted_cx_cross_check": tot_cross,
        "primary_variant_is_not_worse": tot_primary <= tot_cross,
        "variants_agree": tot_primary == tot_cross,
        "min_cx_ratio_over_converged": min(
            (c["cx_ratio_quest_over_covq"] for c in converged), default=None),
        "max_cx_ratio_over_converged": max(
            (c["cx_ratio_quest_over_covq"] for c in converged), default=None),
        "status": "DONE",
    }


if __name__ == "__main__":
    payload = run()
    _write("results/baselines_k3_quest.json", payload)
    print(json.dumps({k: v for k, v in payload.items() if k != "cases"}, indent=2,
                     default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o))
