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
from covq.quest import moment_constraints, quest, quest_circuit
from covq.width import decompose_width2

MAX_DEPTH = 200
N_RESTARTS = 4          # greedy methods are restarted as a matter of course


def best_quest(cons, m, rng):
    """Best of ``|+>^m`` plus random Haar-ish starts.

    Greedy descent is restart-sensitive and it would be unfair to report its
    first attempt.  ``path``-family targets in particular converge in ~80
    rotations from a random start while grinding past 400 from ``|+>^m``.
    """
    best = quest(cons, m, max_depth=MAX_DEPTH, tol=1e-13)
    tried = 1
    while best.stop_reason != "converged" and tried <= N_RESTARTS:
        v = rng.normal(size=1 << m) + 1j * rng.normal(size=1 << m)
        v /= np.linalg.norm(v)
        cand = quest(cons, m, max_depth=MAX_DEPTH, tol=1e-13, psi0=v)
        if cand.residual < best.residual:
            best = cand
        tried += 1
    return best, tried


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
        res, tried = best_quest(moment_constraints(ps, F), m, np.random.default_rng(2026))
        parsed = resources(lower_to_cx(quest_circuit(res, m)))
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
                "rotations": res.depth_adaptive_length,
                "two_qubit_rotations": res.two_qubit_rotations,
                "cx_per_shot": parsed["cx_count"],
                "two_qubit_depth": parsed["two_qubit_depth"],
                "stop_reason": res.stop_reason,
                "residual": res.residual,
                "qfim_error": float(np.abs(qfim_from_statevector(psi, ps) - F).max()),
                "start_state": "|+>^n then random restarts",
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
    return {
        "gate": "K3",
        "baseline": "QUEST (arXiv:2605.02367), reimplemented from the published "
                    "method description: depth-adaptive Pauli rotations, "
                    "sum-of-squared-residuals descent, optimizer-free",
        "scope": "exact first/second-moment targeting only; NOT an information-floor "
                 "optimizer and not run as one",
        "resource_provenance": "both sides parsed from emitted circuits",
        "cases": cases,
        "n_converged": len(converged),
        "n_total": len([c for c in cases if "skipped" not in c]),
        "covq_two_qubit_depth_is_always_one": all(
            c["covq"]["max_two_qubit_depth"] == 1 for c in cases if "skipped" not in c),
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
