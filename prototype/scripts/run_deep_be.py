"""Run the three registered QUEST-bE cases omitted by the v0.5 budget.

The scientific configuration matches the canonical K3 campaign: |+>^m start,
full default Pauli pool, exact one-angle grid 64, joint L-BFGS maxiter 200,
max depth 80, and residual tolerance 1e-13. Each case is run in a separate CI
job so timeout or failure is explicit rather than silently dropped.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from covq.circuits import lower_to_cx, resources
from covq.instances import banded, path_target
from covq.paulis import z_generators
from covq.qfim import qfim_from_statevector
from covq.quest import moment_constraints, quest_circuit, quest_published


def target(name: str):
    if name == "path_m4":
        return path_target(4, 0.45).F, 4
    if name == "path_m5":
        return path_target(5, 0.45).F, 5
    if name == "banded_m5":
        return banded(5, 2, 0.3).F, 5
    raise ValueError(name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("instance", choices=["path_m4", "path_m5", "banded_m5"])
    ap.add_argument("--out-dir", default="results/quest_be_deep")
    args = ap.parse_args()

    F, m = target(args.instance)
    ps = z_generators(m)
    cons = moment_constraints(ps, F)
    started = time.time()
    perf = time.perf_counter()
    result = quest_published(
        cons,
        m,
        variant="bE",
        max_depth=80,
        tol=1e-13,
        grid=64,
        maxiter=200,
    )
    elapsed = time.perf_counter() - perf
    parsed = resources(lower_to_cx(quest_circuit(result, m)))
    psi = result.state / np.linalg.norm(result.state)
    payload = {
        "schema_version": "covq.quest_be_deep/0.1",
        "instance": args.instance,
        "variant": "bE",
        "configuration": {
            "start_state": "|+>^m",
            "pool": "covq.quest.default_pool (all 1q Paulis plus all 2q Paulis)",
            "grid": 64,
            "joint_optimizer": "L-BFGS-B",
            "joint_maxiter": 200,
            "max_depth": 80,
            "tolerance": 1e-13,
        },
        "started_unix": started,
        "elapsed_seconds": elapsed,
        "stop_reason": result.stop_reason,
        "converged": bool(result.converged),
        "residual": float(result.residual),
        "quest_total_pauli_rotations": result.depth_adaptive_length,
        "quest_two_qubit_pauli_rotations": result.two_qubit_rotations,
        "quest_emitted_cx_count": parsed["cx_count"],
        "quest_two_qubit_depth": parsed["two_qubit_depth"],
        "qfim_max_abs_error": float(np.abs(qfim_from_statevector(psi, ps) - F).max()),
        "history": [float(x) for x in result.history],
        "rotations": [
            {"pauli": rot.label(), "angle": float(angle)}
            for rot, angle in result.rotations
        ],
    }
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{args.instance}.json"
    path.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: payload[k] for k in (
        "instance", "elapsed_seconds", "stop_reason", "converged", "residual",
        "quest_total_pauli_rotations", "quest_two_qubit_pauli_rotations",
        "quest_emitted_cx_count", "quest_two_qubit_depth", "qfim_max_abs_error")}, indent=2))


if __name__ == "__main__":
    main()
