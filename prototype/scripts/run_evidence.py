#!/usr/bin/env python3
"""Evidence run for the Task 0A/0B audit.

This is NOT Task 0C v0.1, which is SUPERSEDED (see ../STATUS.md). It exists to
reproduce the numbers quoted in docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md.

Usage (from the ``covq`` directory):

    python3 scripts/run_evidence.py --out results/prototype_evidence.json

Everything it writes is derived from emitted circuits and exact linear algebra.
It does not tune any tolerance, and it stops at the first gate failure.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from covq import baselines as bl  # noqa: E402
from covq import gates as gt  # noqa: E402
from covq import polytope as pol  # noqa: E402
from covq import programs as prg  # noqa: E402
from covq import width as wid  # noqa: E402
from covq.instances import PROTOCOL_SEEDS, task0_suite  # noqa: E402
from covq.paulis import z_generators  # noqa: E402
from covq.qfim import crb_conditioning  # noqa: E402
from covq.topology import build as build_topology  # noqa: E402


def environment() -> dict:
    import scipy

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "seeds": PROTOCOL_SEEDS,
        "qfim_normalisation": "per shot",
        "resource_provenance": "emitted circuits, lowered to {1q, CX} and routed",
    }


def compile_instance(inst, max_m_milp: int = 5, max_m_dense: int = 6) -> dict:
    m = inst.m
    ps = z_generators(m)
    out: dict = {
        "instance": inst.name,
        "family": inst.family,
        "m": m,
        "ground_truth": {k: v for k, v in inst.ground_truth.items()
                         if k != "hidden_decomposition"},
        "meta": inst.meta,
        "target_conditioning": crb_conditioning(inst.F),
    }

    t0 = time.perf_counter()
    dec = pol.exact_decompose(inst.F)
    out["exact_lp"] = {
        "feasible": bool(dec.feasible),
        "support": dec.support,
        "residual": dec.residual if np.isfinite(dec.residual) else None,
        "seconds": time.perf_counter() - t0,
        "certificate": None if dec.certificate is None else {
            "kind": dec.certificate.kind,
            "margin": dec.certificate.margin,
            "verified": bool(dec.certificate.verify(inst.F)),
            "data": dec.certificate.data if dec.certificate.kind == "hypermetric" else "omitted",
        },
    }

    t0 = time.perf_counter()
    d2 = wid.decompose_width2(inst.F)
    out["width2"] = {
        "feasible": bool(d2.feasible),
        "n_branches": d2.n_branches,
        "residual": d2.residual if np.isfinite(d2.residual) else None,
        "seconds": time.perf_counter() - t0,
        "violation": None if d2.violation is None else {
            "kind": d2.violation.kind,
            "subset": list(d2.violation.subset),
            "lhs": d2.violation.lhs,
            "rhs": d2.violation.rhs,
        },
    }
    if m <= 6:
        t0 = time.perf_counter()
        brute = wid.decompose_width_k(inst.F, 2, force_enumeration=True)
        out["width2"]["brute_force_agrees"] = bool(brute.feasible) == bool(d2.feasible)
        out["width2"]["brute_force_seconds"] = time.perf_counter() - t0
    out["min_width"] = wid.min_width(inst.F, max_m_exact=6)

    if not dec.feasible:
        return out

    red = pol.caratheodory_reduce(dec)
    out["caratheodory"] = {"support": red.support,
                           "bound": m * (m - 1) // 2 + 1,
                           "residual": red.residual}

    t0 = time.perf_counter()
    cg = pol.column_generation(inst.F, rng=np.random.default_rng(PROTOCOL_SEEDS[0]))
    out["column_generation"] = {"feasible": bool(cg.feasible), "support": cg.support,
                                "seconds": time.perf_counter() - t0, **cg.info}

    t0 = time.perf_counter()
    fw = pol.frank_wolfe(inst.F, iters=200)
    out["frank_wolfe"] = {"residual": fw.residual, "seconds": time.perf_counter() - t0,
                          **fw.info}

    if m <= max_m_milp:
        t0 = time.perf_counter()
        ms = pol.min_support_decompose(inst.F, max_m=max_m_milp)
        out["min_support_milp"] = {"feasible": bool(ms.feasible), "support": ms.support,
                                   "seconds": time.perf_counter() - t0}

    programs: dict[str, prg.Program] = {
        "labelled_global_cat": prg.labelled_schedule_from_signs(red.signs, red.weights),
        "caratheodory_single_pure_state": prg.single_pure_state_program(red.signs, red.weights),
        "coherent_flag_relaxed": prg.coherent_flag_program(red.signs, red.weights),
        "coherent_flag_literal": prg.coherent_flag_program_literal(red.signs, red.weights),
    }
    if d2.feasible:
        programs["labelled_width2"] = prg.labelled_schedule_from_width(d2)
    try:
        programs["baseline_dense_state_prep"] = bl.dense_state_preparation(
            red.signs, red.weights, max_m=max_m_dense)
    except ValueError:
        pass
    try:
        programs["baseline_fixed_state"] = bl.fixed_state_cat_backend(inst.F)
    except ValueError:
        pass
    programs["baseline_global_ghz"] = bl.global_ghz(m)
    programs["baseline_product_probe"] = bl.product_probe(m)
    had = bl.hadamard_sign_schedule(m)
    if had is not None:
        programs["baseline_hadamard_schedule"] = had
    rand_prog, rand_resid = bl.random_signed_cat_schedule(
        inst.F, max(2, red.support), np.random.default_rng(PROTOCOL_SEEDS[1]))
    programs["baseline_random_cat_schedule"] = rand_prog

    rows = {}
    for name, prog in programs.items():
        F = prog.realised_qfim(ps)
        entry = {
            "qfim_frobenius_error": float(np.linalg.norm(F - inst.F, "fro")),
            "qfim_operator_error": float(np.linalg.norm(F - inst.F, 2)),
            "realises_target": bool(np.linalg.norm(F - inst.F, "fro") < 1e-9),
            **prog.emitted(),
        }
        entry.pop("per_setting", None)
        if prog.kind == "labelled":
            audit = prg.audit_branch_width(prog)
            entry["max_declared_width"] = audit["max_declared_width"]
            entry["max_finest_width"] = audit["max_finest_width"]
            entry["width_declaration_holds"] = audit["all_declared_widths_hold"]
        if prog.metadata.get("gate_counts_are_upper_bound") or \
                prog.metadata.get("deliberately_weak"):
            entry["comparison_caveat"] = (
                "generic state preparation; count is an upper bound on the "
                "baseline's true cost, so any win against it is inadmissible")
        rows[name] = entry
    rows["baseline_random_cat_schedule"]["fit_residual"] = rand_resid
    out["programs"] = rows

    topo_rows = {}
    for topo in ("line", "square_grid", "heavy_hex", "modular_two_cluster"):
        edges = build_topology(topo, m)
        for name in ("labelled_global_cat", "labelled_width2"):
            if name not in programs:
                continue
            em = programs[name].emitted(edges=edges)
            topo_rows[f"{name}@{topo}"] = {
                "two_qubit_gate_count": em["two_qubit_gate_count"],
                "swap_count": em["swap_count"],
                "max_two_qubit_depth": em["max_two_qubit_depth"],
            }
    out["routing"] = topo_rows
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/prototype_evidence.json")
    ap.add_argument("--m", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--skip-gates", action="store_true")
    args = ap.parse_args()

    instances = task0_suite(tuple(args.m))
    report: dict = {"environment": environment(),
                    "n_instances": len(instances),
                    "m_values": args.m}

    if not args.skip_gates:
        print("running correctness gates ...", flush=True)
        t0 = time.perf_counter()
        report["gates"] = gt.run_all(instances, stop_on_failure=True)
        report["gates"]["seconds"] = time.perf_counter() - t0
        for g in report["gates"]["gates"]:
            print(f"  {g['status']:9s} {g['name']}")
        if report["gates"]["summary"]["n_fail"]:
            out = Path(args.out)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(report, indent=2, default=str))
            print(f"\nGATE FAILURE -- report written to {out} and aborting "
                  f"without tolerance changes.", file=sys.stderr)
            return 1

    print("compiling instances ...", flush=True)
    rows = []
    for inst in instances:
        t0 = time.perf_counter()
        rows.append(compile_instance(inst))
        print(f"  {inst.name:28s} {time.perf_counter() - t0:6.2f}s", flush=True)
    report["instances"] = rows

    feasible = [r for r in rows if r["exact_lp"]["feasible"]]
    w2 = [r for r in feasible if r["width2"]["feasible"]]
    disagreements = [r["instance"] for r in rows
                     if r["width2"].get("brute_force_agrees") is False]
    report["summary"] = {
        "n_feasible": len(feasible),
        "n_infeasible": len(rows) - len(feasible),
        "n_width2_feasible": len(w2),
        "width2_bruteforce_disagreements": disagreements,
        "min_width_histogram": _hist([r["min_width"].get("k") for r in feasible]),
        "caratheodory_support_max_ratio": max(
            (r["caratheodory"]["support"] / r["caratheodory"]["bound"]
             for r in feasible if "caratheodory" in r), default=None),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"\nwrote {out}")
    print(json.dumps(report["summary"], indent=2, default=str))
    return 0


def _hist(vals) -> dict:
    h: dict[str, int] = {}
    for v in vals:
        h[str(v)] = h.get(str(v), 0) + 1
    return h


if __name__ == "__main__":
    raise SystemExit(main())
