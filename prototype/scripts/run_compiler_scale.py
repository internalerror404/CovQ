#!/usr/bin/env python3
"""Compiler-only CovQ scale benchmark: cutting-plane SDP master + matching pricing.

No statevector or density-matrix simulation is used. The branch master is solved
as an LP over accumulated Loewner eigenvector cuts. The exact pricing subproblem
is a maximum-weight matching on the hardware graph.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import networkx as nx
import numpy as np
import scipy
from scipy.optimize import linprog


@dataclass(frozen=True)
class Branch:
    matching: tuple[tuple[int, int], ...]
    signs: tuple[int, ...]
    cost: float
    C: np.ndarray  # downstream CFI pullback A^T B A

    @property
    def key(self) -> tuple:
        return self.matching, self.signs


def grid_graph(n: int) -> nx.Graph:
    rows = max(2, int(math.floor(math.sqrt(n))))
    cols = int(math.ceil(n / rows))
    while rows * cols < n:
        cols += 1
    g0 = nx.grid_2d_graph(rows, cols)
    nodes = sorted(g0.nodes())[:n]
    g = g0.subgraph(nodes).copy()
    if not nx.is_connected(g):
        comp = max(nx.connected_components(g), key=len)
        g = g.subgraph(comp).copy()
    return nx.convert_node_labels_to_integers(g, ordering="sorted")


def heavy_hex_graph(n: int) -> nx.Graph:
    """Connected heavy-hex family obtained by subdividing a hexagonal lattice."""
    side = 2
    while True:
        h = nx.hexagonal_lattice_graph(side, side, periodic=False, with_positions=False)
        g = nx.Graph()
        base = {u: i for i, u in enumerate(sorted(h.nodes(), key=str))}
        g.add_nodes_from(base.values())
        nxt = len(base)
        for u, v in sorted(h.edges(), key=lambda e: (str(e[0]), str(e[1]))):
            w = nxt
            nxt += 1
            g.add_edge(base[u], w)
            g.add_edge(w, base[v])
        if g.number_of_nodes() >= n:
            break
        side += 1
    order = list(nx.bfs_tree(g, min(g.nodes())).nodes())
    keep = order[:n]
    sub = g.subgraph(keep).copy()
    if not nx.is_connected(sub):
        raise RuntimeError("BFS prefix should be connected")
    return nx.convert_node_labels_to_integers(sub, ordering="default")


def bfs_groups(g: nx.Graph, d: int) -> list[list[int]]:
    order = list(nx.bfs_tree(g, min(g.nodes())).nodes())
    groups = [list(x) for x in np.array_split(np.asarray(order, dtype=int), d)]
    if any(len(x) == 0 for x in groups):
        raise ValueError("too many downstream groups")
    return groups


def design_matrix(g: nx.Graph, d: int) -> tuple[np.ndarray, list[list[int]]]:
    groups = bfs_groups(g, d)
    A = np.zeros((g.number_of_nodes(), d), dtype=float)
    for k, group in enumerate(groups):
        A[group, k] = 1.0 / math.sqrt(len(group))
    return A, groups


def branch_from_matching(A: np.ndarray, matching: Iterable[tuple[int, int]], signs: Iterable[int], c0: float, edge_cost: float) -> Branch:
    matching = tuple(sorted(tuple(sorted((int(i), int(j)))) for i, j in matching))
    signs = tuple(int(s) for s in signs)
    C = np.asarray(A.T @ A, dtype=float).copy()
    for (i, j), s in zip(matching, signs):
        ai, aj = A[i], A[j]
        C += s * (np.outer(ai, aj) + np.outer(aj, ai))
    return Branch(matching, signs, c0 + edge_cost * len(matching), C)


def product_branch(A: np.ndarray, c0: float) -> Branch:
    return Branch((), (), c0, A.T @ A)


def pricing_oracle(g: nx.Graph, A: np.ndarray, Y: np.ndarray, c0: float, edge_cost: float) -> tuple[Branch, float, float]:
    AY = A @ Y
    weighted = nx.Graph()
    weighted.add_nodes_from(g.nodes())
    edge_sign: dict[tuple[int, int], int] = {}
    for i, j in g.edges():
        qij = float(AY[i] @ A[j])
        w = 2.0 * abs(qij) - edge_cost
        if w > 0.0:
            weighted.add_edge(i, j, weight=w)
            edge_sign[tuple(sorted((i, j)))] = 1 if qij >= 0 else -1
    matching_set = nx.algorithms.matching.max_weight_matching(weighted, maxcardinality=False, weight="weight")
    matching = tuple(sorted(tuple(sorted(e)) for e in matching_set))
    signs = tuple(edge_sign[e] for e in matching)
    match_weight = float(sum(weighted.edges[e]["weight"] for e in matching))
    base = float(np.trace(A @ Y @ A.T)) - c0
    return branch_from_matching(A, matching, signs, c0, edge_cost), base + match_weight, match_weight


def unique_unit_vectors(d: int) -> list[np.ndarray]:
    out = []
    for k in range(d):
        v = np.zeros(d)
        v[k] = 1.0
        out.append(v)
    out.append(np.ones(d) / math.sqrt(d))
    return out


def compile_floor(g: nx.Graph, A: np.ndarray, G_req: np.ndarray, c0: float = 1.0, edge_cost: float = 0.002, tol_slack: float = 2e-7, tol_price: float = 2e-7, max_iterations: int = 500) -> dict:
    d = A.shape[1]
    branches = [product_branch(A, c0)]
    branch_keys = {branches[0].key}
    cuts = unique_unit_vectors(d)
    oracle_calls = duplicate_pricing = lp_calls = 0
    t_start = time.perf_counter()
    final = None

    for iteration in range(1, max_iterations + 1):
        M = np.array([[float(v @ b.C @ v) for b in branches] for v in cuts])
        rhs = np.array([float(v @ G_req @ v) for v in cuts])
        res = linprog(np.array([b.cost for b in branches]), A_ub=-M, b_ub=-rhs, bounds=[(0.0, None)] * len(branches), method="highs", options={"presolve": True})
        lp_calls += 1
        if not res.success:
            return {"status": "master_infeasible", "message": res.message, "iterations": iteration, "wall_clock_s": time.perf_counter() - t_start}
        n = np.asarray(res.x, dtype=float)
        Fdown = sum(float(x) * b.C for x, b in zip(n, branches))
        slack = (Fdown - G_req + (Fdown - G_req).T) / 2.0
        eigvals, eigvecs = np.linalg.eigh(slack)
        min_slack = float(eigvals[0])

        eta = -np.asarray(res.ineqlin.marginals, dtype=float)
        Y = np.zeros((d, d), dtype=float)
        for weight, v in zip(eta, cuts):
            if weight > 0:
                Y += float(weight) * np.outer(v, v)
        Y = (Y + Y.T) / 2.0
        dual_bound = float(np.sum(eta * rhs))

        candidate, violation, match_weight = pricing_oracle(g, A, Y, c0, edge_cost)
        oracle_calls += 1
        added_branch = False
        if violation > tol_price:
            if candidate.key not in branch_keys:
                branches.append(candidate)
                branch_keys.add(candidate.key)
                added_branch = True
            else:
                duplicate_pricing += 1

        added_cut = False
        if min_slack < -tol_slack:
            v = eigvecs[:, 0]
            pivot = int(np.argmax(np.abs(v)))
            if v[pivot] < 0:
                v = -v
            if not any(abs(float(v @ old)) > 1.0 - 1e-9 for old in cuts):
                cuts.append(v.copy())
                added_cut = True

        if not added_branch and not added_cut:
            primal = float(res.fun)
            gap = primal - dual_bound
            final = {
                "status": "optimal", "iterations": iteration,
                "wall_clock_s": time.perf_counter() - t_start,
                "objective": primal, "dual_bound": dual_bound,
                "primal_dual_gap": gap,
                "relative_gap": gap / max(1.0, abs(primal)),
                "min_contract_slack": min_slack,
                "pricing_violation": float(violation),
                "matching_weight_last": match_weight,
                "column_count": len(branches),
                "active_column_count": int(np.count_nonzero(n > 1e-9)),
                "loewner_cut_count": len(cuts), "oracle_calls": oracle_calls,
                "lp_calls": lp_calls, "duplicate_pricing_events": duplicate_pricing,
                "total_exposure": float(n.sum()),
                "max_matching_size": max(len(b.matching) for b in branches),
                "mean_active_matching_size": float(sum(x * len(b.matching) for x, b in zip(n, branches)) / max(n.sum(), 1e-15)),
            }
            break

    if final is None:
        final = {"status": "iteration_limit", "iterations": max_iterations, "wall_clock_s": time.perf_counter() - t_start, "column_count": len(branches), "loewner_cut_count": len(cuts), "oracle_calls": oracle_calls, "lp_calls": lp_calls}
    return final


def run_case(topology: str, m: int, seed: int, d: int, gamma: float, c0: float, edge_cost: float) -> dict:
    _ = np.random.default_rng(seed)
    g = grid_graph(m) if topology == "grid" else heavy_hex_graph(m)
    A, groups = design_matrix(g, d)
    result = compile_floor(g, A, gamma * np.eye(d), c0=c0, edge_cost=edge_cost)
    return {
        "topology": topology,
        "topology_definition": "rectangular nearest-neighbor grid induced prefix" if topology == "grid" else "edge-subdivided hexagonal lattice (heavy-hex family), connected BFS prefix",
        "m": m, "n_edges": g.number_of_edges(), "max_degree": max(dict(g.degree()).values()),
        "is_bipartite": nx.is_bipartite(g), "seed": seed,
        "downstream_dimension": d, "group_sizes": [len(x) for x in groups],
        "gamma": gamma, "base_shot_cost": c0, "edge_activation_cost": edge_cost,
        **result,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="results/compiler_scale")
    ap.add_argument("--sizes", default="50,100,200,300,500")
    ap.add_argument("--seeds", default="2026,3407,9181")
    ap.add_argument("--d", type=int, default=8)
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--c0", type=float, default=1.0)
    ap.add_argument("--edge-cost", type=float, default=0.002)
    args = ap.parse_args()

    sizes = [int(x) for x in args.sizes.split(",") if x]
    seeds = [int(x) for x in args.seeds.split(",") if x]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for topology in ("grid", "heavy_hex"):
        for m in sizes:
            for seed in seeds:
                print(f"running {topology} m={m} seed={seed}", flush=True)
                rec = run_case(topology, m, seed, args.d, args.gamma, args.c0, args.edge_cost)
                records.append(rec)
                print(json.dumps({k: rec.get(k) for k in ("status", "wall_clock_s", "column_count", "oracle_calls", "min_contract_slack", "relative_gap")}), flush=True)

    payload = {
        "schema_version": "covq.compiler_scale/0.1",
        "scientific_scope": "compiler-only; no state or channel simulation",
        "algorithm": "Loewner eigenvector cutting planes + LP restricted master + exact maximum-weight-matching pricing",
        "environment": {
            "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
            "networkx": nx.__version__, "platform": platform.platform(),
            "processor": platform.processor(), "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        },
        "registration": {"sizes": sizes, "topologies": ["grid", "heavy_hex"], "seeds": seeds, "downstream_dimension": args.d, "gamma": args.gamma, "base_shot_cost": args.c0, "edge_activation_cost": args.edge_cost},
        "records": records,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["record_hash"] = hashlib.sha256(canonical.encode()).hexdigest()
    (out_dir / "compiler_scale.json").write_text(json.dumps(payload, indent=2) + "\n")

    fields = sorted({k for r in records for k in r})
    with (out_dir / "compiler_scale.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in r.items()})

    summary = []
    for topology in ("grid", "heavy_hex"):
        for m in sizes:
            rows = [r for r in records if r["topology"] == topology and r["m"] == m]
            ok = [r for r in rows if r["status"] == "optimal"]
            rec = {"topology": topology, "m": m, "n_runs": len(rows), "n_optimal": len(ok)}
            for key in ("wall_clock_s", "column_count", "oracle_calls", "loewner_cut_count", "iterations", "relative_gap", "min_contract_slack"):
                vals = np.array([r[key] for r in ok], dtype=float) if ok else np.array([])
                if vals.size:
                    rec[f"{key}_median"] = float(np.median(vals))
                    rec[f"{key}_min"] = float(vals.min())
                    rec[f"{key}_max"] = float(vals.max())
            summary.append(rec)
    (out_dir / "compiler_scale_summary.json").write_text(json.dumps({"summary": summary}, indent=2) + "\n")
    print(f"wrote {out_dir}")


if __name__ == "__main__":
    main()
