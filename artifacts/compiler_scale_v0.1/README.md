# CovQ Compiler Scale Artifact v0.1

This artifact records the compiler-only scaling campaign added in manuscript v0.10.

## Registered surface

- sizes: `m = {50, 100, 200, 300, 500}`
- topologies: rectangular nearest-neighbor grid and an edge-subdivided hexagonal-lattice heavy-hex family
- downstream dimension: `d = 8`
- contract: eight disjoint normalized group modes with `G_req = 1.4 I_8`
- cost model: base shot `c0 = 1`, active edge `ce = 0.002`
- repeats: three single-threaded timing repeats (`2026`, `3407`, `9181` are repeat identifiers)
- excluded: statevector simulation, density-matrix simulation, device execution, routing, and correlated noise

## Algorithm

The solver alternates:

1. an LP restricted master over currently generated branch columns;
2. minimum-eigenvector cuts for the downstream Loewner constraint; and
3. exact maximum-weight-matching pricing on the hardware graph.

Termination requires both primal contract feasibility and no positive reduced-cost branch. The final dual matrix is the positive combination of accumulated eigenvector cuts, so the last pricing call certifies the full branch library rather than only the restricted master.

## Result

All 30 runs terminated with optimal status. At `m=500`, median wall clocks were `2.393 s` (grid, 40 columns/oracle calls) and `1.665 s` (heavy-hex, 35 columns/oracle calls). The largest recorded relative primal-dual gap magnitude was below `3.3e-14` and the largest contract-slack residual magnitude was below `5.6e-14`.

No asymptotic exponent is inferred from five sizes; the nonmonotone column counts reflect the finite graph boundaries and their commensurability with the eight group modes.

## Files

- `compiler_scale.json`: canonical full record
- `compiler_scale_summary.json`: per-topology/per-size median and range
- `compiler_scale.csv`: flat table
- implementation: `prototype/scripts/run_compiler_scale.py`
