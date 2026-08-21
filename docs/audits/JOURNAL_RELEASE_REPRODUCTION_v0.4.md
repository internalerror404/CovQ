# J4 — release reproduction report

**Submission gate: `PASS`.**

A disagreement beyond the registered tolerance pauses submission. It does not reopen the scientific scope (`REGISTRATION.md`).

## What was compared

- Development records at `HEAD~1` versus a clean-tree rerun.
- 30 headline quantities, bound by JSON pointer in `results/table_manifest.json`.
- Relative tolerance `1e-06`.

| outcome | count |
|---|---|
| unchanged | 29 |
| changed | 0 |
| new in this release | 0 |
| changed by the QUEST fidelity correction | 1 |
| absent | 0 |

## Rows that are not identical

| quantity | verdict | development | release | rel. diff |
|---|---|---|---|---|
| `n10_noiseless_covq` | EXPECTED_CHANGE | 1.7854844662521132 | 1.7854705518721636 | 7.79e-06 |

## Why the QUEST-derived quantities moved

The v0.1 baseline omitted the joint angle-reoptimisation phase and was therefore not QUEST. Replacing it with the published algorithm necessarily moves every QUEST-derived number; nothing else may move.

| quantity | reason |
|---|---|
| `k3_converged` | path_m5 no longer hits the depth cap once angles are reoptimised |
| `k3_max_cx_ratio` | greedy overestimated QUEST depth; 135.6 was an artefact |
| `k3_min_cx_ratio` | unchanged families re-measured under the published algorithm |
| `n10_noiseless_quest` | same target state, far fewer gates |
| `quest_two_qubit_rotations` | 44 under terminal greedy, 3 under QUEST-tE |
| `n10_noiseless_covq` | deployable template switched from the sampled to the exact binomial channel; the two agree to 7.8e-6 here |

## Provenance

- All records clean: `True`
- Single source commit across records: `True`
- Commit(s): `9b7f97015344`
- Quantities bound to a record: 30, unresolved: 0

## Manuscript cell audit

Tables not matching: **1**. Prose claims not matching: **3**.

| table | verdict | differing cells |
|---|---|---|
| `quest_exact` | MISMATCH | 17 |
| `pilot` | MATCH | 0 |
| `noise_transition` | MATCH | 0 |
| `n10` | EXPECTED_CHANGE | 6 |
| `fixed_setting` | MATCH | 0 |

The one non-matching table is intended: the N10 table gains a deployable column, because the frozen `f = 0.02` fixed-pilot exposure becomes the primary CovQ figure and the matched-analyzer value is demoted to a diagnostic ceiling.

## Preserved negative results

These must survive into the manuscript unchanged.

| result | required status |
|---|---|
| fixed-setting greedy support selection | `FAILS`, 7.24 % worst excess |
| QUEST at zero two-qubit noise | more shot-efficient than the schedule |
| `shot_scaled_edge_compile` | `LIMITED` |
| `K1b` end-to-end compiler collision | `NOT_ESTABLISHED` |
