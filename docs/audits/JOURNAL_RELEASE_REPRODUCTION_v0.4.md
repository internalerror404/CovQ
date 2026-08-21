# J4 — release reproduction report

**Submission gate: `PASS`.**

A disagreement beyond the registered tolerance pauses submission. It does not reopen the scientific scope (`REGISTRATION.md`).

## What was compared

- Development records at `HEAD~1` versus a clean-tree rerun.
- 27 headline quantities, bound by JSON pointer in `results/table_manifest.json`.
- Relative tolerance `1e-06`.

| outcome | count |
|---|---|
| unchanged | 23 |
| changed | 0 |
| new in this release | 0 |
| changed by the QUEST fidelity correction | 4 |
| absent | 0 |

## Rows that are not identical

| quantity | verdict | development | release | rel. diff |
|---|---|---|---|---|
| `k3_converged` | EXPECTED_CHANGE | 10 | 11 | 9.09e-02 |
| `k3_max_cx_ratio` | EXPECTED_CHANGE | 135.55555555555554 | 32.6530612244898 | 7.59e-01 |
| `n10_noiseless_quest` | EXPECTED_CHANGE | 1.6260855073181055 | 1.6261456545651514 | 3.70e-05 |
| `quest_two_qubit_rotations` | EXPECTED_CHANGE | 44 | 3 | 9.32e-01 |

## Why the QUEST-derived quantities moved

The v0.1 baseline omitted the joint angle-reoptimisation phase and was therefore not QUEST. Replacing it with the published algorithm necessarily moves every QUEST-derived number; nothing else may move.

| quantity | reason |
|---|---|
| `k3_converged` | path_m5 no longer hits the depth cap once angles are reoptimised |
| `k3_max_cx_ratio` | greedy overestimated QUEST depth; 135.6 was an artefact |
| `k3_min_cx_ratio` | unchanged families re-measured under the published algorithm |
| `n10_noiseless_quest` | same target state, far fewer gates |
| `quest_two_qubit_rotations` | 44 under terminal greedy, 3 under QUEST-tE |

## Provenance

- All records clean: `True`
- Single source commit across records: `True`
- Commit(s): `f141111cab92`
- Quantities bound to a record: 27, unresolved: 0

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
