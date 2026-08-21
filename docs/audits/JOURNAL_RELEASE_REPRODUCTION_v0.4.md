# J4 — release reproduction report

**Submission gate: `PASS`.**

A disagreement beyond the registered tolerance pauses submission. It does not reopen the scientific scope (`REGISTRATION.md`).

## What was compared

- Development records at `HEAD` versus a clean-tree rerun.
- 27 headline quantities, bound by JSON pointer in `results/table_manifest.json`.
- Relative tolerance `1e-06`.

| outcome | count |
|---|---|
| unchanged | 27 |
| changed | 0 |
| new in this release | 0 |
| absent | 0 |

Every compared quantity reproduced exactly.

## Provenance

- All records clean: `True`
- Single source commit across records: `True`
- Commit(s): `dc2213544ffd`
- Quantities bound to a record: 27, unresolved: 0

## Manuscript cell audit

Tables not matching: **0**. Prose claims not matching: **0**.

| table | verdict | differing cells |
|---|---|---|
| `quest_exact` | MATCH | 0 |
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
