# CovQ v0.4 — QUEST Citation and Reimplementation Fidelity Audit

## Disposition

**Citation metadata: verified.**

- arXiv:2605.02367v1
- Anjali Mahapatra and Gururaj Kadiri
- *Quantum State Engineering Under Multiple Expectation-Value Constraints*
- QUEST expands to *Quantum Unitary Engineering of States to Target*.

**Benchmark fidelity: not yet verified; the current code is not a literal implementation of any complete published QUEST variant.**

This is a release blocker for manuscript wording such as “faithful QUEST reimplementation.” It does not affect the CovQ theorems or compiler results.

## What the paper specifies

The paper defines four QUEST variants:

| insertion location | operator selection | name |
|---|---|---|
| terminal | gradient | QUEST-tG |
| terminal | exact | QUEST-tE |
| best position | gradient | QUEST-bG |
| best position | exact | QUEST-bE |

Every published variant has two phases per iteration:

1. insert one Pauli rotation;
2. jointly reoptimize all accumulated angles using classical optimization.

The numerical study uses L-BFGS for that joint optimization.

For the exact insertion variants, the one-angle cost is represented as

\[
C(\theta)=\alpha+\beta\cos(2\theta)+\gamma\sin(2\theta)
          +\delta\cos(4\theta)+\epsilon\sin(4\theta),
\]

recovered from five cost samples (four new samples because the zero-angle cost is shared), and minimized over a dense grid. QUEST-bE searches all insertion positions; QUEST-tE searches only the terminal position.

## What the current CovQ baseline does

`prototype/src/covq/quest.py`:

- appends only at the terminal position;
- selects the rotation by a direct exact fringe fit;
- never jointly reoptimizes earlier angles;
- uses dense-grid plus golden-section refinement;
- defaults to \(|+\rangle^{\otimes n}\), with random restart handling outside the core routine;
- describes itself as optimizer-free and “faithful.”

The direct three-point expectation fit is mathematically capable of recovering the same one-angle trigonometric information in a statevector simulator. That is not the central problem. The load-bearing difference is the omitted full-angle optimization phase; the current implementation is therefore best described as a **QUEST-inspired terminal greedy Pauli-path baseline**, not QUEST-tE or QUEST-bE.

## Required correction

Choose one of two defensible paths.

### Path A — faithful benchmark rerun (preferred)

Implement and freeze at least QUEST-tE as published:

1. terminal exact insertion;
2. published one-angle cost reconstruction or an algebraically proven equivalent;
3. joint L-BFGS optimization of all accumulated angles after every insertion;
4. explicit operator pool, initial-state rule, restart budget, tolerance, and depth cap;
5. emitted-circuit resource parsing.

For a stronger comparator, also implement QUEST-bE and report the better published variant under a frozen compute budget. The paper reports bE as the fastest exact variant in one of its studies, so omitting best-position insertion must be disclosed.

Rerun only the QUEST-dependent scientific records, then perform the final clean full emission so every canonical artifact identifies the same release SHA.

### Path B — relabel and narrow

Retain the current numbers but rename the baseline everywhere to, for example:

> terminal greedy Pauli-rotation expectation-targeting baseline inspired by QUEST

Delete “faithful QUEST reimplementation,” “QUEST required,” and any unqualified method-level comparison. The results then characterize this implementation only, not QUEST.

Path B preserves the paper but weakens K3 materially. Path A is recommended.

## Manuscript language pending correction

Until Path A passes, replace or quarantine:

- “Against a faithful QUEST reimplementation …”
- “QUEST required factors of 2.5–135.6 …”
- “QUEST uses 44 two-qubit rotations …”
- “QUEST is more shot efficient at zero noise …”

Those numerical statements remain valid for the current emitted baseline circuit, but not yet for a verified implementation of the cited QUEST algorithm.

## Release rule

The current clean rerun may finish unchanged and should be preserved as a reproducibility checkpoint. Do not create the submission tag from it. Apply the benchmark-fidelity correction on the release branch, rerun the affected records, then run one final clean canonical emission before tagging.
