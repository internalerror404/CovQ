# K3 v0.2 — QUEST as published, and what that costs the comparison

Supersedes the QUEST half of `K3_QUEST_BASELINE_AND_ADAPTIVE_READOUT_v0.1.md`.
Cause: `QUEST_FIDELITY_AUDIT_v0.4.md`. Evidence: `prototype/src/covq/quest.py`
(`quest_published`), `results/baselines_k3_quest.json`,
`results/noise/quest_operational_comparison.json`.

## What was wrong

Published QUEST (arXiv:2605.02367; Mahapatra and Kadiri) runs **two phases per
iteration**: insert one Pauli rotation, then **jointly reoptimise every accumulated
angle** by classical optimisation — L-BFGS in the paper's study. The implementation
behind v0.1 did only the first. It appended at the terminal position and never revisited
an earlier angle.

That is not a detail. Later rotations routinely make earlier angles suboptimal, so
without reoptimisation the depth needed to reach a target is a large overestimate. On the
registered `m = 4` contract target the difference is **50 rotations / 44 two-qubit gates
versus 4 / 3**.

Everything in v0.1 that quoted a QUEST resource count is therefore withdrawn. The
implementation is retained as `quest`, relabelled a *terminal greedy Pauli-path
expectation-targeting baseline*, and is never again reported as QUEST.

`quest_published` implements the algorithm with both phases, `variant='tE'` (terminal
exact) and `variant='bE'` (best-position exact). The joint phase uses an adjoint
gradient — one forward and one backward sweep regardless of how many observables or
angles — which is what makes reoptimising at every insertion affordable. The gradient is
checked against finite differences in the test suite.

## Corrected exact-target comparison

| instance | CovQ settings | CovQ E[cx] | QUEST rot | QUEST cx | QUEST 2q depth | cx ratio | v0.1 ratio |
|---|---|---|---|---|---|---|---|
| `matching_m3` | 2 | 0.800 | 1 | 2 | 2 | 2.50 | 2.50 |
| `path_m3` | 3 | 0.900 | 2 | 4 | 4 | **4.44** | ~~135.56~~ |
| `toeplitz_m3` | 4 | 0.823 | 2 | 4 | 4 | 4.86 | 4.86 |
| `matching_m4` | 2 | 1.600 | 2 | 4 | 2 | 2.50 | 2.50 |
| `path_m4` | 3 | 1.350 | 14 | 24 | 20 | **17.78** | ~~117.04~~ |
| `toeplitz_m4` | 5 | 1.338 | 3 | 6 | 6 | 4.48 | 4.48 |
| `matching_m5` | 2 | 1.600 | 2 | 4 | 2 | 2.50 | 2.50 |
| `path_m5` | 3 | 1.800 | 28 | 46 | 36 | **25.56** | ~~197.78, capped~~ |
| `toeplitz_m5` | 9 | 1.868 | 4 | 8 | 8 | 4.28 | 4.28 |
| `star_m4` | 4 | 0.900 | 5 | 8 | 8 | **8.89** | ~~124.44~~ |
| `banded_m5` | 5 | 1.470 | 26 | 48 | 36 | **32.65** | ~~80.27~~ |

**11 of 11 converge**, up from 10 — `path_m5` no longer hits the depth cap, so the
"one capped instance" caveat is also withdrawn. The ratio range is **2.50 – 32.65**,
not 2.5 – 197.8.

The matching and Toeplitz families are unchanged: QUEST already found short circuits
there, so joint reoptimisation had nothing to recover. Everything the greedy routine
found *hard* was an artefact of the greedy routine.

**What survives.** Two regimes still exist, and CovQ's two-qubit depth is still **1 in
every instance** against QUEST's 2–36 — that is the matching theorem, not a measurement,
and it does not depend on the baseline's quality. What does not survive is the
order-of-magnitude framing: the separation is one order, not two, and the largest ratio
is 32.65×.

## Corrected equal-accounting noise campaign

QUEST's preparation for the N10 target drops from 44 two-qubit gates to **3**, which is
most of the campaign's mechanism. CovQ is scored on its emitted-readout CFI under the
frozen `f = 0.02` pilot; QUEST gets its noisy QFI as an optimistic bound.

| edge depol | CovQ deployable | CovQ oracle | product | QUEST bound | winner | v0.1 QUEST |
|---|---|---|---|---|---|---|
| 0.00 | 1.785 | 1.766 | 3.255 | **1.626** | QUEST | 1.626 |
| 0.01 | 1.822 | 1.802 | 3.255 | **1.697** | QUEST | ~~2.876~~ |
| 0.02 | 1.860 | 1.839 | 3.255 | **1.769** | QUEST | ~~5.122~~ |
| 0.05 | **1.980** | 1.957 | 3.255 | 2.001 | CovQ | ~~32.80~~ |
| 0.10 | **2.207** | 2.180 | 3.255 | 2.469 | CovQ | ~~835.8~~ |
| 0.20 | **2.799** | 2.759 | 3.255 | 3.893 | CovQ | ~~351 089~~ |

**The crossover moves from 1 % to between 2 % and 5 %, and the advantage at 20 % collapses
from `1.27e5×` to `1.39×`.**

The qualitative regime statement survives: one-state preparation wins while coherence is
cheap; the depth-one schedule wins once it is not. The quantitative claim does not, and
any sentence quoting `3.51e5` must go.

Two things still make the surviving CovQ wins conservative rather than generous: CovQ is
measured on an attained classical Fisher matrix while QUEST is credited with a quantum
bound it has no compiled readout to reach, and CovQ pays its pilot while QUEST pays no
calibration at all.

## Manuscript sentences that must change

| location | current | required |
|---|---|---|
| abstract | "faithful QUEST reimplementation" | QUEST-tE as published, with the joint reoptimisation phase |
| abstract | "factors of $2.5$–$135.6$" | `2.5`–`32.7` |
| abstract | "ten converged instances… one harder path instance reached the QUEST depth cap" | eleven converged; no cap reached |
| abstract | "at $20\%$… 2.759 versus $3.51\times10^5$" | `2.799` versus `3.893` |
| abstract / §results | "less exposure at every tested point from $1\%$… onward" | from `5 %` onward; QUEST wins at `0`, `1 %`, `2 %` |
| §results-n10 | "44 two-qubit rotations in the QUEST preparation" | `3` |
| §results-n10 | "QUEST/CovQ exposure ratio is approximately $1.27\times10^5$" | `1.39` |
| §results-quest-exact | "On the ten converged targets the ratio ranges from $2.5$ to $135.6$" | eleven targets, `2.5` to `32.7` |
| §limitations | "faithful reimplementation, not official QUEST source" | keep, but say which variant (tE) and that bE agrees |

## What this does not touch

No CovQ theorem, certificate, compiler result, readout result, estimator result, or noise
theorem depends on the QUEST baseline. The pair-width characterisation, the cost identity,
the conic dual, the block-local pricing theorem, attainability, the estimator, the pinching
proposition and the fixed-setting negative result are all unaffected.
