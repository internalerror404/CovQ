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

## Emitted-gate accounting, fixed in both directions

Three unit errors were found and corrected; one of them favoured CovQ and two favoured
QUEST.

1. **Rotations are not gates.** A two-qubit Pauli rotation lowers to **two** CX. The
   `m = 4` QUEST state uses 3 two-qubit *rotations* and **6 emitted CX**; describing it as
   "3 two-qubit gates" was wrong. Rotations and emitted gates are now separate record
   fields and are never conflated.
2. **The emitted circuit did not prepare its own input.** `quest_circuit` emitted only the
   rotations, omitting `|+⟩^m`, so simulating it started from `|0…0⟩` and produced a
   different state from the one the result claims. CX count unaffected — the preparation
   is single-qubit — but an artifact that does not prepare its own input is not an emitted
   program.
3. **Noise was charged per rotation, not per gate.** QUEST paid one two-qubit
   depolarizing event per Pauli rotation while CovQ paid one per Bell pair. Since a
   rotation is two CX, **QUEST was under-billed by exactly 2×**. Both arms now run through
   one code path that charges once per emitted CX.

Lowering symmetry is now explicit and identical on both arms: same native set `{1q, CX}`,
near-zero-angle removal, adjacent-rotation combination, all-to-all routing with no SWAPs,
one noise event per emitted CX. The first two rules are no-ops on the registered
solutions — no angle falls below tolerance and no two adjacent rotations share a Pauli —
which is worth being able to state rather than leaving silently absent.

## Corrected equal-accounting noise campaign

| edge depol | CovQ deployable | CovQ oracle | product | QUEST ceiling | ratio | reading |
|---|---|---|---|---|---|---|
| 0.00 | 1.785 | 1.766 | 3.255 | 1.626 | 0.911 | optimistic QUEST ceiling lower |
| 0.01 | 1.822 | 1.802 | 3.255 | 1.768 | 0.971 | optimistic QUEST ceiling lower |
| **0.02** | 1.860 | 1.839 | 3.255 | 1.917 | **1.031** | CovQ beats the ceiling |
| 0.05 | 1.980 | 1.957 | 3.255 | 2.441 | 1.233 | CovQ beats the ceiling |
| 0.10 | 2.207 | 2.180 | 3.255 | 3.699 | 1.676 | CovQ beats the ceiling |
| 0.20 | 2.799 | 2.759 | 3.255 | 8.961 | **3.201** | CovQ beats the ceiling |

Crossover on the registered grid: `q_edge ∈ (0.01, 0.02]`. Margin at 20 %: **3.20×**.

**The low-noise rows are not operational QUEST victories** and must not be described as
such. QUEST is credited with its noisy SLD QFI and has no compiled attaining readout, so
`F_C,CovQ < F_Q,QUEST` proves nothing operational. The reverse inequality is the strong
statement: a real QUEST measurement cannot exceed its own QFI ceiling, so
`F_C,CovQ > F_Q,QUEST` at 2 %, 5 %, 10 % and 20 % is a genuine CovQ win.

**Two independent channels, because seed spread is not enough.** Pilot-seed variance
measures the sampler's own noise and nothing else — it cannot detect a biased expectation,
because every seed would share the bias. The deployable template is therefore computed a
second way, deterministically: the exact channel enumerates the full binomial pilot support

  `E[g] = Σ_{k₁,k₂} Bin(k₁; n, p₀) Bin(k₂; n, p₁) · g(φ − atan2(−ĉ(k₁), d̂(k₂)))`

over all `201 × 201` outcomes at `n = 200` pilot shots per phase, and samples nothing. This
is legitimate for exactly the reason the scalar template is: an analyzer angle can only
scale a cat block's CFI, never rotate it, so the whole object is one scalar times the
matched-quadrature matrix.

| `q_edge` | Monte Carlo | exact enumeration | relative disagreement |
|---|---|---|---|
| 0.00 | 1.785484 | 1.785471 | `7.8e-6` |
| 0.01 | 1.821954 | 1.821938 | `8.8e-6` |
| 0.02 | 1.859557 | 1.859537 | `1.1e-5` |
| 0.05 | 1.979638 | 1.979607 | `1.6e-5` |

The exact channel's only free parameter is its angle grid, and it is converged: `4096`
versus `16384` points shift `κ` by at most `1.5e-6`, an order below the channel
disagreement. Sampling noise therefore dominates the residual, which is what one wants —
the deterministic channel is the reported value and the sampled one is the cross-check.

**Total uncertainty on the deployable exposure is `1.6e-5` relative.** The narrowest
crossing margin, `3.1 %` at `q_edge = 0.02`, exceeds it by roughly `1900×`. The QUEST arm
is a density-matrix simulation with no sampling at all, so it contributes none. The
crossover is a measured feature of the registered instance, not a numerical artefact, and
"CovQ requires less exposure at every tested point from 2 % onward" is a claim the numbers
support.

## bE coverage — a real limitation on the upper ratio

`bE` agrees with `tE` **exactly** where both ran: 40 emitted CX in total across the 8
instances covered, identical instance by instance. But best-position insertion costs
`(pool size) × (current depth)` per iteration and does not finish on deep targets, so
coverage is **8 of 11**, and the three uncovered instances are

  `path_m4` (17.78), `path_m5` (25.56), `banded_m5` (32.65)

— that is, **exactly the three largest ratios**. This is not a coincidence: they are deep,
which is both why `bE` is unaffordable there and where `bE` would have the most room to
improve on `tE`. So the *upper* end of the `2.50 – 32.65` range rests on `tE` alone, and a
referee is entitled to treat `32.65` as an upper bound on the separation rather than a
measured one. The shallow half of the table is cross-validated; the headline maximum is
not.

## Manuscript sentences that must change

| location | current | required |
|---|---|---|
| abstract | "faithful QUEST reimplementation" | QUEST-tE as published, with the joint reoptimisation phase |
| abstract | "factors of $2.5$–$135.6$" | `2.5`–`32.7` |
| abstract | "ten converged instances… one harder path instance reached the QUEST depth cap" | eleven converged; no cap reached |
| abstract | "at $20\%$… 2.759 versus $3.51\times10^5$" | `2.799` versus `8.961`, a `3.20×` ratio |
| abstract / §results | "less exposure at every tested point from $1\%$… onward" | from `2 %` onward; the optimistic QUEST ceiling is lower at `0` and `1 %` |
| §results-n10 | "44 two-qubit rotations in the QUEST preparation" | 3 two-qubit Pauli rotations, **6 emitted CX** |
| §results-n10 | "QUEST/CovQ exposure ratio is approximately $1.27\times10^5$" | `3.20` |
| §results-quest-exact | "On the ten converged targets the ratio ranges from $2.5$ to $135.6$" | eleven targets, `2.5` to `32.7` |
| §limitations | "faithful reimplementation, not official QUEST source" | keep, but say which variant (tE) and that bE agrees |

## What this does not touch

No CovQ theorem, certificate, compiler result, readout result, estimator result, or noise
theorem depends on the QUEST baseline. The pair-width characterisation, the cost identity,
the conic dual, the block-local pricing theorem, attainability, the estimator, the pinching
proposition and the fixed-setting negative result are all unaffected.
