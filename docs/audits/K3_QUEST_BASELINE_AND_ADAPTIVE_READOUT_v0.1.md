> ## ⚠ RETRACTED IN PART — superseded by `K3_QUEST_BASELINE_v0.2.md`
>
> Every QUEST number in this document was produced by a baseline that is **not
> QUEST**. Published QUEST runs two phases per iteration — insert a rotation,
> then *jointly reoptimise all accumulated angles*. The implementation used here
> omitted the second phase entirely, which understated the baseline by roughly
> an order of magnitude. The corrected results are in
> `K3_QUEST_BASELINE_v0.2.md`; the fault is diagnosed in
> `QUEST_FIDELITY_AUDIT_v0.4.md`.
>
> The adaptive-recentering half of this document (C10e / M5) is unaffected and
> stands.

# K3 — the QUEST baseline, and where the operating point comes from

Evidence: `prototype/src/covq/quest.py`, `prototype/src/covq/estimator.py`,
`results/baselines_k3_quest.json`, `results/measurements/adaptive_recentering.json`.

## What was compared, and what was deliberately not

QUEST is arXiv:2605.02367, *Quantum State Engineering Under Multiple Expectation-Value
Constraints* — "Quantum Unitary Engineering of States to Target". It builds the state as a
**depth-adaptive sequence of Pauli rotations**, each chosen to descend a
**sum-of-squared-residuals** cost, with no high-dimensional optimizer. Reimplemented here
from the published method description, not ported from released source; read every number
below with that caveat attached.

Per the v0.3 handoff, QUEST is a baseline for **exact first/second-moment targeting** and
is *not* an information-floor optimizer. It receives one target point; the floor compiler
receives a Loewner constraint `AᵀFA ⪰ G_req` and searches over every `F` meeting it.
Forcing QUEST into that role would need a separately specified outer optimization over
targets, which was not attempted. **K3 is run on exact targets only.**

Both sides receive the identical constraint set — `⟨P_i⟩ = 0`, `⟨P_iP_j⟩ = F_ij` — and
both sides' two-qubit counts are parsed from emitted circuits (gate C9), never from a
formula.

## Two fairness corrections that changed the numbers

Both are recorded because both moved the result in QUEST's favour.

1. **`|0…0⟩` is a trap for this constraint class.** It is a joint `Z` eigenstate, so every
   `⟨Z_i⟩` is maximally wrong and every `Z`-type pool element is a no-op on it. Greedy
   descent stalled on `path`, `toeplitz` and `star` targets it clears easily from `|+⟩^m`
   — the origin of the moment space, where all first and second moments vanish. A first
   table generated from `|0…0⟩` showed QUEST failing outright at residual `10⁻¹`. That
   table was wrong and was discarded.
2. **Path targets are restart-sensitive.** `path_m4` grinds past 400 rotations from
   `|+⟩^m` but converges in 83 from a random start. Reported numbers are best-of-restarts.

Before drawing any conclusion from a stall, the target was checked for single-pure-state
realisability: the Carathéodory-sparse single state hits every one of these targets to
`10⁻¹⁶`. So a stall is always an optimization outcome here, never evidence that the
target is unreachable by one pure state.

## Result

| instance | CovQ settings | CovQ E[cx]/shot | CovQ 2q depth | QUEST rot | QUEST cx | QUEST 2q depth | QFIM err | stop | cx ratio |
|---|---|---|---|---|---|---|---|---|---|
| `matching_m3` | 2 | 0.800 | 1 | 1 | 2 | 2 | 1.2e-11 | converged | 2.5× |
| `path_m3` | 3 | 0.900 | 1 | 61 | 122 | 122 | 2.3e-07 | converged | 135.6× |
| `toeplitz_m3` | 4 | 0.823 | 1 | 2 | 4 | 4 | 1.6e-09 | converged | 4.9× |
| `matching_m4` | 2 | 1.600 | 1 | 2 | 4 | 2 | 7.4e-09 | converged | 2.5× |
| `path_m4` | 3 | 1.350 | 1 | 86 | 158 | 142 | 1.5e-07 | converged | 117.0× |
| `toeplitz_m4` | 5 | 1.338 | 1 | 3 | 6 | 4 | 3.7e-09 | converged | 4.5× |
| `matching_m5` | 2 | 1.600 | 1 | 2 | 4 | 2 | 7.4e-09 | converged | 2.5× |
| `path_m5` | 3 | 1.800 | 1 | 200 | 356 | 296 | 1.9e-04 | max_depth | 197.8× |
| `toeplitz_m5` | 9 | 1.868 | 1 | 4 | 8 | 6 | 4.9e-09 | converged | 4.3× |
| `star_m4` | 4 | 0.900 | 1 | 56 | 112 | 94 | 2.2e-07 | converged | 124.4× |
| `banded_m5` | 5 | 1.470 | 1 | 59 | 118 | 80 | 2.0e-07 | converged | 80.3× |

10 of 11 converged. `path_m5` reached residual `1.9e-04` at the 200-rotation cap after
five starts; it is converging, slowly, not failing.

**Free cross-check.** CovQ's `E[cx]/shot` column is measured from the emitted schedule,
and it equals `Σ_e |F_e|` to machine precision in all 11 instances (max deviation
`4.4e-16`). That is the pair-activation cost identity of Theorem 2 / gate C6, confirmed
here independently of the test that was written for it.

## Reading it honestly

**Two regimes, and only one of them is a real separation.**

- `matching` and `toeplitz` targets are *easy* for QUEST: 1–4 rotations, 2–8 CX. CovQ's
  expected 2q cost is lower (2.5×–4.9×) but these are the same order of magnitude. On a
  device where switching settings costs anything material, a 2.5× edge in expected gates
  against a 2-to-9× disadvantage in settings is **not** a demonstrated win. The v0.3 cost
  declaration `C` includes a setting cost precisely for this reason, and it is unpriced
  here.
- `path`, `star` and `banded` targets separate by 80×–198×. That is a real gap and it is
  structural, not a tuning artefact: these targets need many edges at moderate weight, and
  QUEST must realise all of them coherently in one state while CovQ spreads them across
  branches and pays only `Σ_e |F_e|` in expectation.

**The one unconditional CovQ advantage is depth.** CovQ's two-qubit depth is **1 in every
instance**, and that is a theorem rather than a measurement: a branch is a matching, so
its Bell pairs are disjoint by construction and execute in parallel. QUEST ranges 2–296.
For a noisy device this matters more than the gate count.

**What K3 does not license.** No claim of advantage over QUEST on information floors —
QUEST was not run on that problem and the handoff forbids forcing it there. No claim of
end-to-end hardware advantage: settings, calibration, and setting-switch costs are
unpriced on both sides.

## Adaptive recentering: closing the M2 objection

M2 established that the matched readout must be phased to the operating point. The obvious
objection is that the operating point is what we are trying to estimate. A two-stage
protocol resolves it, and the pilot turns out to be nearly free.

Stage 1 spends a pilot fraction split between analyzers `A = 0` and `A = π/2`, which
together determine each block phase **including its sign** — one setting cannot, since
`cos` is even. Stage 2 re-phases every block to `A_B = φ̂_B − π/2`. The final estimate is a
single MLE over all counts from both stages, so no pilot data is discarded.

`m = 3`, `θ = (0.6, −0.4, 0.9)`, 800 replicates, against the **oracle** bound `F⁺/N`:

| pilot fraction | N | max abs bias | worst efficiency | naive `1/(1−f)` | median margin |
|---|---|---|---|---|---|
| 0.05 | 20 000 | 1.6e-04 | 1.083 | 1.053 | 0.9947 |
| 0.10 | 20 000 | 9.9e-04 | 1.061 | 1.111 | 0.9975 |
| 0.20 | 20 000 | 5.6e-04 | 1.112 | 1.250 | 0.9987 |
| 0.40 | 20 000 | 6.1e-04 | 1.061 | 1.667 | 0.9994 |

Efficiency stays near 1.06–1.11 **independently of the pilot fraction**, well below the
discard-the-pilot penalty. Two reasons: the pilot counts enter the final likelihood, and
the pilot settings are themselves informative — at least one of `|sin φ|`, `|cos φ|` is
always `≥ 1/√2`, so no pilot shot is wasted the way a genuinely singular setting would be.

**A finding that is not incidental.** The final MLE seeded at an arbitrary point lands in
the wrong mode often enough to inflate the variance ratio to 10–2400×. Seeded from the
pilot — least squares on the design matrix against the measured block phases — it is
efficient. The periodic likelihood is genuinely multimodal, so the pilot stage is not
merely a convenience for phasing the analyzer; it is what makes the estimator identifiable
at all. That is a second, independent reason `Readout` cannot be optional in the output
package.

## Status changes

`QUEST_BASELINE: ABSENT → RUN_EXACT_MODE_ONLY`.
`ADAPTIVE_RECENTERING_POLICY: ABSENT → DONE`.
`K3` is answered for exact targets and remains open for information floors.
