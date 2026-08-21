# Eq (111) — the noisy operational information-floor compiler

Evidence: `prototype/src/covq/noise.py`, gates N1–N6/N8 in `prototype/tests/test_covq.py`,
records under `results/noise/`. The six earlier measurement records are untouched.

Implemented in the reformulated shape, not as a penalised norm-to-target heuristic:

  `min_{Π,n,R} C_total(Π,n,R)   s.t.   Aᵀ F_C^𝒩(Π,n,R) A ⪰ G_req`,

with `F_C^𝒩 = Σ_b n_b F_{C,b}^{𝒩,R_b}` the CFI of the **declared readout** under the
noise model. The mixed-state QFIM enters only as the reference upper bound
`F_C^𝒩 ⪯ F_Q^𝒩`; generator covariance is a diagnostic and nothing more — the pinching
proposition already showed it can be exactly constant while the program is dead.

## The result: block-local noise preserves matching-based compilation

**Theorem (verified, N4 + N5).** If the channel and the readout are block-local, every
branch template splits exactly as

  `F_{C,M,σ}^𝒩 = F_0^𝒩 + Σ_{e∈M} ΔF_{e,σ_e}^𝒩`,

and for a dual matrix `Q = AYAᵀ` the pricing subproblem collapses to

  `⟨Q,F_0^𝒩⟩ − c_0 + max_{M∈MATCH(H)} Σ_{e∈M} [ max_σ ⟨Q,ΔF_{e,σ}^𝒩⟩ − c_e ]`,

still a maximum-weight matching — the same oracle as ideal Eq (65), with reweighted
edges. **Noise changes the templates and the weights; it does not change the
combinatorics.**

| gate | check | result |
|---|---|---|
| N1 | noiseless limit reproduces the ideal signed matching | `4.4e-16` |
| N2 | dephased cat at quadrature gives `F_C = F_Q = v² ssᵀ` | exact, `v = (1−2p)^k` |
| N3 | complete pinching: covariance preserved, `F_Q` and `F_C` → 0 | `Cov` drift `0`, both `< 1e-12` |
| N4 | edge additivity of the noisy branch CFI | **exactly `0.0`**, m = 3,4,5 |
| N5 | pricing oracle vs exhaustive signed matchings | 12 cases, 0 disagreements, `4.4e-16` |
| N6 | primal/dual on the convex exposure problem | gap `< 1e-15` at every noise level |
| N8 | asymmetric readout needs a third analyzer point | two-point fit provably off quadrature |

Everything is computed on one- and two-qubit subsystems rather than the full `2^m` state.
That is the computational payoff of block-locality, and it is why the noisy compiler costs
no more than the ideal one.

### One structural refinement the ideal case does not need

**The pairless branch is a separate column.** It activates no pair, so it never waits
through a pair preparation and carries no idle dephasing — which puts it on a *different
reference template* from the one the edge weights are measured against. Every non-empty
pair-width branch has two-qubit depth exactly one and therefore imposes the same idle
wait, so the additive decomposition must be taken against the **idling** product template.

This is not a technicality. Using the idle-free reference inside pricing broke N5 in
**15 of 18 cases** before it was fixed, with errors up to `0.71`, and the pairless branch
won 19 of 24 priced instances once both candidates were considered. With
`idle_dephasing = 0` the two references coincide and this collapses to one matching
problem; with idle noise present, dropping the extra candidate loses the optimum whenever
idling costs more than the best edge gains.

### Declared scope

Holds for noise acting *within* blocks. **Not** claimed for crosstalk, correlated branch
noise, route collisions, or coherent errors coupling distinct blocks — those break the
factorisation the argument rests on, and are excluded from this campaign by construction.

## The compiler changes its answer, and there is a threshold

`m = 4`, `G_req = 0.6·J + 0.6·I`, `c_e = 0.1`, `c_0 = 1.0`, dephasing 0.02, idle 0.01:

| edge depolarization | primal cost | gap | settings | entangled |
|---|---|---|---|---|
| 0.00 | 2.11928 | 2.1e-16 | 3 | 3 |
| 0.05 | 2.34823 | 0 | 3 | 3 |
| 0.10 | 2.61639 | 1.7e-16 | 3 | 3 |
| 0.15 | 2.93325 | 1.5e-16 | 3 | 3 |
| **0.20** | **3.25521** | 1.4e-16 | **1** | **0** |
| 0.30 | 3.25521 | 1.4e-16 | 1 | 0 |

Past a threshold between `0.15` and `0.20` the pair primitive stops paying for itself and
the compiler abandons entanglement entirely; the cost then plateaus, because no edge is
used and edge noise no longer enters. A noise-aware objective that returned the ideal
schedule everywhere would be evidence it was doing no work. This one does not.

## Pricing the settings caveat

K3 left "CovQ needs 2–9 settings against QUEST's one" as an unpriced hedge. With batched
execution `N_switch = q − 1`, so both totals are affine in the shot count and

  `N* = (q_C − q_Q)·c_setup / (c_Q − c̄_C)`.

Crossover shots, cost unit = one two-qubit gate:

| instance | q_C | E[cx]_C | q_Q | cx_Q | c_setup=10 | 100 | 1000 | 10 000 |
|---|---|---|---|---|---|---|---|---|
| `matching_m3` | 2 | 0.800 | 1 | 2 | 8 | 83 | 833 | 8 333 |
| `path_m3` | 3 | 0.900 | 1 | 122 | 0 | 2 | 17 | 165 |
| `toeplitz_m3` | 4 | 0.823 | 1 | 4 | 9 | 94 | 944 | 9 441 |
| `matching_m4` | 2 | 1.600 | 1 | 4 | 4 | 42 | 417 | 4 167 |
| `path_m4` | 3 | 1.350 | 1 | 158 | 0 | 1 | 13 | 128 |
| `toeplitz_m4` | 5 | 1.338 | 1 | 6 | 9 | 86 | 858 | 8 580 |
| `matching_m5` | 2 | 1.600 | 1 | 4 | 4 | 42 | 417 | 4 167 |
| `toeplitz_m5` | 9 | 1.868 | 1 | 8 | 13 | 130 | 1 305 | 13 047 |
| `star_m4` | 4 | 0.900 | 1 | 112 | 0 | 3 | 27 | 270 |
| `banded_m5` | 5 | 1.470 | 1 | 118 | 0 | 3 | 34 | 343 |

Even at `c_setup = 10⁴` gate-equivalents — a deliberately pessimistic setup charge — the
worst crossover is **13 047 shots**, and the easy families cross at a few hundred. Real
metrology campaigns run `10⁴`–`10⁶` shots. So the settings penalty amortizes in every
tested case, and the K3 hedge should be replaced by this boundary rather than repeated.

**Optimization caveat, kept separate.** `q = |{b : n_b > 0}|` is a cardinality/fixed-charge
cost. Pricing it inside the objective makes the problem mixed-integer conic and voids the
oracle-polynomial claim. The convex per-use compiler above and this amortization boundary
are therefore reported as two distinct objects; no `λ_q q` term was added to the conic
program.

## What is still open

- **N7** deployable pilot-recentered arm inside the *noisy* compiler (the ideal-case
  version is closed as C10e).
- **N9** fixed-setting-cost solver against exhaustive support enumeration.
- **N10** a single campaign giving CovQ, QUEST, sparse-Carathéodory and product arms the
  same shots, noise map and accounting.
- QUEST under noise. Giving it an unpriced optimal POVM and calling the result operational
  would be exactly the asymmetry the handoff warns against; the fair conservative
  comparison hands QUEST its noisy QFI as an optimistic bound while CovQ reports emitted
  readout CFI. Not yet run.
- Crosstalk, coherent drift, Clifford-frame noise, coherent flags — all outside scope.
