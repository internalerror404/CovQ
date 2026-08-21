# CovQ v0.12 review-response artifact

This directory records the targeted major revision prompted by two external reviews.

## What is new

1. **Universal twofold ceiling.** Pair width obeys `0 <= F <= 2I`, hence ideal pair-width exposure can never beat product exposure by more than a factor of two for any downstream map or floor.
2. **Mode-correct finite campaigns.** Proposition 12.3 is scoped to normalized fixed exposure; the new shot-scaled corollary samples `T` branches with `p_b=n_b/T` and shows absolute allocation margin grows as `L sqrt(T log(d/eta))` while relative margin shrinks as `T^{-1/2}`.
3. **Explicit cost model.** The flagship table is identified as `c0=1, ce=0` exposure accounting. Separate curves reoptimize the complete signed-matching pool over `ce/c0`.
4. **Contract breadth.** A deliberately multidirectional floor reduces the noiseless product/CovQ ratio from `1.823` to `1.199` and selects the pairless branch by `q=0.20`.
5. **Model sensitivity.** Additional points bracket the registered crossover in `(0.0125,0.015]` for per-CX depolarization and `(0.015,0.0175]` for per-CX amplitude damping. These are model-dependent brackets, not hardware thresholds.
6. **Specification-matched comparator.** Using the same branches, labels, and probabilities as CovQ, a Pauli-rotation construction emits exactly two CX per active pair versus one for CovQ. The separate single-state QUEST result remains a different specification.
7. **Numerical validation wording.** Deterministic binomial enumeration versus Monte Carlo is described only as numerical implementation validation, not physical-model uncertainty.

## What was not run

No remaining bE run was performed. The fully cross-validated one-state QUEST range remains `2.50-8.89x`; the three deeper tE-only values remain explicit upper bounds on separation from bE.

## Production correction

The v0.12 PDF removes the duplicated Section 8, restores Section 5.8, fixes equation and table numbering, repairs the contents and severed text, shortens the illustrative application, adds the public repository URL, and states that the archival DOI is pending rather than inventing one.

`revision_summary.json` is the compact machine-readable record. Full CSVs, figures, revision sources, and PDF validation reports are distributed with the v0.12 delivery bundle.
