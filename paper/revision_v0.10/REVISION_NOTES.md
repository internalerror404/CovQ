# CovQ manuscript protection and gain pass — v0.10

This revision implements the high-leverage protection pass and two research additions without expanding the paper into a real-device, routed, correlated-noise, or width-three heuristic study.

## Protection pass

- repairs the empty Section 1.6, Section 14.9, and Appendix F headings;
- repairs severed sentences in boxed text and corrects the pilot-sweep table caption;
- replaces the Marchenko–Pastur comparison with a 100,000-draw Gaussian parametric bootstrap, including a simultaneous 95% envelope across the two registered operating points;
- adds explicit positioning against distributed quantum sensing (Proctor–Knott–Dunningham, Eldredge et al., Rubio–Dunningham, and Shettell–Markham);
- leads the abstract and introduction with the resource phase diagram;
- states the pair-width spectral envelope `0 <= F <= 2 I` in the abstract.

## Headline hardening already incorporated

The deployable readout template is no longer defended only by variation across Monte Carlo pilot seeds. It is evaluated through two independent channels: exact enumeration of the full binomial pilot support and Monte Carlo pilot draws. Their largest relative disagreement across the registered noise grid is `6.0e-5`; the narrowest `3.09%` crossover margin is more than `512x` larger. This controls expectation bias that seed-to-seed spread alone could not reveal.

## Gain 1 — compiler-only scale

Adds a state-simulation-free campaign at `m = {50,100,200,300,500}` on rectangular grids and a connected heavy-hex-family construction. The solver uses a restricted LP master, minimum-eigenvector Loewner cuts, and exact maximum-weight-matching pricing. It reports wall clock, column count, oracle calls, gap, and contract slack. All 30 runs terminated with optimal status.

## Gain 2 — finite randomized campaigns

Adds a matrix-Bernstein proposition. For branch pullbacks `X_b = A^T F_{C,b} A`, pair width gives `0 <= X_b <= 2 ||A||_op^2 I`. If the compiler purchases margin `delta`, the empirical average of `N` randomized branch executions satisfies the original Fisher floor except with probability at most

`d exp[-N delta^2 / (2 L^2 + (2/3) L delta)]`, with `L = 2 ||A||_op^2`.

A sufficient sample count is

`N >= [8 ||A||_op^4 / delta^2 + 4 ||A||_op^2 / (3 delta)] log(d/eta)`.

## Free wins

- promotes the exact pilot-information decomposition to a lemma;
- names and freezes the comparator as **CovQ-QUEST Reference Reimplementation v1** and ships a machine-readable finite-difference gradient pin;
- moves the full unit-diagonal feasibility background from the main paper to Appendix G.

## Deliberately not included

The refined `q_CX` grid and amplitude-damping replication were not run in this revision. The paper therefore retains the careful statement “at every tested point from 2% onward” and does not upgrade it to a channel-robust threshold claim.

## Source notes

`covq_replacement_pages.tex` contains the fully typeset replacement pages. `assemble_covq_v010.py` assembles them with the cleaned v0.8 PDF, updates page numbers, metadata, and bookmarks, and inserts the compiler-scale page. `patch_two_channel_evidence.py` replaces the N10 discussion with the exact-versus-Monte-Carlo validation. The final rendered PDF is distributed separately because this repository interface stores the text revision sources and machine-readable evidence.
