# CovQ-QUEST Reference Reimplementation v1

This artifact names and freezes the QUEST comparator used by the CovQ manuscript.

## Scope

- **Method:** QUEST-tE (terminal exact insertion) with the published second phase: joint L-BFGS reoptimization of every accumulated angle after each insertion.
- **Cross-check:** QUEST-bE (best insertion position) where the registered compute budget permits it.
- **Input:** exact first- and second-moment targets only.
- **Not claimed:** QUEST is not treated as an information-floor optimizer.
- **Initialization:** \(|+\rangle^{\otimes m}\), included in the emitted self-contained circuit.
- **Accounting:** logical Pauli rotations, weight-two rotations, emitted CX count, and CX depth are separate fields; noise is applied once per emitted CX.

The active implementation is `prototype/src/covq/quest.py`. The canonical exact-target record is `results/baselines_k3_quest.json`. The earlier terminal-greedy routine that omitted joint angle reoptimization remains only as a named ablation and must not be reported as QUEST.

## Finite-difference gradient pin

The joint-angle adjoint gradient is pinned by `test_quest_published_gradient_matches_finite_differences` in `prototype/tests/test_covq.py`:

- deterministic seed: `3`
- three-qubit fixture
- three mixed-weight Pauli rotations: `X0Y1`, `Z2`, and `X1X2`
- centered finite-difference step: `1e-6`
- componentwise absolute tolerance: `1e-6`

See `fd_gradient_pinning.json` for the machine-readable record.

## Provenance

- QUEST correction commit: `f141111cab92a2ac0f7c1cf64466453c6a7c6194`
- clean CovQ release commit used as the artifact anchor: `8cca6ef106b16a37a891077da88adbe8666d4c30`
- canonical result source commit recorded by the release artifact: `daa5587ee9e6a9d4791ea0855b8f38a68ccff336`

This is a reimplementation from the published method description, not official QUEST source.
