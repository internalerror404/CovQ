# CovQ Manuscript-to-Repository Integration Contract

## Purpose

The coding repository should implement the scientific object defined by manuscript v0.2, not the superseded v0.1 prototype. This note gives the minimum handoff contract between the paper and the codebase.

## Recommended repository layout

```text
CovQ/
├── README.md
├── STATUS.md
├── REGISTRATION.md
├── pyproject.toml
├── configs/
│   ├── experiment_v0.2.yaml
│   ├── tolerance_policy_v0.2.yaml
│   └── benchmark_manifest_v0.2.json
├── src/covq/
│   ├── semantics/
│   ├── matching/
│   ├── full_sector/
│   ├── floors/
│   ├── circuits/
│   ├── verification/
│   ├── baselines/
│   └── reporting/
├── tests/
├── results/
│   ├── canonical/
│   ├── failures/
│   └── generated/
├── paper/
│   ├── main.tex
│   ├── references.bib
│   ├── sections/
│   ├── figures/
│   └── generated/
├── docs/
│   ├── audits/
│   └── decisions/
└── archive/v0.1/
```

## Required semantic modes

The public IR should distinguish at least:

```text
semantics:
  labeled_schedule
  coherent_flag
  dephased_flag
  unlabeled_mixture_control

compile_mode:
  exact_pair
  pair_projection
  fixed_shot_information_floor
  shot_scaled_information_floor
  full_sign_polytope_control

frame:
  physical_Z
  logical_Z_with_clifford_wrapper
```

A coherent flag must never be reported as an interchangeable implementation of labeled branch width. Its flag dimension, Schmidt rank, controlled preparation, global entanglement, and joint-readout resources belong to a separate ledger.

## Core mathematical modules

### `matching/`

Must implement:

- graph support validation;
- signed edge extraction `x = offdiag(F)|_E`;
- degree inequalities;
- odd-set/blossom separation for nonbipartite graphs;
- local-only certification for bipartite graphs;
- matching decomposition;
- signed Bell schedule construction;
- independent certificate verification.

For every accepted exact pair target, verify:

```text
F_reconstructed = I + A_H(x)
abs(x) in MATCH(H)
activation_probability[e] = abs(x[e])
expected_pair_cost = sum_e c[e] * abs(x[e])
```

### `semantics/`

Must independently calculate:

```text
F_CQ       = sum_r p_r F_r
F_coherent = sum_r p_r F_r + Cov_p(mu_r)
F_dephased = sum_r p_r F_r
```

Tests must include nonzero and unequal branch means; otherwise the coherent correction is not exercised.

### `floors/`

Implement the two distinct convex programs:

1. fixed-shot normalized QFIM `F = I + A_H(x) >= G_star`;
2. shot-scaled total Fisher `F_tot = t I + A_H(y) >= G_star` with the homogenized matching cone.

Never repair an infeasible fixed-shot requirement by silently rescaling it. Return a typed infeasibility result or switch modes only through an explicit user/config decision.

### `full_sector/`

This is a control and baseline module, not the headline compiler. It should provide:

- exact small-instance sign-polytope membership/decomposition;
- global infeasibility certificates where available;
- sparse one-setting pure-state realization;
- labeled global-cat realization;
- coherent sparse data-state realization;
- explicit flag-data purification.

## Mandatory analytic controls

Use both three-generator targets from the manuscript:

```text
F_plus:
  diagonal: [1, 1, 1]
  off_diagonal: [0.4, 0.4, 0.4]
  expected: globally feasible, pair-width infeasible

F_minus:
  diagonal: [1, 1, 1]
  off_diagonal: [-0.4, -0.4, -0.4]
  expected: positive definite, globally infeasible
```

The first forces the code to distinguish resource rejection from physical rejection. The second forces the code to distinguish positive semidefiniteness from sign-correlation feasibility. Neither may be silently projected in exact mode.

## Status vocabulary

Use a machine-checkable status field:

```text
DONE
PENDING
FAILS
HOLLOW
ABSENT
SUPERSEDED
BACKGROUND
PRIOR_ART
```

Suggested `STATUS.md` headline entries:

```text
FULL_UNIT_DIAGONAL_FEASIBILITY       BACKGROUND
LABELED_CQ_ADDITIVITY                DONE
COHERENT_FLAG_FORMULA                DONE
PAIR_WIDTH_MATCHING_THEOREM          DONE
PAIR_ACTIVATION_COST_IDENTITY        DONE
WIDTH3_OPTIMIZATION_BOUNDARY         DONE
TASK_0C_V0_1                         SUPERSEDED
EXACT_COMPILER_CAMPAIGN_V0_2         PENDING
BASELINE_CAMPAIGN_V0_2               PENDING
NOISE_CAMPAIGN_V0_2                  PENDING
APPLICATION_CAMPAIGN_V0_2            PENDING
```

Here `DONE` for a theorem means a proof and theorem-level test fixture exist; it does not mean performance experiments have run.

## Freeze-before-run requirements

Before the first scientific campaign:

1. freeze benchmark generation, seed derivation, tolerances, and target IDs;
2. freeze exact correctness gates C1--C16 from the manuscript;
3. freeze baseline interfaces and matched resource budgets;
4. hash source, configuration, dependencies, and generated target packs;
5. record which hypotheses are directional and which are descriptive;
6. prohibit runtime tolerance relaxation and target regeneration after a failure.

Pilot runs may determine directional thresholds, but pilots must use disjoint target IDs and seeds. Once thresholds are frozen, the registered pack is immutable.

## Canonical record requirements

Every result record should include:

```text
schema_version
run_id
target_id
target_hash
config_hash
source_commit
seed
semantic_mode
compile_mode
generator_frame
hardware_graph_hash
solver_status
certificate_type
certificate_payload
qfim_target
qfim_realized
residuals
branch_schedule
emitted_circuit_hashes
logical_resources
physical_resources
measurement_contract
correctness_dependencies
status
failure_reason
```

Circuit-derived resources must be regenerated from emitted artifacts. Analytical estimates may be stored separately but cannot populate physical gate, depth, routing, or ancilla fields.

## Paper integration

Only canonical `DONE` records with recognized manifests and passed dependencies may generate:

```text
paper/generated/exact_results.tex
paper/generated/baseline_results.tex
paper/generated/noise_results.tex
paper/generated/scaling_results.tex
paper/generated/application_results.tex
paper/generated/result_manifest.json
```

The current manuscript intentionally contains pending result templates. Replace them through deterministic generation scripts; do not edit measured cells by hand.

## Minimum repository milestone before changing the abstract

The abstract's empirical sentence may be revised only after all of the following exist:

- matching and coherent-flag correctness gates pass independently;
- both infeasibility controls pass;
- sparse pure-state and QUEST interfaces are operational;
- emitted circuit resource regeneration passes;
- the exact campaign manifest is sealed;
- a canonical result table is reproducible from a clean clone.

Until then, retain the present statement that empirical compiler and hardware advantages are unclaimed.
