# CovQ — registration

`SCIENTIFIC_SCOPE = FROZEN` as of 2026-08-21.
`NEW_EXPERIMENTS = PROHIBITED_UNLESS_RELEASE_RERUN_FAILS`.

Submission identity: **theory and certified prototype**. This is not a hardware
demonstration and must not be positioned as one.

## Registered and immutable

No target, seed, tolerance, grid, restart rule, cost model, or baseline rule below may
change. A release rerun that disagrees beyond the registered tolerance **pauses
submission**; it does not reopen the scientific scope.

| item | value |
|---|---|
| pair-width verification | 160 instances, `m ∈ {3,4,5,6}`, seeds from `PROTOCOL_SEEDS` |
| bipartite corollary | 225 instances |
| width-three boundary | 6 graphs |
| frozen controls | `F(+)_ij = +0.4`, `F(−)_ij = −0.4`, `m = 3` |
| K3 exact-target families | matching / path / toeplitz at `m ∈ {3,4,5}`, star `m=4`, banded `m=5` |
| K3 QUEST depth cap | 200 rotations, 4 random restarts after `|+⟩^m` |
| noise sweep | edge depolarization `{0, .01, .02, .05, .08, .10, .12, .15, .20, .30}` |
| base noise | per-qubit dephasing `0.02`, idle dephasing `0.01` |
| floor target | `G_req = 0.6 J + 0.6 I`, `m = 4` |
| deployable pilot | `f = 0.02`, phases `{0, π/2}`, retained in likelihood, globally fixed |
| estimator MC | 4000 replicates, shots `{500, 5000, 50000}` |
| SLD regression | `‖F_spectral − F_SLD‖_max ≤ 1e-10` on full-rank states; PSD floor `−1e-10` |
| LP feasibility tolerance | `1e-7` for the deployable re-price (HiGHS floor), `1e-8` elsewhere |

## Not claimed

- real hardware demonstration
- universal noise advantage
- exact fixed-setting optimization at scale
- an optimized sparse-state baseline
- coherent-flag resource optimality
- routed or crosstalk-aware compilation
- an attaining noisy readout for QUEST

## Preserved negative results

These must survive into the manuscript. Removing them would misrepresent the work.

| result | status |
|---|---|
| fixed-setting greedy support selection | `FAILS` — 7.24 % worst excess |
| QUEST at zero two-qubit noise | `QUEST_BETTER` — 1.626 vs 1.766 shots |
| `shot_scaled_edge_compile` | `LIMITED` — stalls once entanglement activates |
| `K1b` end-to-end compiler collision | `NOT_ESTABLISHED` |

## Release gates

`J0` scope freeze · `J1` independent SLD regression · `J2` frozen-pilot N10 ·
`J3` clean pipeline · `J4` reproduction diff · `J5` documentation ·
`J6` machine-generated paper · `J7` proof and citation audit · `J8` archive.
