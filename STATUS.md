# CovQ — STATUS

**`SCIENTIFIC_SCOPE = FROZEN`** (2026-08-21). `NEW_EXPERIMENTS = PROHIBITED_UNLESS_RELEASE_RERUN_FAILS`. Registration: `REGISTRATION.md`. Submission identity: theory and certified prototype — **not** a hardware demonstration.

Last updated: 2026-08-21. Decision of record: `docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md`.
Current documents: `charter/*_v0.3.*`. Proof of record:
`docs/proofs/THEOREM_PAIR_WIDTH_AND_INFORMATION_FLOOR_v0.1.md`.

**The spine.** information requirement → exact pair-width geometry → certified minimum-cost
program → hardware-native Bell-pair schedule. The `k = 2` theorem is the engine; the
information-floor compiler is the reason anyone needs the engine.

## Journal release gates (v0.4)

| gate | status | note |
|---|---|---|
| `J0` scope freeze | `DONE` | `REGISTRATION.md`; targets, seeds, tolerances, grids, restart rules, cost models and baseline rules all pinned. |
| `J1` independent SLD regression | `DONE` | `prototype/tests/test_mixed_qfim_sld_regression.py`. Lyapunov solve vs spectral Eq (110): `1.33e-15` on 30 full-rank complex states, against the registered `1e-10`. Rank-deficient states are singular by construction and are checked as a linearly convergent limit instead. |
| `J2` frozen-pilot N10 | `DONE` | `f = 0.02`, phases `{0, π/2}`, pilot retained, globally fixed. Deployable overhead `1.10 %` → `1.44 %`; never below the oracle ceiling. |
| `J3` clean pipeline | `PASS` | All 14 canonical records regenerated from a clean tree: `dirty=false`, `source_commit=dc22135`, 0 non-conforming. Took three passes — the first two were dirtied by my own concurrent edits, which is what the flag is for. |
| `J4` reproduction diff | `PASS` | 27 headline quantities compared, **0 changed, 0 absent**. Also a determinism check: two independent reruns of the same tree agree exactly. `docs/audits/JOURNAL_RELEASE_REPRODUCTION_v0.4.md`. |
| `J5` documentation | `DONE` | `prototype/README.md` rewritten (it described an 85-test tree with no noise, estimator or QUEST work), `REGISTRATION.md` added, `STATUS.md` refreshed. |
| `J6` machine-generated paper | `PASS / FIGURES_ABSENT` | LaTeX source vendored at `paper/CovQ_Paper_v0.4.tex`. All five result tables regenerate from records into `paper/generated/`; four **MATCH** the source to 0.5 % relative and N10 is `EXPECTED_CHANGE` (it gains the deployable column). Seven prose-quoted figures also bound and matching. 27/27 quantities resolved, all records clean, single source commit. **`make_figures.py` and the `.bib` never arrived**, so Figures 3–4 cannot be regenerated or checksummed. |
| `J7` proof and citation audit | `PARTIAL` | `docs/audits/J7_PROOF_AUDIT_v0.4.md`. Self-audit only — it does not discharge the requirement for an independent reader. Four items flagged. Citations **not verifiable**: `arxiv.org` is blocked by the egress proxy. |
| `J8` archive | `PARTIAL` | `results/RELEASE_MANIFEST.json`: SHA-256 over every source and result file plus an environment lock (70 files, `6f2b2a45896c6826…`). Tag `covq-v0.4-journal-release` created locally at `455965e8ef23` but **the git proxy in this environment refuses tag refs** — branch pushes succeed, `refs/tags/*` does not. Recreate it after clone with `git tag -a covq-v0.4-journal-release 455965e8ef23`. **DOI minting is external**; figure and PDF checksums need those artifacts in-repo. |

## Disposition

| Item | Status |
|---|---|
| `THEOREM_A_FULL_FEASIBILITY` | `BACKGROUND` |
| `SIGNED_CAT_PRIMITIVE` | `BACKGROUND` |
| `LABELED_FLAG_ADDITIVITY` | `CORRECT / SEMANTICS` |
| `COHERENT_FLAG_FORMULA` | `CORRECTED` |
| `THEOREM_C_K2_MATCHING` | `PAPER_GRADE_CANDIDATE` |
| `THEOREM_D_NORM_STABILITY` | `BACKGROUND` |
| `THEOREM_E_FULL_WIDTH_HARDNESS` | `PRIOR_ART` |
| `TASK_0C_V0_1` | `SUPERSEDED` |
| `TASK_0B_5` | `NEXT` |

## Additions to the disposition

| Item | Status | Note |
|---|---|---|
| `THEOREM_C_K2_MATCHING_HARDWARE` | `PAPER_GRADE_CANDIDATE` | The `MATCH(H)` form, stronger than the `K_m` form. Implemented and exercised. |
| `WIDTH_2_VS_3_COMPLEXITY_BOUNDARY` | `OPEN / LIKELY_COROLLARY` | Partition into Triangles is classically NP-complete, so the reduction should go through in a paragraph. Re-scoped: the risk is no longer that it is false but that it is folklore in the bounded-cluster partitioning literature. |
| `TASK_0B_5_STEP_1_PRIOR_ART` | `DONE` | `docs/audits/PRIOR_ART_SWEEP_TASK0B5_v0.1.md`. K2 does not fire, conditional on framing. |
| `K2_MATCHING_TRANSLATION` | `NOT_FIRED / CONDITIONAL` | Survives only as identification-plus-achievability, never as new polyhedral combinatorics. |
| `SUPPORT_FUNCTION_RECOVERS_K_PRODUCIBILITY` | `VERIFIED` | `max{1ᵀF1 : F ∈ Q^prog_{m,k}} = ⌊m/k⌋k² + r²`, exact in all 15 cases with `m ≤ 7, k ≤ 3`; closed form `m + 2⌊m/2⌋` at `k = 2` up to `m = 12`. |
| `COHERENT_FLAG_WIDTH_LOOPHOLE` | `CONFIRMED_NUMERICALLY` | Conditional data width is not a resource for coherent flags. See below. |
| `INFORMATION_FLOOR_COMPILER` | `PROMOTED_TO_MAIN_PROBLEM` | `min Cost(Π) s.t. AᵀF_ΠA ⪰ G_req`. Theorem 2, closed. Convex; LP over `MATCH(H)` + Loewner cutting planes; certificate-producing. |
| `THEOREM_1_PAIR_WIDTH` | `PROVED (2A–2D)` | Necessity, magnitude sufficiency, hardware form, constructive emission. Every branch has two-qubit depth exactly 1. |
| `THEOREM_2_CONTRACT` | `PROVED (2E)` | Objective **is** the expected Bell-pair bill via `t_e = |F_e|` and down-closedness. |
| `THEOREM_3_CERTIFICATES` | `PROVED (2F)` | Feasible / lower bound / infeasible, each independently checkable. |
| `ODD_EVEN_MODE_CEILING` | `PROVED + VERIFIED` | Width-2 caps `u_SᵀFu_S` at 2 (even `k`) and `2 − 1/k` (odd `k`). Blossom as physics. Exact for `k = 2..8`. |
| `COMMON_MODE_CURVE` | `PROVED + VERIFIED` | `cost*(γ) = m(γ−1)/2`; compiled cost matches to 1e-8 with zero optimality gap. |
| `SIGN_PARITY_AT_WIDTH_3` | `NOTED` | Inside a 3-block `(s_is_j)(s_js_k)(s_is_k) = +1`, so magnitudes cannot decide feasibility at `k ≥ 3`. Explains why `k = 2` is clean, and where the width-3 difficulty enters. |
| `K9_FLOOR_ADDS_NOTHING` | `NEW / LIVE` | If every benchmark contract is met by the product probe or one fixed state, the composite claim is empty. |
| `SUPPLEMENTARY_SWEEP_FLOOR_INTERFACE` | `DONE` | Fisher-information LMI constraints in minimum-cost sensor design are established. The interface alone is **not** claimable. |
| `SUBUNIT_DIAGONAL_SECTOR` | `REJECTED_FOR_V1` | Needs joint `(μ, M)` moment feasibility, not `M ∈ Q_m`. Not audited. |
| `K1_SPLIT` | `DONE` | K1a specification collision: **FIRED**. K1b end-to-end compiler collision: **NOT ESTABLISHED**. |
| `QUEST_BASELINE_V0_1` | `RETRACTED` | Not QUEST: it omitted the joint angle-reoptimisation phase. All v0.1 QUEST resource numbers withdrawn. See `docs/audits/QUEST_FIDELITY_AUDIT_v0.4.md`. |
| `QUEST_BASELINE` | `RUN_EXACT_MODE_ONLY` | arXiv:2605.02367, reimplemented from the published method description. `docs/audits/K3_QUEST_BASELINE_AND_ADAPTIVE_READOUT_v0.1.md`, `results/baselines_k3_quest.json`. **QUEST-tE as published** (insert, then joint L-BFGS reoptimisation). `docs/audits/K3_QUEST_BASELINE_v0.2.md`. **11/11 converge**, ratio range **2.50–32.65** (was 2.5–197.8 against the flawed baseline). `matching`/`toeplitz` unchanged at 2.5–4.9×; `path`/`star`/`banded` fall to 8.9–32.7×. CovQ 2q depth is 1 in every instance, by the matching theorem, and does not depend on the baseline. `bE` verified to agree with `tE`. Not run on information floors. |
| `SPARSE_PURE_STATE_BASELINE` | `IMPLEMENTED` | Carathéodory-sparse single state, plus three realisation-matched siblings. |
| `MANUSCRIPT_INTEGRATION_MD` | `RECEIVED` | Arrived 2026-08-21 with `CovQ_Repository_Handoff_v0.3.md` and `CovQ_Experiment_Protocol_v0.3.yaml`. Unblocks the record schema, status vocabulary and freeze-before-run rules. The gate table is C1–C14 (v0.3), not C1–C16; `C2a/b/c` already match the prototype. |
| `GATE_ALIAS_MAP` | `DECLARED` | `M1→C10a`, `M1_schedule→C10b`, `M2→C10c`, `M3_qfi→C15a`, `M3_sep→C15b`, `M3_dephase→C15c`. Composite `C10` is split so a readout failure cannot invalidate unrelated matrix-pullback code. Commit `4a261a8` artifacts are **not** retroactively renamed. |
| `FROZEN_CONTROLS` | `FROZEN` | `F(+)_ij = 0.4` globally feasible / pair infeasible; `F(-)_ij = -0.4` PD / globally infeasible. Both verified with the three-way classifier. |
| `SHOT_SCALED_EDGE_MODE` | `LIMITED` | Agrees with the branch form to 1e-15 when both converge, but stalls once entanglement activates. Needs an interior-point SDP backend. |
| `MANUSCRIPT_V0_2` | `AUDITED` | `docs/audits/MANUSCRIPT_V0_2_NUMERICAL_AUDIT.md`. Thm 5.3, Cor 5.4/5.5, Thm 7.1, Eq (65), Thm 8.2, Eq (63) all verified numerically. |
| `SHOT_SCALED_FLOOR` | `IMPLEMENTED / SUPERSEDES_A_BASED` | Manuscript Problem 8.1 with the Eq-(65) matching oracle. Strong duality exact. The `A`-based version is retained as the downstream-projection mode. |
| `THEOREM_3_DUAL_EXPLICIT` | `DONE` | Conic dual derived and implemented for both formulations; Farkas rays certify infeasibility; separation oracle is an explicit interface. |
| `K9_ANSWERED` | `CLOSED_FORM` | Pair entanglement beats the product probe iff `c_e < 2c₀/m`; max advantage exactly `2×`. With `c_e = c₀` the floor compiler is trivial — must be stated in §8. |
| `WIDTH_3_HARDNESS` | `PROVED_IN_MANUSCRIPT` | Partition Into Triangles, §7.2; verified on 6 graphs. Supersedes the repo's `OPEN` status. |
| `SECTION_11_MEASUREMENT` | `AUDITED / PARTIAL` | `docs/audits/SECTION_11_MEASUREMENT_COMPILATION_AUDIT_v0.1.md`. Local attainability closed operationally for cat-block schedules; minimum-cost measurement synthesis remains the frontier §11.5 declares. |
| `LOCAL_READOUT_ATTAINS_QFIM` | `PROVED + VERIFIED` | A product of single-qubit equatorial measurements attains every branch QFIM exactly (`< 3.6e-15`), for block width `k = 2, 3, 4`. Label retention then gives `CFI_schedule = F`. No ancilla, no joint measurement. |
| `READOUT_IS_NOT_OPTIONAL` | `NEW / VERIFIED` | The fixed all-X readout returns **identically zero** information at `θ = 0`, and loses every mixed-sign block along the whole uniform ray. §11.1(ii) is load-bearing; `Readout` in the Eq (114) output tuple cannot be defaulted. |
| `EQ_110_IMPLEMENTED` | `VERIFIED` | Reduces to `4 Cov` on pure states to `4.4e-16`; matches the closed form `F_Q = v²` for the depolarised probe exactly. |
| `Z_CORRELATORS_BLIND_TO_DEPHASING` | `NEW / VERIFIED` | All generators are `Z`-diagonal, so `Z` dephasing leaves the whole covariance matrix bit-for-bit invariant while the SLD QFIM decays to zero. At `p = 0.5` a Bell block reports `4 Cov = 2` and carries no information at all. Strengthens the §11.3(a) warning into a concrete statement. |
| `MEASUREMENT_GATE_NUMBERING` | `DEFERRED` | New gates are `M1`–`M3`, not `C13`–`C16`. The contract that fixes the C-series numbering has not been received; renumbering is deferred rather than guessed. |
| `MIN_COST_MEASUREMENT_SYNTHESIS` | `OPEN` | M1 exhibits *an* attaining readout, not the cheapest under a declared setting cost. Eq (112)'s estimator layer and §11.4's noisy objective (111) are also unimplemented. |
| `DIRECT_Z_BLOCK_READOUT_EXISTENCE` | `DONE` | General-`k` proof + `k = 2,3,4` validation. |
| `DIRECT_Z_READOUT_2Q_GATE_COST` | `OPTIMAL_ZERO` | All but one qubit per block in `X`; one pivot carries the analyzer angle. |
| `DIRECT_Z_READOUT_ANCILLA_COST` | `OPTIMAL_ZERO` | No ancilla, no joint Bell measurement, no entangling readout. |
| `SCHEDULE_LEVEL_CFI_EQUALS_QFIM` | `DONE` | `F_C^schedule = F_Π`; worst gap `2.2e-15` over 48 cases. |
| `READOUT_SINGULARITY_CONTROL` | `DONE` | Failure set `φ_B − A_B ∈ πZ`; zero **first-order score**, model nonregular there — not distinguishability-free. |
| `MATCHED_QUADRATURE_READOUT` | `DONE` | Exact root `α* = atan2(−c, d)` from two parity evaluations; regularity margin `η_ro = 1` by construction. Replaced a grid argmax that returned an arbitrary detuning. |
| `MIXED_STATE_QFIM_IMPLEMENTATION` | `DONE` | Eq (110); reduces to `4 Cov` on pure states to `4.4e-16`; reproduces `F_Q = v²` exactly. |
| `COVARIANCE_AS_NOISY_QFI` | `REFUTED` | Pinching onto the generator eigenbasis preserves `Cov(P)` exactly (`4.9e-15`) while `F_Q ≡ 0` (`4.7e-28`). General proposition, verified on random commuting families. |
| `NOISY_READOUT_GAP_AT_P_0_2` | `RETRACTED` | Was a grid tie-break artefact of this prototype, not physics. `F_C(v,δ) = v²sin²δ/(1−v²cos²δ)` reproduces the erroneous value to 15 digits; quadrature attains `F_Q` at every visibility. |
| `GENERAL_MIN_COST_MEASUREMENT_SYNTHESIS` | `PARTIAL` | Optimal in 2q-gate and ancilla count. Unpriced: local basis changes, setting changes, calibration, robustness, Clifford unwrapping, coherent-flag joint readout. |
| `CLIFFORD_FRAME_READOUT_COST` | `ABSENT` | |
| `COHERENT_FLAG_JOINT_READOUT` | `ABSENT` | |
| `ESTIMATOR_LAYER_EQ_112` | `DONE` | Branch-conditioned MLE over block-parity counts. Sufficient statistic is one parity count per `(branch, block)`; its Fisher matrix equals the schedule QFIM exactly (`0.0`). Consistent, and saturates `F_Π⁺/N` on the identifiable quotient — efficiency eigenvalues in `[0.992, 1.039]` against a Marchenko–Pastur band of `[0.946, 1.056]` at `N = 5·10⁴`, `R = 4000`. Gate alias `C10d` / `M4`. `results/measurements/estimator_efficiency.json`. |
| `NOISE_AWARE_OBJECTIVE_EQ_111` | `ABSENT` | M3 supplies the measuring stick; no compiler optimises it. |
| `ADAPTIVE_RECENTERING_POLICY` | `DONE` | Two-stage protocol; pilot split between `A = 0` and `A = π/2` fixes each block phase including sign. Efficiency `1.06–1.11` against the **oracle** `F⁺/N` at every pilot fraction from 0.05 to 0.40, far below the `1/(1−f)` discard penalty, because the pilot enters the final likelihood. Gate alias `C10e` / `M5`. |
| `LIKELIHOOD_MULTIMODALITY` | `NEW / VERIFIED` | The periodic likelihood is genuinely multimodal: an arbitrarily seeded MLE inflates the variance ratio to 10–2400×. Seeded from the pilot it is efficient. The pilot stage is therefore what makes `θ` identifiable at all — a second, independent reason `Readout` cannot be optional. |
| `EQ_111_NOISE_AWARE_COMPILER` | `DONE / SCOPED` | `docs/audits/EQ_111_NOISE_AWARE_COMPILER_v0.1.md`. Implemented as the noisy operational floor `min C s.t. AᵀF_C^𝒩A ⪰ G_req`, not a norm-to-target penalty. Convex; primal/dual gap `< 1e-15` at every noise level. |
| `BLOCK_LOCAL_PRICING_THEOREM` | `PROVED + VERIFIED` | Block-local noise reweights edges but preserves the max-weight matching oracle. Edge additivity exact (`0.0`); oracle vs brute force 12 cases, 0 disagreements. **Scope: excludes crosstalk, correlated branch noise, route collisions, coherent inter-block errors.** |
| `PAIRLESS_BRANCH_IS_A_SEPARATE_COLUMN` | `NEW / STRUCTURAL` | It activates no pair, so it never idles, so it sits on a different reference template. Using the idle-free reference inside pricing broke N5 in 15/18 cases (errors to `0.71`). Collapses to one matching problem iff `idle_dephasing = 0`. |
| `NOISE_KILLS_ENTANGLEMENT_ABOVE_THRESHOLD` | `NEW / VERIFIED` | Between edge depolarization `0.15` and `0.20` the compiler abandons pairs entirely and the cost plateaus. The noise-aware objective genuinely changes the answer. |
| `SETTING_COST_AMORTIZATION` | `DONE` | `N* = (q_C−q_Q)c_setup/(c_Q−c̄_C)`. Worst crossover `13 047` shots at `c_setup = 10⁴` gate-equivalents; easy families cross at a few hundred. The K3 unpriced-settings hedge is replaced by this boundary. Cardinality cost kept **out** of the conic program. |
| `ZERO_CONTRAST_EDGE_CASE` | `FIXED` | `atan2(−c,d)` is undefined at zero contrast; now reported as `arbitrary_zero_information` with a null margin instead of a spurious `η_ro = 1`. |
| `ASYMMETRIC_READOUT_THREE_POINT` | `DONE` | Outcome-dependent confusion offsets the fringe; two points misread it. Three-point fit at `0, 2π/3, 4π/3` separates offset from harmonic. Gate N8. |
| `N7_DEPLOYABLE_READOUT` | `DONE` | Oracle analyzer angles are an upper bound. Deployable pilot arm retains `0.9932` (light) / `0.9751` (heavy) of it at the optimal pilot fraction, and exactly `1.0000` at unit visibility. Two separable losses: pilot-angle error (vanishes in pilot size) and off-quadrature pilot shots (∝ pilot fraction). |
| `N9_FIXED_SETTING_COST` | `FAILS` | Greedy forward selection is **not** exact: `7.24 %` excess at `λ_q = 0`, exact only once the cardinality penalty makes one setting optimal — i.e. outside the multi-setting regime. Exhaustive enumeration required on small instances; greedy labelled heuristic. |
| `N10_EQUAL_ACCOUNTING_CAMPAIGN` | `DONE / WEAKENED` | Rerun against QUEST-tE, whose preparation drops from 44 two-qubit gates to **3**. **QUEST now wins at 0, 1 % and 2 %**; CovQ wins from **5 %**, and the margin at 20 % is **1.39×, not 1.27e5×**. The qualitative regime statement survives; the quantitative one does not. CovQ still beats product-only everywhere, and is still scored on attained CFI against QUEST's unattained bound. Sparse-Carathéodory is included but its column reflects **this prototype's** 232-CX dense preparation, an upper bound, not the method. |
| `EQ_110_COMPLEX_ELEMENT_BUG` | `FIXED` | The Eq (110) contraction used `gs[j].conj().T`, which Hermiticity collapses to `gs[j][a,b]`, computing `Σ w·z²` rather than `Σ w·|z|²`. Correct only for real matrix elements — every earlier fixture, signed cat states included. Surfaced by the QUEST arm as a non-PSD QFIM. Fixed to `gs[j].T`; all previously reported numbers unchanged; tests now include random complex pure states and a PSD check. |
| `QUEST_UNDER_NOISE` | `RUN / UPPER_BOUND_ONLY` | Run in the conservative form. QUEST's column remains a **bound**, not an attainable number, until someone compiles a readout for it. Its preparation needs 44 two-qubit gates against a CovQ branch's ≤ 2, which is the entire mechanism of the separation. |
| `PROTOTYPE` | `PARTIAL` | `prototype/`. Built before this decision arrived; retained as the audit's evidence, not as Task 0C. See `prototype/README.md`. |

## Corrections carried into v0.2

Three claims from the first-pass review are wrong or overstated and are corrected in
`docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md`:

1. **The relaxed ℓ1-rank hardness result is a cone result.** For a convex decomposition
   `p ≥ 0`, `Σ p_r = 1`, so `‖p‖₁ = 1` identically and relaxed rank is vacuous on the
   convex hull. It was also attributed to a protocol arm that does not exist: v0.1 has
   Frank–Wolfe sparse approximation, not ℓ1-rank minimisation. These are different objects
   and must not be conflated.
2. **K1 was reported as fired; it is half-fired.** The specification collision with QUEST
   is real and the abstract must be rewritten. The end-to-end compiler collision — target
   plus hardware circuits plus *feasibility certificates* — is not established.
3. **`F_ii ≤ 1` was proposed as the fix for the four degeneracies; it is not safe.**
   Feasibility for a subunit diagonal is a joint first-and-second-moment condition
   `(μ, M) ∈ conv{(s, ss^T)}`, not `M ∈ Q_m`. The information-floor formulation achieves
   the same un-degeneration without leaving the audited sector.

A fourth correction is the prototype's own, and it retracts a number this repository
would otherwise have published:

4. **The "relaxed coherent-flag normal form is ~12× cheaper" measurement was not a
   like-for-like comparison.** The cheap construction buys its gate count by moving the
   information into flag–data coherence. Dephasing the flag sends its QFIM to zero, while
   the literal controlled-cat form is unaffected. They do not compile the same object.
   Gate `C2c` now separates them; see the table below.

## Measured (prototype, m = 4, seed 2026)

| program | cond. data width | flag qubits | Schmidt rank | CX | tr F coherent | tr F dephased | coherence load-bearing |
|---|---|---|---|---|---|---|---|
| `sign_purification` | 1 | 4 | 8 | 188 | 4.000 | 0.000 | **yes** |
| `flagged_cat_relaxed` | 1 | 3 | 8 | 50 | 4.000 | 0.000 | **yes** |
| `flagged_cat_literal` | 4 | 2 | 4 | 334 | 4.000 | 4.000 | no |

Read: the first two realise the target with conditional data width 1 — which is exactly why
width is not a valid resource for a general coherent-flag backend — and neither survives
flag dephasing. Only the third is interchangeable with a labelled schedule.

## Gates (prototype, m ∈ {3,4,5})

`C1 C2 C2c C3 C4 C5 C6 C7 C8 C9 C11` PASS · `C10 C12` MEASURED · 0 FAIL.

126 tests, ~34 s, numpy + scipy only.

`C10` and `C12` are reported as MEASURED by construction: a gate that cannot fail is not a
gate, and both of those state facts that hold for every matrix or every cat schedule.

## The one relationship worth quoting

The standard `k`-producibility bound is exactly the support function of `Q^prog_{m,k}` in the
single direction `b = 1`:

```
max { 1ᵀ F 1 : F ∈ Q^prog_{m,k} }  =  ⌊m/k⌋k² + r² ,     r = m − ⌊m/k⌋k
```

Verified exactly for every `(m, k)` with `m ≤ 7`, `k ≤ 3`. So the metrology literature's
witness is one direction; the characterisation gives all of them, and at `k = 2` the complete
facet description. This is the safe way to position the result against that literature, and
it doubles as an independent physics check on the construction.

## Not claimed

- No advantage over QUEST, generic moment matching, or any VQA. None has been run.
- No "first" of any kind.
- No novelty for the information-floor *interface* on its own: minimum-cost design under
  Fisher-information LMI constraints is established classically. Only the composite object is
  claimed.
- No strongly polynomial claim for the contract problem: the matching half inherits strong
  combinatorics, the contract half carries a semidefinite constraint.
- No new polyhedral combinatorics. Edmonds at `k = 2` and Partition into Triangles at
  `k ≥ 3` are both classical; the quantum content is the identification and the constructive
  circuit emission.
- No minimum-entanglement claim for the Clifford-conjugated generator arm: width there is
  frame-relative, and the prototype measures cases in both directions (a Z-frame width-3
  branch that is physically a product state, and Z-frame width-`m` branches that stay
  width-`m`).
- No gate-count comparison against `generic_full_state_preparation`; that baseline is
  deliberately weak and its counts are upper bounds.
