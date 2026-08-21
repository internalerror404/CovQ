# CovQ — STATUS

Last updated: 2026-08-21. Decision of record: `docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md`.
Current documents: `charter/*_v0.3.*`. Proof of record:
`docs/proofs/THEOREM_PAIR_WIDTH_AND_INFORMATION_FLOOR_v0.1.md`.

**The spine.** information requirement → exact pair-width geometry → certified minimum-cost
program → hardware-native Bell-pair schedule. The `k = 2` theorem is the engine; the
information-floor compiler is the reason anyone needs the engine.

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
| `QUEST_BASELINE` | `ABSENT` | Not run. No advantage over QUEST is claimable until it is. |
| `SPARSE_PURE_STATE_BASELINE` | `IMPLEMENTED` | Carathéodory-sparse single state, plus three realisation-matched siblings. |
| `MANUSCRIPT_INTEGRATION_MD` | `NOT_RECEIVED / BLOCKING` | The canonical v0.2 implementation contract (gates C1–C16, JSON record schema, status vocabulary, manifest and DONE rules, freeze-before-run) has not reached the repository. No manuscript table may be populated until it does. |
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
