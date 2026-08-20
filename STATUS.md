# CovQ — STATUS

Last updated: 2026-08-20. Decision of record: `docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md`.

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
| `INFORMATION_FLOOR_SEMANTICS` | `ADOPTED` | `min Cost(Π) s.t. F_Π ⪰ G_⋆`. Replaces exact matching as the compiler's primary mode. |
| `SUBUNIT_DIAGONAL_SECTOR` | `REJECTED_FOR_V1` | Needs joint `(μ, M)` moment feasibility, not `M ∈ Q_m`. Not audited. |
| `K1_SPLIT` | `DONE` | K1a specification collision: **FIRED**. K1b end-to-end compiler collision: **NOT ESTABLISHED**. |
| `QUEST_BASELINE` | `ABSENT` | Not run. No advantage over QUEST is claimable until it is. |
| `SPARSE_PURE_STATE_BASELINE` | `IMPLEMENTED` | Carathéodory-sparse single state, plus three realisation-matched siblings. |
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

65 tests, ~2 s, numpy + scipy only.

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
- No new polyhedral combinatorics. Edmonds at `k = 2` and Partition into Triangles at
  `k ≥ 3` are both classical; the quantum content is the identification and the constructive
  circuit emission.
- No minimum-entanglement claim for the Clifford-conjugated generator arm: width there is
  frame-relative, and the prototype measures cases in both directions (a Z-frame width-3
  branch that is physically a product state, and Z-frame width-`m` branches that stay
  width-`m`).
- No gate-count comparison against `generic_full_state_preparation`; that baseline is
  deliberately weak and its counts are upper bounds.
