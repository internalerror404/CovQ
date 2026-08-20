# CovQ Research Charter v0.3

**Supersedes** `archive/v0.2/` (drafted before the principal optimisation theorem was closed).
**Basis** `docs/proofs/THEOREM_PAIR_WIDTH_AND_INFORMATION_FLOOR_v0.1.md` (2A–2F closed),
`docs/audits/PRIOR_ART_SWEEP_TASK0B5_v0.1.md`, `docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md`,
`docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md`.

## Working title

**CovQ: Compiling Fisher-Information Contracts to Pair-Entangled Quantum Programs**
*Hardware Matching Polytopes, Certified Minimum-Cost Schedules, and Information Floors*

## Thesis

> CovQ compiles a Fisher-information contract: given a downstream information floor and a
> hardware graph, it emits the least-cost bounded-entanglement quantum program that
> certifiably satisfies the contract. For pair-entangled labelled programs the feasible set
> has an exact matching-polytope characterisation, which is what makes compilation
> constructive and certifiable.

One spine, four links:

```
information requirement → exact pair-width geometry → certified minimum-cost program
                                                    → hardware-native Bell-pair schedule
```

The `k = 2` theorem is the mathematical engine. The information-floor compiler is the reason
anyone needs the engine. Neither is presented without the other.

## The compiler

```
(𝒢, H, A, G_req, 𝒞) ↦ Π ,        min_Π Cost_𝒞(Π)   s.t.   Aᵀ F_Π A ⪰ G_req
```

`𝒢` independent commuting Pauli generators · `H` the native interaction graph ·
`A ∈ R^{m×d}` the downstream parameter combinations · `G_req ⪰ 0` the required Fisher floor
in downstream space · `𝒞` a declared resource and noise cost model · `Π` a **labelled**
circuit schedule.

Exact matching `F_Π = F_⋆` remains as a compiler *mode* and a verification benchmark. It
stops being the reason the compiler exists. Four things forced that change: prescribing an
exact unit-diagonal QFIM collapses Loewner dominance to equality; a referee reasonably asks
why a user wants one particular QFIM rather than sufficient information for a task; exact
matching obscures the compiler's actual freedom; and generic expectation targeting already
occupies most of the "partial state specification" motivation.

## Contributions

**C1 — Certified information-floor compilation (Theorem 2).** For non-negative native edge
costs the contract problem is the convex program

```
min Σ_e c_e t_e   s.t.  AᵀFA ⪰ G_req,  diag F = 1,  F_ij = 0 (ij ∉ E),
                        −t_ij ≤ F_ij ≤ t_ij,  t ∈ MATCH(H)
```

whose objective **is** the expected per-shot native Bell-pair bill, because `MATCH(H)` is
down-closed and so `t_e = |F_e|` at any optimum. Solved by LP over the matching polytope plus
Loewner cutting planes: convex, certificate-producing optimisation to prescribed precision.
*Not* a strongly polynomial algorithm — the matching half inherits strong combinatorics, the
contract half carries a semidefinite constraint, and the paper says so.

**C2 — Exact pair-width geometry (Theorem 1), the engine.**

```
F ∈ Q^lab_{H,2}  ⟺  diag F = 1,  F_ij = 0 for ij ∉ E,  (|F_ij|)_{ij∈E} ∈ MATCH(H)
```

Constructive, and every branch has **two-qubit depth exactly 1** — a matching is by
definition a set of gates that all run in parallel. Polynomial-time decidable by Edmonds
separation, against NP-hard membership at unbounded width. The quantum step is that signs are
independently realisable on a matching, so magnitudes alone decide feasibility; the remark in
2B shows this fails at width 3 by a sign-parity obstruction, so it is a genuine feature of
the pair level rather than the first case of a pattern.

**C3 — Certificates.** Three outcomes, each independently checkable: a feasible `F` with its
matching decomposition, emitted circuits and satisfaction margin; a lower bound on achievable
cost from any valid relaxation; or an infeasibility certificate for the requested floor. This
was one of the original white-space claims and it is now discharged.

**C4 — The odd/even collective-mode ceiling.** For `u_S = 1_S/√k`, width-2 caps `u_SᵀFu_S` at
`2` for even `k` and `2 − 1/k` for odd `k`. Edmonds' blossom constraint restated as physics:
an odd collective mode always leaves one qubit unmatched and provably cannot reach the
pair-entanglement ceiling. This is the operational meaning of branch width.

**C5 — Relation to the `k`-producibility witnesses.** The standard entanglement-depth bound
`⌊m/k⌋k² + r²` is the support function of the reachable body in one direction. We give the
complete facet description at `k = 2`. Witness theory asks whether a measured scalar rules out
small depth; CovQ asks whether an entire requested contract is achievable at small width, and
emits the program.

**C6 — The complexity boundary.** Width 2 reduces to maximum-weight matching; width 3
conjecturally NP-hard via Partition into Triangles. Its role is to explain why the pair-width
compiler is a *privileged tractable regime*, not to claim a complexity result.

**C7 — Conditioning-aware downstream stability** on the identifiable quotient, and negative
results: regimes where product or small-width programs are preferable, and where the generator
frame rather than the target sets the entanglement cost.

## Background propositions (correct, not contributions)

BP1 feasibility `= conv{ssᵀ}` (Pitowsky; Huber–Marić plus the amplitude realisation).
BP2 the signed-cat primitive. BP3 labelled additivity `F_Π = Σ_r p_r F_r`, a discrete optimal
design. BP4 submultiplicativity, replaced by C7. BP5 membership and decomposition-rank
hardness (Caprara et al.; their relaxed ℓ1-rank result is a *cone* result, vacuous on the
convex hull where `‖p‖₁ ≡ 1`). BP6 attainability: `Im⟨∂_iψ|∂_jψ⟩ = 0` identically for
commuting Hermitian generators — local at the operating point, POVM may depend on it, singular
targets restricted to their identifiable support, and hardware runs report an attainable
*classical* Fisher matrix unless the measurement is compiled.

## Program semantics

**Labelled schedule** — the sole subject of the bounded-width theorem. Per-shot
`F_Π = Σ_r p_r F_r`; width is a valid resource.

**Coherent accessible flag** — a separate resource model, not an interchangeable backend.
`F_Ψ = Σ_r p_r F_r + Cov_p(μ_r)`, equal to the labelled schedule iff branch means *coincide*.
Conditional branch width is **not** a resource here: `|Ψ_p⟩ = Σ_s √p(s)|s⟩_f|s⟩_d` realises any
feasible target with every conditional branch a product eigenstate, while dephasing the flag
sends the QFIM to zero. Report flag dimension, Schmidt rank across the flag|data cut,
controlled-preparation cost and joint readout instead. Programs on opposite sides of the
dephasing test are **not comparable on gate count**.

**Unlabelled mixed state** — forbidden as a substitute in general, though measured to cost
*nothing* for signed-cat schedules, since `Z_i|C_s⟩ = s_i|C⁻_s⟩` maps the support into the
kernel.

## Resource model

Reported separately, never merged: settings · `k_logical` · `k_physical` · Clifford two-qubit
cost · two-qubit gate count and depth · SWAP count · ancilla count · flag dimension and
Schmidt rank for coherent programs · expected per-shot native edge cost `Σ_e c_e|F_e|` ·
predicted noisy QFIM error.

## Novelty boundary

**Occupied.** Expectation-value-constrained synthesis (QUEST). Correlation, cut, matching and
bounded-cluster partition polytopes. QFIM entanglement witnessing and depth criteria. Optimal
experimental design over additive Fisher information, including Fedorov–Wynn. **And — found
in the supplementary sweep — minimum-cost sensor selection and placement with Fisher-information
LMI constraints, i.e. information as a constraint rather than an objective, together with
Loewner-monotone design criteria and Loewner-constrained spectral design.** The
information-floor *interface* therefore has a clear classical shadow and must not be presented
as novel on its own.

**Claimed — the composite object, and only that.**

> We found no prior work that combines an information-floor compiler interface with an exact
> bounded-entanglement QFIM achievable set and constructive hardware-native circuit emission.

**Not claimed.** Any "first". Any advantage over QUEST or any VQA until one is run. The
discovery of matching or clique-partitioning polytopes. Any new polyhedral combinatorics —
Edmonds at `k = 2` and Partition into Triangles at `k ≥ 3` are classical, and the sufficiency
decomposition is off-the-shelf. Any minimum-entanglement statement for the Clifford-conjugated
arm. Any strongly polynomial claim for the contract problem.

**Two presentation hazards, pre-empted in the introduction rather than after the theorem.**
(i) Matchings are the standard model for *parallel two-qubit gate scheduling*, where a matching
is a time slice inside one circuit; here it is a *branch of a randomised schedule*. (ii) The
clique-partitioning polytope is over 0/1 transitive edge indicators; `Q^lab_{H,2}` is over
signed correlations in `[−1,1]` and equals the correlation polytope at `k = m`.

## Kill criteria

K1a specification collision with QUEST — **fired**; abstract rewritten. K1b end-to-end
compiler collision *with equivalent feasibility certificates* — not established, still live.
K2 the width-2 theorem is a known matching statement with no quantum content — **not fired,
conditional on the identification-plus-achievability framing**. K3 no material measured
advantage over QUEST. K4 no useful bounded-width tradeoff — **cleared**. K6 the
width-2/width-3 boundary is false or already known and nothing replaces it — live, and
"already known" is now the likelier branch. K7 only spacetime targets work. K8 only a VQA
result survives. **K9 (new): the information-floor interface adds nothing over classical
Fisher-constrained design once the quantum feasible set is fixed** — i.e. if every benchmark
contract is met by the product probe or by a fixed state, the composite claim is empty.

## Benchmarks that make the floor non-trivial

The floor formulation must avoid contracts that `F = I` already satisfies.

1. **Common mode.** `A = 1/√m`, `G_req = [γ]`: closed-form curve `cost*(γ) = m(γ−1)/2` for
   `1 ≤ γ ≤ 1 + 2⌊m/2⌋/m`. Product probe reaches only `γ = 1`.
2. **Odd/even collective modes.** Ceiling `2` vs `2 − 1/k` — a sharp, analytic, blossom-driven
   separation.
3. **Overlapping multi-direction contracts.** Two 3-modes sharing a qubit bind jointly at
   `γ* ≈ 1.3304` against `5/3` for one alone or two disjoint — the compiler must balance
   directions instead of concentrating information in one GHZ-like mode.
4. Infeasible-with-certificate, boundary, structured, Clifford-transformed, and one
   application contract.

## Application mode (one late benchmark)

Paper 1 separates the forward map from the programmable probe information, `F_θ = Jᵀ W_ρ J`.
CovQ's application mode is `F_Π = W_ρ`, `A = J`, contract `Jᵀ F_Π J ⪰ G_{θ,req}`: given a
required floor on the identifiable spacetime parameters, compile the cheapest probe program
whose link-level QFIM meets it. More natural than extracting one exact `W_⋆`, because many
`W_ρ` meet the same metric-level requirement and choosing among them is the compiler's whole
purpose. The paper stays quantum-computing-only because `A` is generic.

## Paper structure

1. Introduction — Fisher-information contracts, not target states
2. Related work and novelty boundary (optimal design, Fisher-constrained sensor design,
   matching polytopes, gate scheduling, QFIM witnesses)
3. Commuting-Pauli programs, the Fisher IR, and the resource vector
4. Background: feasibility, signed cats, labelled additivity, attainability
5. **Exact pair-width geometry on a hardware graph** (the engine)
6. **Certified information-floor compilation** (the contract, and the certificates)
7. Operational consequences: the odd/even ceiling, the common-mode curve, `k`-producibility
8. Complexity: width two versus width three
9. Exact, sparse, and approximate algorithms
10. Generator-aware Clifford normalisation and frame-relative width
11. Routing as a priced backend
12. Coherent flags as a separate resource model
13. Conditioning-aware stability on the identifiable quotient
14. Experimental methodology
15. Results
16. Application benchmark from Lorentzian tomography
17. Limitations and open problems
18. Conclusion

## v1 scope exclusions

Noncommuting generators. Arbitrary mixed-state QFI. Universal speedup claims. QAOA as primary
solver. Fault-tolerant synthesis. Spacetime as the main object. **Subunit-diagonal realised
probes** — feasibility there is the joint condition `(μ, M) ∈ conv{(s, ssᵀ)}`, not `M ∈ Q_m`,
and the prescribed-moment literature is unaudited; the information floor delivers the same
un-degeneration without leaving the audited sector.
