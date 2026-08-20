# DECISION 001 — Pivot to labelled pair-width compilation

**Date:** 2026-08-20
**Status:** Adopted
**Supersedes:** `archive/v0.1/` (charter, protocol, agent prompt, paper skeleton)
**Trigger:** Task 0A/0B audit, `docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md`

## Decision

CovQ pivots. It does not terminate.

| | |
|---|---|
| Task 0C v0.1 | **STOP** — recorded `SUPERSEDED`, not silently modified |
| CovQ as a paper | **PIVOT** |
| Theorems A / B / E as headline contributions | **RETIRE** to background |
| Bounded-width compilation | **CONTINUE** |
| The `k = 2` matching characterisation | **PROMOTE** to the mathematical spine |
| Coherent flags | **REDEFINE** as a separate resource model, not an interchangeable backend |

This is the charter's own stop rule operating as designed. v0.1 required a pivot if no
theorem survived beyond the correlation-polytope translation, and separately required a
useful bounded-width result. The first condition fired; the second is where the surviving
result lives.

## Why each line

**Retire A, B, E.** Theorem A is Pitowsky's correlation polytope plus the observation that
commuting-Z QFI depends only on computational-basis probabilities, so
`p(s) ↦ Σ_s √p(s)|s⟩` transports the classical characterisation verbatim. Huber and Marić
already give the polytope and its vertices. Theorem E's membership and decomposition-rank
hardness is Caprara et al. Theorem D is submultiplicativity. Each is true; none is a
contribution.

**Continue bounded width.** It is the one axis on which v0.1 asked a question that the
prior art does not answer, and the `k = 2` level turns out to be exactly solvable.

**Promote `k = 2`.** For the labelled, zero-mean, physical-Z backend:

> `F ∈ Q^lab_{m,2}` iff `(|F_ij|)_{i<j} ∈ MATCH(K_m)`.

Both directions are proved in the audit and the sufficiency direction is constructive: it
emits the schedule. The hardware-native form is stronger and is what the compiler actually
needs:

> `Q^lab_{H,2} = { F : diag F = 1, F_ij = 0 for ij ∉ E, (|F_ij|)_{ij∈E} ∈ MATCH(H) }`.

This maps a target QFIM onto a schedule of parallel Bell pairs on physical edges, with
routing priced separately as a backend that enlarges the admissible edge set. That is
recognisably a compiler theorem rather than a relabelling of Bernoulli correlations.

**Redefine coherent flags.** The exact identity is
`F_Ψ = Σ_r p_r F_r + Cov_p(μ_r)`, so equality with the labelled schedule needs the branch
means to *coincide*, not to vanish. Worse, conditional branch width is not a resource for
coherent flags at all: `|Ψ_p⟩ = Σ_s √p(s)|s⟩_f|s⟩_d` realises **any** feasible target with
every conditional data branch a product eigenstate. Nothing became cheap; the cost moved
into flag–data entanglement, flag dimension, Schmidt rank, and joint readout. One width
number cannot be shared across the two backends.

Chosen treatment: **labelled programs are the sole subject of the bounded-width theorem.**
Coherent flags keep a semantics section with their own resource vector. The prototype's
gate `C2c` decides which side a given program is on by dephasing its flag.

## Consequences

1. The target sector stays unit-diagonal for the *realised* probe. The `F_ii ≤ 1` proposal
   is rejected for v1: feasibility there is a joint `(μ, M) ∈ conv{(s, ss^T)}` condition,
   not `M ∈ Q_m`, and the prescribed-moment literature has not been audited.
2. Compiler semantics become an **information floor**:
   `min_Π Cost(Π)` subject to `F_Π ⪰ G_⋆`, with `G_⋆` positive semidefinite and possibly
   subunit-diagonal. This un-degenerates Loewner dominance without leaving the audited
   sector, and it answers the motivation question a referee asks first: a user declares the
   minimum information downstream estimation needs, and the compiler finds the cheapest
   physical QFIM that dominates it. Exact matching survives as a verification mode and a
   benchmark family, not as the only semantics.
3. Theorem D is replaced by conditioning-aware stability on the identifiable quotient,
   reporting rank, `λ⁺_min`, `κ⁺`, kernel principal angles, and `‖F⁺ − F⋆⁺‖`.
4. Optimal experimental design is **embraced in related work**, not hidden. A labelled
   schedule is a discrete design; Frank–Wolfe is Fedorov–Wynn. CovQ supplies the quantum
   design space and the exact bounded-entanglement feasible sets, not the optimiser.
5. Novelty gate K1 splits into K1a (specification collision with QUEST — **fired**, rewrite
   the abstract) and K1b (end-to-end compiler collision — **not established**).

## Open, and required before the strongest version of the paper

The `k = 2` theorem clears the charter's K4 gate. It does not by itself carry a major
quantum-computing paper, because the underlying combinatorics are classical matching
theory. The paper needs a complexity transition:

- width 2: linear optimisation reduces to maximum-weight matching (Edmonds, polynomial);
- width 3: linear utility optimisation NP-hard.

Candidate reduction: Partition into Triangles. Give graph edges utility 1 and non-edges 0;
a width-3 partition of `3q` vertices achieves utility `3q` exactly when the vertices
partition into `q` triangles. This needs a formal reduction and a careful comparison
against the size-constrained clique-partitioning polytope literature.

**Safe novelty framing:** we identify the commuting-Pauli QFIM program class with an
established bounded-cluster polyhedral structure, prove exact quantum achievability with
emitted circuits, and locate the width-two tractable boundary. **Not:** we discovered
matching or clique-partitioning polytopes.

## Execution order

Task 0B.5, in this order, before any further prototype work:

1. Prior-art sweep: matching polytopes, size-constrained clique partitioning, QFIM
   k-producibility bounds, hardware-native Bell-pair scheduling.
2. Formal proof of the `k = 2` theorem for both `K_m` and a hardware graph `H`.
3. Prove or abandon the width-2/width-3 boundary.
4. Resolve coherent-flag resource semantics (decision above; write it up).
5. Rewrite charter, protocol, paper abstract. *(Draft in `charter/`, this commit.)*
6. Only then, the minimal prototype.

## Deviation recorded

Step 6 was partially executed before this decision arrived, and is retained rather than
discarded because it supplies the audit's numerical evidence — including the brute-force
verification of the `k = 2` theorem, the demonstration that Edmonds' odd-set facets are
*active* rather than redundant, and the measurement that retracts a claim the first-pass
review would otherwise have carried into the charter. It is labelled `PARTIAL`, is not
Task 0C, and does not close step 6. See `prototype/README.md`.
