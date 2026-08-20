# CovQ Research Charter v0.2

**Supersedes** `archive/v0.1/CovQ_QFIC_Research_Charter_v0.1.md`.
**Basis** `docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md`, `docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md`.

## Working title

**CovQ: Pair-Entangled Compilation of Commuting-Pauli Information Requirements**
*Hardware Matching Polytopes, Labelled Schedules, and Certified Information Floors*

## Thesis

CovQ compiles commuting-Pauli **information requirements** into backend-specific quantum
programs under explicit entanglement, coherence, settings, and hardware resources. Its first
exact regime is **labelled pair-entangled compilation**, characterised by a signed hardware
matching polytope.

The compiler receives an information floor rather than a state, a unitary, or an exact
covariance:

```
(𝒢, G⋆, H, ε, 𝒞)  ↦  Π ,        min Cost(Π)  s.t.  F_Π ⪰ G⋆ .
```

`𝒢` is an independent commuting Pauli generator program, `G⋆ ⪰ 0` a required information
matrix with possibly subunit diagonal, `H` a hardware coupling graph, `ε` a tolerance, `𝒞` a
declared resource and noise cost model, and `Π` a **labelled** circuit schedule.

The motivating question is not "why would anyone prescribe a covariance matrix". It is:
a user declares the minimum information downstream estimation requires, and the compiler
returns the cheapest physical QFIM that dominates it.

Exact matching, `F_Π = F⋆`, is retained as a verification mode and a benchmark family. It is
no longer the compiler's semantics, because inside a unit-diagonal target sector
`F ⪰ F⋆` forces `F = F⋆` and the three notions collapse into one.

Spacetime tomography is one late benchmark. The paper must remain coherent if all spacetime
material is deleted.

## Formal setting

`U_θ = exp[−(i/2) Σ_i θ_i P_i]` with `P_i² = I` commuting and independent. For a pure probe,
`F_ij = Re⟨P_i P_j⟩ − ⟨P_i⟩⟨P_j⟩`, hence `F_ii = 1 − ⟨P_i⟩²`: unit diagonal and zero mean are
the same condition. Realised probes are kept unit-diagonal; the *requirement* `G⋆` is not.

After simultaneous Clifford diagonalisation to `P_i = Z_i`, a zero-mean pure state
`Σ_z √p_z |z⟩` has `F = E_p[zz^T]`.

## Contributions, in the order they will be defended

**C1 — the width-two characterisation (the spine).** For the labelled, zero-mean,
physical-Z backend,

```
F ∈ Q^lab_{m,2}   ⟺   (|F_ij|)_{i<j} ∈ MATCH(K_m)
```

and, hardware-natively, with no routing inside the primitive backend,

```
Q^lab_{H,2} = { F : diag F = 1, F_ij = 0 for ij ∉ E, (|F_ij|)_{ij∈E} ∈ MATCH(H) } .
```

Sufficiency is constructive: it emits a schedule of parallel Bell pairs on physical edges.
Feasibility is polynomial-time by Edmonds separation, against NP-hard membership at `k = m`,
with explicit separating certificates and at most `O(m²)` branches by Carathéodory.

**C1a — corollary tying the result to the metrology literature.** The standard
`k`-producibility bound is the support function of `Q^prog_{m,k}` at `b = 1`:
`max{1ᵀF1 : F ∈ Q^prog_{m,k}} = ⌊m/k⌋k² + r²`, verified exactly for `m ≤ 7`, `k ≤ 3`, and
equal to `m + 2⌊m/2⌋` in closed form at `k = 2`. The characterisation gives every other
direction, and at `k = 2` the complete facet description.

**C2 — the complexity boundary.** *Currently unproved, and expected to be a corollary rather
than a theorem.* Width 2 reduces to maximum-weight matching (Edmonds). Target: width 3 linear
utility optimisation is NP-hard, via Partition into Triangles — which is classically
NP-complete, so the reduction should go through in a paragraph. The live risk is therefore
not that it is false but that it is **already folklore** in the bounded-cluster partitioning
literature. Settle that before claiming it.

**C3 — routing as a priced backend.** Routing enlarges the admissible edge set at an
explicit SWAP/path cost. Reported against `Q^lab_{H,2}` on the native graph, never folded
into it.

**C4 — the information-floor compiler.** `min Cost(Π)` s.t. `F_Π ⪰ G⋆`. For the width-2
hardware-native backend this is a matching-polytope-constrained semidefinite problem with
certificates.

**C5 — coherent-flag semantics as a *separate* resource model.** Not an interchangeable
backend. See below.

**C6 — conditioning-aware downstream stability**, on the identifiable quotient.

**C7 — negative results and regime boundaries**, including where product or small-width
programs are preferable and where the frame choice, not the target, sets the entanglement
cost.

## Background propositions (correct, not contributions)

- **BP1 (was Theorem A).** By symmetric Bernoulli correlation-polytope theory plus the
  pure-state amplitude realisation, the feasible zero-mean unit-diagonal QFIMs are exactly
  `conv{ss^T}`. Inherited from Pitowsky and from Huber–Marić.
- **BP2 (was the signed-cat proposition).** `|C_s⟩ = (|s⟩+|−s⟩)/√2` realises `F = ss^T`.
  A compiler primitive.
- **BP3 (labelled additivity).** Fisher information is additive over independent runs with
  known labels, so `F_Π = Σ_r p_r F_r` per shot. A semantics lemma and a unit test. It is a
  discrete optimal experimental design; the branch circuits are design points and the branch
  probabilities the design measure.
- **BP4 (was Theorem D).** `‖AᵀΔFA‖ ≤ ‖A‖²‖ΔF‖` is submultiplicativity. Replaced by C6.
- **BP5 (was Theorem E).** Membership and decomposition-rank hardness over correlation
  polytopes is Caprara et al. Cited, not proved. Note their relaxed ℓ1-rank result is a
  *cone* result and is vacuous on the convex hull, where `‖p‖₁ = 1` identically.
- **BP6 (attainability).** `Im⟨∂_iψ|∂_jψ⟩ = −(1/4) Im⟨ψ|P_iP_j|ψ⟩ = 0` identically for
  commuting Hermitian generators. Load-bearing for every "certified performance" phrase and
  absent from v0.1. Stated with its qualifications: attainability is local at the operating
  point; the optimal POVM may depend on it; singular targets are restricted to their
  identifiable support; labelled programs permit branch-conditioned measurements while
  coherent programs need a joint flag–data measurement; hardware runs report an explicit
  attainable *classical* Fisher matrix unless the measurement is actually compiled.

## Program semantics

**Labelled schedule.** Sample branch `r` with probability `p_r`, run `U_r`, retain `r` in the
classical record. Per-shot QFIM `Σ_r p_r F_r`. **This is the sole subject of the
bounded-width theorem.**

**Coherent accessible flag.** `|Ψ⟩ = Σ_r √p_r |r⟩_f |ψ_r⟩_d`, generators on data only. The
exact identity is

```
F_Ψ = Σ_r p_r F_r + Cov_p(μ_r) ,
```

with equality to the labelled schedule **iff the branch means coincide**, not merely vanish.

**The width loophole, and why coherent flags are a separate model.** For any feasible `F`,
`|Ψ_p⟩ = Σ_s √p(s)|s⟩_f|s⟩_d` realises it with every conditional data branch a *product
eigenstate* — conditional width 1 — while dephasing the flag sends the QFIM to zero. Nothing
became cheap; the resource moved into flag–data entanglement, flag dimension, Schmidt rank
across the flag|data cut, controlled-preparation cost, and joint readout. Therefore:

> Conditional data-branch width is **not** a valid resource measure for a general
> coherent-flag backend, and one width number must not be shared across the two backends.

Coherent-flag programs are reported with their own resource vector, and gate C2c decides
which side of the dephasing test a program is on. Programs on opposite sides are **not
comparable on gate count**.

**Unlabelled mixed state.** Still not interchangeable in general. Measured finding: for a
signed-cat schedule the gap is *exactly zero*, because `Z_i|C_s⟩ = s_i|C⁻_s⟩` maps the
support into the kernel. The gap is real elsewhere — `(|C_s⟩, |C⁻_s⟩)` at equal weight has
covariance `ss^T` and QFI zero.

## Resource model

Reported separately, never merged:

`k_logical` (branch width in the Z frame) · `k_physical` (product width of the emitted
physical branch) · Clifford two-qubit cost · branch count / settings · two-qubit gate count
and depth · SWAP count · ancilla count · flag dimension and Schmidt rank for coherent
programs · predicted noisy QFIM error.

A globally entangling Clifford turns a width-2 normalised state into a globally entangled
physical state, so "minimum physical entanglement" is valid for the physical-Z arm and not
automatically for the Clifford-conjugated arm. Measured cases run in **both** directions,
which makes the frame choice a compiler degree of freedom rather than only a caveat.

## Novelty boundary

**Occupied.** Expectation-value-constrained pure-state synthesis (QUEST). Partial /
don't-care state preparation. General state preparation and unitary synthesis. Cat, graph
and stabiliser synthesis. Bernoulli correlation and cut polytopes. QFIM entanglement
witnessing and dimensionality. Approximate Carathéodory. Matching and size-constrained
clique-partitioning polytopes. Optimal experimental design over additive Fisher information,
including Fedorov–Wynn exchange algorithms.

**Claimed.** The identification of the labelled commuting-Pauli QFIM program class with a
signed *hardware* matching polytope; exact quantum achievability with emitted circuits;
the width-two tractable boundary; feasibility and infeasibility certificates connected to
emitted programs; the information-floor compiler mode.

**Not claimed.** Any "first". Any advantage over QUEST or any VQA until one is run. The
discovery of matching or clique-partitioning polytopes. Any minimum-entanglement statement
for the Clifford-conjugated arm. **Any new polyhedral combinatorics**: Edmonds at `k = 2` and
Partition into Triangles at `k ≥ 3` are both classical, and the sufficiency direction uses an
off-the-shelf strongly-polynomial matching decomposition. The quantum content is the
identification of the QFIM feasible set with that structure, the exact achievability, and the
constructive circuit emission. The paper says so in its own introduction rather than letting a
referee say it first.

**Two presentation hazards to pre-empt in the introduction.** (i) Matchings are the standard
model for *parallel two-qubit gate scheduling*, where a matching is a time slice inside one
circuit; here a matching is a *branch of a randomised schedule*. (ii) The clique-partitioning
polytope is over 0/1 transitive edge indicators; `Q^prog_{m,k}` is over signed correlations in
`[-1,1]` and equals the correlation polytope at `k = m`. Both distinctions belong before the
theorem, not after it.

## Kill criteria

K1a specification collision — **already fired**; the abstract is rewritten accordingly.
K1b end-to-end compiler collision, i.e. prior work accepting the target and emitting hardware
circuits *with equivalent feasibility certificates* — not established; still live.
K2 the width-two theorem turns out to be a known statement about matching polytopes in
disguise with no quantum content.
K3 no material advantage over QUEST or generic moment matching on any declared resource,
measured.
K4 the width hierarchy yields no useful or provable resource tradeoff. *(Cleared by C1.)*
K5 freedom over the QFIM-equivalence class never beats compiling one fixed state.
K6 the width-2/width-3 boundary is false or already known, **and** nothing replaces it.
   *Sweep finding: "already known" is now the likelier branch. Step 3 of Task 0B.5 is
   re-scoped from "prove or abandon" to "prove, then establish whether it is folklore".*
K7 only spacetime-derived targets work.
K8 only a VQA result survives.

## Contribution checklist

At least: one exact feasibility and normal-form theorem on the hardware graph; one
complexity separation or a stated, defended failure to obtain one; one exact certified
algorithm; one scalable approximation algorithm with certificates; emitted hardware circuits
with measured resource accounting; generic, structured, infeasible, boundary and application
benchmark families; four realisation-matched baselines plus QUEST; and a negative result or
regime boundary.

## Working abstract v0.2

> Quantum compilers receive a target unitary, statevector, or circuit. We study a different
> specification: an *information requirement* for a commuting-Pauli parameter program, given
> as a positive semidefinite matrix that the realised quantum Fisher information matrix must
> dominate. The compiler may choose any labelled schedule of circuits whose per-shot Fisher
> information meets the requirement, and it is scored on settings, entanglement width,
> two-qubit gates and depth, routing, and ancillas. Our main result characterises exactly
> which requirements are met by *pair-entangled* schedules on a given hardware graph: a
> target is realisable by branches of entanglement width two on `H` precisely when its matrix
> of off-diagonal magnitudes, supported on the edges of `H`, lies in the matching polytope of
> `H`. The characterisation is constructive — it emits a schedule of parallel Bell pairs on
> physical edges — and decidable in polynomial time by Edmonds separation, in contrast to the
> NP-hardness of membership at unbounded width. Violated blossom inequalities give explicit
> certificates that a requirement provably needs more than pair entanglement, and we exhibit
> families of perfectly realisable Fisher matrices that require width three for exactly that
> reason. We give exact and certified-approximate algorithms, map programs to restricted
> hardware graphs under gate, depth and routing costs, and compare against four
> realisation-matched single-state and flagged baselines. We separate labelled schedules from
> coherently flagged programs by an explicit dephasing test, and show that conditional branch
> width is not a resource measure for the latter.

## Paper structure

1. Introduction — information requirements, not states
2. Related work and novelty boundary (incl. optimal design, matching polytopes)
3. Commuting-Pauli programs, the Fisher IR, and the resource vector
4. Background: feasibility, signed cats, labelled additivity, attainability
5. **Pair-entangled compilation and the hardware matching polytope**
6. **Complexity: width two versus width three**
7. Certificates: blossom obstructions and uncertainty relations
8. Exact, sparse, and approximate algorithms
9. Generator-aware Clifford normalisation and frame-relative width
10. Routing as a priced backend
11. Coherent flags as a separate resource model
12. Information floors and conditioning-aware stability
13. Experimental methodology
14. Results
15. Application benchmark from Lorentzian tomography
16. Limitations and open problems
17. Conclusion

## v1 scope exclusions

Noncommuting generator QFIMs. Arbitrary mixed-state QFI compilation. Universal speedup
claims. QAOA as primary solver. Fault-tolerant synthesis. Spacetime as the main object.
**Subunit-diagonal realised probes**: feasibility there is a joint `(μ, M) ∈ conv{(s, ssᵀ)}`
condition, not `M ∈ Q_m`, and the prescribed-moment literature has not been audited. The
information floor delivers the same un-degeneration without leaving the audited sector.
