# Task 0B.5 step 1 — prior-art sweep

**Question this sweep exists to answer:** is the `k = 2` theorem a contribution or a
translation? It has to be settled before the formal proof is written (step 2), because the
answer changes what the proof is *for*.

**Verdict: K2 does not fire — conditionally.** The theorem is not a known statement in
disguise. But its combinatorial half is entirely classical, so it survives only if presented
as *identification plus quantum achievability*, never as new polyhedral combinatorics. The
framing sentence in `DECISION_001` is the one that has to be used.

**New risk surfaced, on K6.** See A2 and §3: the width-3 hardness is probably provable, and
probably a short corollary of a classical NP-complete problem rather than a theorem. That is
still worth stating, but it must not be claimed as new.

---

## A1 — Matching polytopes and Edmonds' theorem

**Status: fully classical. Used as-is.**

Edmonds' description — degree constraints, odd-set (blossom) constraints, non-negativity —
is standard, with polynomial-time optimisation and separation. More pointedly for us: every
point of the matching polytope decomposes into a convex combination of integral matchings,
and there are *strongly polynomial* algorithms that construct the decomposition by
maintaining laminar families of tight odd cuts. Bipartite case is Birkhoff–von Neumann;
general graphs are Edmonds.

**Consequence.** The sufficiency direction of our theorem is not an algorithmic contribution.
The decomposition it needs is off-the-shelf. What the theorem contributes is that the QFIM
feasibility question *is* that decomposition question, and that each matching in the
decomposition is realised by an explicit parallel-Bell-pair circuit.

**Risk: MEDIUM-HIGH as a presentation hazard.** Any claim phrased as "we characterise a new
polytope" gets rejected on sight. Phrase it as: the feasible set is the matching polytope
pulled back through the absolute-value map, and the pullback is exactly achievable.

## A2 — Size-constrained clique partitioning

**Status: adjacent, not the same object.**

The clique partitioning polytope is the convex hull of characteristic vectors of clique
partitions — **0/1 edge indicators subject to transitivity**. Its facet theory (Grötschel–
Wakabayashi, Oosten et al., chorded-cycle and 𝑘-cycle facets, and 2024–2026 work on new
facets and on redundancy of transitivity constraints) is a mature literature, as is the
bounded-clique-number variant.

`Q^prog_{m,k}` is **not** that polytope:

| | clique partitioning | `Q^prog_{m,k}` |
|---|---|---|
| variables | 0/1 edge indicators | signed correlations in `[-1, 1]` |
| constraint | transitivity | positive semidefiniteness + block realisability |
| at `k = m` | the whole clique | the **correlation polytope** `Q_m`, not a partition polytope |
| vertices | partitions | `⊕_B s_B s_B^T` over partitions **and sign patterns** |

The overlap is real but confined: "partition into cliques of size ≤ 2" *is* a matching, so at
our level the combinatorial skeleton is the classical one. The step that is not inherited is
the reduction to **magnitudes** — that the signs impose no constraint beyond `(|F_ij|)`
lying in the matching polytope. That is what has to be proved, and it is what makes the
characterisation an *iff* rather than a relaxation.

**Risk: MEDIUM.** Mitigation: a paragraph in §2 that states the table above explicitly, and a
citation to the bounded-cluster partition-polytope line so no referee thinks it was missed.

## A3 — QFIM `k`-producibility bounds

**Status: our object is a strict refinement of theirs, and this is an asset.**

The metrology literature gives two things. First, **scalar** bounds: for a `k`-producible
state of `N` particles the QFI in a collective direction is at most
`⌊N/k⌋k² + r²` with `r = N − ⌊N/k⌋k`, which is the standard entanglement-depth witness.
Second, **matrix witness lower bounds** — a moment matrix that lower-bounds the quantum
Fisher matrix, used for *detection* of `k`-separability.

Neither characterises the **achievable set** of QFI matrices under a width constraint.
Detection asks "is this state at least `k`-entangled"; we ask "which matrices are reachable
by a width-`k` program". Different quantifier.

**Verified relationship.** The known bound is exactly the support function of `Q^prog_{m,k}`
in the single direction `b = 1`:

```
max { 1ᵀ F 1 : F ∈ Q^prog_{m,k} }  =  ⌊m/k⌋k² + r²
```

Checked by exhaustive extreme-point enumeration for every `(m, k)` with `m ≤ 7`, `k ≤ 3`:
exact agreement in all 15 cases. At `k = 2` the matching polytope reproduces it in closed
form as `m + 2⌊m/2⌋`, agreeing for `m` up to 12 without any enumeration.
(`prototype/tests/test_covq.py::test_support_function_at_b_ones_reproduces_k_producibility_bound`.)

**Consequence for positioning.** Say: *the standard `k`-producibility bound is the value of
our polytope's support function in one direction; we compute it in every direction, and at
`k = 2` we give the complete facet description.* That is defensible, it credits the prior
work properly, and it doubles as an independent physics check on the whole construction.

**Risk: LOW.**

## A4 — Matchings in quantum gate scheduling

**Status: same combinatorial object, different role. A presentation hazard, not a collision.**

Matchings are the standard model for parallel two-qubit gate scheduling: two gates run
concurrently iff they act on disjoint qubits, so assigning gates to time slices is edge
colouring. This is well established in the compilation literature.

The role differs:

| | gate scheduling | CovQ |
|---|---|---|
| what a matching is | a **time slice inside one circuit** | a **branch of a randomised schedule** |
| what is being decomposed | a fixed circuit's gate list | a target **QFIM** |
| objective | depth / makespan | feasibility, then resource cost |
| the convex combination | none — it is a partition of a gate set | the design measure `p_r` |

**Risk: MEDIUM as presentation.** A referee who reads "matching decomposition" in a
compilation paper will assume edge colouring. Pre-empt it in the introduction, in one
sentence, before the theorem.

## A5 — Target QFIM → circuit synthesis

**Status: no direct prior art located, on a second independent search.**

Searches for target-QFIM synthesis, observable-covariance circuit compilation, and inverse
QFIM/state-synthesis returned nothing that accepts a QFIM and emits circuits. The literature
splits into QFIM *evaluation* (and the observation that it is typically ill-conditioned, so
inverting it introduces error — which supports our C6/conditioning reformulation) and generic
circuit synthesis from unitaries or states.

**Consequence.** The pipeline framing survives; K1b remains not established. This is
consistent with the v0.2 audit and does not change it.

---

## Implications for the gates

```
K2  no theorem beyond a known polytope translation   NOT FIRED (conditional on framing)
K6  width-2/width-3 boundary false or already known  LIVE — and now more likely "already
                                                     known" than "false"
K1b end-to-end compiler collision                    still NOT ESTABLISHED
```

**On K6.** Partition into Triangles is a classical NP-complete problem, and partitioning into
cliques of size ≤ 3 is hard by the same token. The reduction sketched in `DECISION_001` is
therefore likely to go through — and to go through in a paragraph. Expect a **corollary**,
not a theorem. Two consequences:

1. Step 3 of Task 0B.5 should be re-scoped from "prove or abandon" to "prove, then determine
   whether it is already stated in the bounded-cluster partitioning literature". The failure
   mode is no longer that it is false; it is that it is folklore.
2. The paper's complexity content will lean on classical results at **both** ends: Edmonds at
   `k = 2`, Partition into Triangles at `k ≥ 3`. The quantum content is the identification and
   the constructive circuit emission, not the combinatorics. This should be stated to the
   reader plainly rather than discovered by them.

That is a real reduction in the paper's ambition and the charter should absorb it. It is not
a kill: an exact, constructive, hardware-native tractability boundary for QFIM compilation is
still worth publishing in a compilation venue. But it is not a complexity-theory contribution
and must not be sold as one.

## What step 2 should now prove

Only the parts that are actually ours:

1. **Necessity** — a width-2 branch has off-diagonal support on a matching, so `(|F_ij|)` is
   dominated by a convex combination of matching incidence vectors, and the matching polytope
   is down-closed. *(Ours.)*
2. **The magnitude reduction** — signs impose no constraint beyond magnitudes. *(Ours; the
   step that makes it an iff.)*
3. **Sufficiency** — cite Edmonds for the decomposition, then give the circuit for each
   matching. *(Construction ours, decomposition inherited.)*
4. **Hardware-native version** on `MATCH(H)`. *(Ours; the compiler-relevant form.)*
5. **Support-function corollary** recovering `⌊m/k⌋k² + r²` at `b = 1`. *(Ties it to the
   physics literature; verified numerically.)*

## Supplementary sweep — the information-floor interface

Run after the decision to promote the information floor to the main problem, to check the
interface itself rather than the `k = 2` theorem. Five terms: information-matrix lower bound
in optimal design, Loewner-constrained design, QFIM domination, minimum-resource Fisher target,
semidefinite information-covering design.

**Finding: the interface has a clear classical shadow. Do not claim it in isolation.**

- **Loewner-monotone design is standard.** Common optimality criteria are monotone with respect
  to the Loewner ordering, and spectral design with a prior information matrix updated by
  rank-one design vectors under a norm budget is an established formulation.
- **Fisher information as a *constraint* rather than an objective is established.** Sparse
  sensor selection routinely imposes performance as LMI constraints — Fisher-information lower
  bounds in localisation — and minimises cost or sensor count subject to them, via SDP
  relaxations, reweighted ℓ1, and greedy methods. That is structurally the same interface as
  `min cost s.t. AᵀFA ⪰ G_req`.
- **On the quantum side**, QFIM equals the minimum covariance over purifications, and the
  matrix Cramér–Rao statement `Σ ⪰ F⁻¹/M` is standard. Nothing found compiles *to* a
  requirement.

**Consequence.** The earlier phrasing "no prior-art shadow at all" was wrong and is retracted.
The defensible claim is the composite object only:

> We found no prior work that combines an information-floor compiler interface with an exact
> bounded-entanglement QFIM achievable set and constructive hardware-native circuit emission.

This also motivates the new kill gate **K9**: if every benchmark contract is satisfied by the
product probe or by one fixed state, the composite claim is empty and the quantum feasible set
is doing no work. The three analytic benchmarks in the charter exist to make K9 answerable —
the common-mode curve, the odd/even ceiling, and the overlapping-mode trade-off all have
closed-form separations from the product probe.

**Supplementary sources.** [Optimal spectral design with prior
information](https://arxiv.org/html/2605.27837v1) ·
[minimal upper bounds in the Loewner order](https://ar5iv.labs.arxiv.org/html/2606.18173) ·
[spatially constrained sensor placement (LMI Fisher constraints)](https://www.emergentmind.com/topics/spatially-constrained-sensor-placement) ·
[Fisher-information-based sensor placement, analytic results and benchmarks](https://arxiv.org/html/2602.02981) ·
[optimal sensor placement using modified Fisher information matrices](https://journals.sagepub.com/doi/full/10.1177/15501477211023022) ·
[fundamental and operational limitations of Fisher-information-based quantum metrology](https://arxiv.org/html/2603.08306)

## Sources

Matching polytopes and decomposition: [Edmonds' matching polytope theorem
(EGRES)](http://lemon.cs.elte.hu/egres/open/Edmonds'_matching_polytope_theorem) ·
[EPFL DISOPT notes](https://www.epfl.ch/labs/disopt/wp-content/uploads/2018/09/matching_polytope.pdf) ·
[Schrijver, polyhedral combinatorics](https://homepages.cwi.nl/~lex/files/paris_for_smf.pdf) ·
[symmetric Birkhoff–von Neumann decomposition algorithms](https://inria.hal.science/hal-04877502v2/document)

Clique partitioning: [Grötschel–Wakabayashi, facets](https://www.semanticscholar.org/paper/Facets-of-the-clique-partitioning-polytope-Gr%C3%B6tschel-Wakabayashi/2e429debdfb9e63c2820dff890433865674c5fa5) ·
[Oosten et al., facets and patching facets](https://onlinelibrary.wiley.com/doi/abs/10.1002/net.10004) ·
[chorded cycle facets](https://arxiv.org/pdf/2411.03407) ·
[cluster deletion and clique partitioning with bounded clique number](https://arxiv.org/pdf/2505.00922) ·
[redundancy of transitivity constraints](https://arxiv.org/pdf/2605.01481)

QFIM and entanglement depth: [Pezzè–Smerzi, quantum metrology from a quantum information
science perspective](https://arxiv.org/pdf/1405.4878) ·
[uncertainty relations between QFI and entanglement monotones (arXiv:2501.14595, retitled)](https://arxiv.org/html/2501.14595v2) ·
[QFI and multipartite entanglement in spin chains](https://arxiv.org/pdf/2307.02407)

Gate scheduling as edge colouring: [two-step approach to scheduling quantum
circuits](https://arxiv.org/pdf/1708.00023) ·
[compilation for dynamically field-programmable qubit arrays](https://arxiv.org/pdf/2405.15095)

Cut and correlation polytopes: [lecture notes on the cut polytope and SOS
cones](https://homes.cs.washington.edu/~jrl/teaching/cse599Isp21/notes/lecture12.pdf)
