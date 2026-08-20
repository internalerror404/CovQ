# Task 0A/0B — Novelty and theorem audit, v0.2

**Scope:** the four v0.1 documents in `archive/v0.1/`.
**Outcome:** pivot. See `docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md`.
**Evidence:** `prototype/` (see §9). Every number quoted below is reproduced by it.

v0.2 supersedes the first-pass review. Three of its conclusions were wrong or overstated
and are corrected here in §7; a fourth correction, §5.3, is the prototype's own and
retracts a measurement the first pass would have carried into the charter.

---

## 1. Four structural collapses in the declared target sector

With `U_θ = exp[−(i/2) Σ_i θ_i P_i]`, `P_i² = I` commuting, and a pure probe,

```
F_ij = Re⟨P_i P_j⟩ − μ_i μ_j ,      μ_i = ⟨P_i⟩ ,
```

so `F_ii = 1 − μ_i²`.

**L1 — the sector is one condition written twice.** `F_ii = 1 ⟺ μ_i = 0`. "Zero-mean *and*
unit-diagonal" over-specifies.

**L2 — the zero-mean hypothesis in the labelled normal form is never an assumption.** Every
branch obeys `F_{r,ii} ≤ 1`. If `F_⋆ = Σ_r p_r F_r` has `F_{⋆,ii} = 1`, then every
positive-weight branch has `F_{r,ii} = 1`, hence `μ_{r,i} = 0`. Theorem B's stated
hypothesis is implied by its conclusion's setting.

**L3 — gate C2 tests a definition.** For a labelled backend `F_Π` *is defined* as the branch
average, so `‖F_Π − Σ_r p_r F_r‖ < 10⁻¹²` cannot fail. Re-specified in §6.

**L4 — Loewner dominance degenerates to exact matching.** If `F ⪰ F_⋆` and
`diag F = diag F_⋆ = 1`, then `D = F − F_⋆ ⪰ 0` with `diag D = 0`; positive semidefiniteness
gives `|D_ij|² ≤ D_ii D_jj = 0`, so `D = 0`. Inside the declared sector, exact matching,
approximate matching, and one-sided information guarantees are not three semantics. They
are one.

The consequence is not cosmetic: the sector chosen for the "initial exact theory" is
precisely the sector in which the distinctive machinery has nothing to distinguish.

---

## 2. The coherent-flag formula, and the loophole it opens

For an orthogonally flagged pure state `|Ψ⟩ = Σ_r √p_r |r⟩_f |ψ_r⟩_d` with generators acting
only on data,

```
F_Ψ  =  Σ_r p_r F_r  +  Cov_p(μ_r) .
```

Equality with the labelled schedule holds **iff the branch-mean vectors coincide** on the
positive-weight support — `μ_r = μ̄` for all `r` with `p_r > 0`. Zero means are sufficient,
not necessary. v0.1's Theorem B states the sufficient condition as if it were the criterion.

"Retain the flag through readout" is ambiguous and the ambiguity is load-bearing. Measuring
or dephasing the flag in `{|r⟩}` kills the covariance term and returns the classical–quantum
labelled program. Preserving coherence and permitting a joint flag–data POVM does not.

### The loophole

Fix any feasible `F ∈ conv{ss^T}` and a zero-mean sign distribution `p` with
`E_p[s] = 0`, `E_p[ss^T] = F`. Prepare

```
|Ψ_p⟩ = Σ_s √p(s) |s⟩_f |s⟩_d .
```

Every conditional data branch `|s⟩_d` is a product eigenstate: `F_s = 0`, conditional data
width 1. Yet `F_{Ψ_p} = Cov_p(s) = F`. Dephasing the flag sends the QFIM to exactly zero,
because `ρ_data` is diagonal in the Z basis and therefore commutes with every generator.

> **Conditional data-branch width is not a valid resource measure for a general
> coherent-flag backend.**

Otherwise every feasible target has coherent "width one". Nothing became cheap; the resource
moved into flag–data entanglement, flag dimension, Schmidt rank across the flag|data cut,
controlled-preparation cost, and joint readout. A single branch-width number shared across
labelled and coherent implementations is unsound.

**Adopted treatment:** labelled programs are the sole subject of the bounded-width theorem;
coherent flags get a semantics section with their own resource vector.

---

## 3. Novelty: what collides and what does not

### 3.1 QUEST — K1 splits

QUEST (arXiv:2605.02367) maps `{O_i, τ_i} ↦ |ψ⟩` with `⟨O_i⟩ = τ_i` under a feasibility
promise, and includes a worked example matching both `Z_i` and `Z_iZ_j` moments. By L1, a
unit-diagonal target QFIM *is* `C(m,2) + m` linear constraints on commuting Pauli
expectations. So the v1 target is a special case of QUEST's input language, on an especially
easy observable family.

**Dead as a positioning claim:** "traditional state preparation specifies a state; CovQ
specifies only a weaker information geometry." A QFIM is mathematically weaker than a
wavefunction, but partial observable specification is not new.

**Not dead as a compiler claim.** QUEST assumes feasibility rather than certifying it, is
not an exact cut-polytope membership algorithm, formulates no bounded branch width, gives no
matching-polytope characterisation, does not optimise a labelled schedule over a hardware
graph, and emits no infeasibility certificates.

```
K1a  specification / abstract collision      : FIRED
K1b  identical end-to-end compiler collision : NOT ESTABLISHED
```

The first-pass review reported K1 as fired outright. That overstates it. The abstract and
introduction must be rewritten; the compiler need not be abandoned.

### 3.2 Correlation-polytope theory — A and E demote to background

Huber and Marić (arXiv:1706.06182) already characterise symmetric Bernoulli correlations as
a polytope affinely isomorphic to a cut polytope, with its vertices. The quantum step
`p(s) ↦ Σ_s √p(s)|s⟩` is immediate, because commuting-Z QFI depends only on
computational-basis probabilities. Theorem A therefore reads:

> **Background Proposition 1.** By symmetric Bernoulli correlation-polytope theory together
> with the pure-state amplitude realisation, the feasible zero-mean unit-diagonal QFIMs are
> exactly `conv{ss^T : s ∈ {±1}^m}`.

Caprara et al. (arXiv:2605.02896) establish hardness for membership and decomposition rank
over correlation polytopes and cones, noting membership hardness was already known. So:
full feasibility membership — prior art; minimum exact setting count — hardness prior art;
generic sparse convex approximation — prior art. Exact algorithms at small `m` remain
legitimate implementation contributions, and target-specific tractable subclasses remain
open and valuable.

### 3.3 QFIM k-producibility

`Q^prog_{m,k}` is likely already present in the QFIM entanglement-depth literature as an
*outer bound on states*. The defensible distinction is narrow but real: QFI convexity puts
`k`-producible mixed states *below* the hull, whereas a labelled program *attains* it. Frame
as "we convert a known witness relaxation into an achievability statement, with emitted
circuits". Do not frame as a new resource hierarchy.

### 3.4 Optimal experimental design — embrace it

`F_Π = Σ_r p_r F_r` is a discrete optimal design: branch circuits are design points, branch
probabilities are design weights, and the protocol's Frank–Wolfe arm is Fedorov–Wynn. There
is no citation to Kiefer–Wolfowitz, Fedorov, Pukelsheim, or Boyd–Vandenberghe §7.5 anywhere
in v0.1. The correct related-work sentence is that CovQ supplies the quantum design space,
the exact bounded-entanglement feasible sets, and emitted hardware programs to which
classical optimal-design machinery then correctly applies.

### 3.5 Bibliography correction

arXiv:2501.14595 has been retitled to *Uncertainty relations between quantum Fisher
information and entanglement monotones*. Cite the current title and version, and confirm
whether the dimensionality result relied on survived revision.

---

## 4. The surviving theorem

**Theorem (width two).** For the labelled, zero-mean, physical-Z backend,

```
F ∈ Q^lab_{m,2}   ⟺   (|F_ij|)_{i<j} ∈ MATCH(K_m) ,
```

equivalently `Σ_{j≠i} |F_ij| ≤ 1` for every `i`, and
`Σ_{{i,j}⊆S} |F_ij| ≤ (|S|−1)/2` for every odd `S ⊆ [m]`.

*Necessity.* A width-2 branch has off-diagonal support on a matching with entries in
`[−1,1]`, so `|F_ij| ≤ Σ_r p_r |c^r_ij| 1[ij ∈ M_r]`, dominated componentwise by
`Σ_r p_r χ_{M_r} ∈ MATCH(K_m)`; the matching polytope is down-closed.

*Sufficiency, constructively.* Edmonds writes `x = Σ_M λ_M χ_M`. Realise matching `M` by a
two-qubit signed cat on each edge carrying `sgn(F_ij)`, and `|+⟩` on every unmatched qubit.
Then `Σ_M λ_M 1[ij ∈ M] sgn(F_ij) = F_ij`.

**Hardware-native form — the one the compiler needs.** For `H = (V, E)` with no routing
inside the primitive backend,

```
Q^lab_{H,2} = { F : diag F = 1,  F_ij = 0 for ij ∉ E,  (|F_ij|)_{ij∈E} ∈ MATCH(H) } .
```

This is a direct map from a target QFIM to a schedule of *parallel Bell pairs on physical
edges*. Routing enters as a distinct backend that enlarges the admissible edge set at an
explicit SWAP/path cost, priced rather than assumed.

**What it delivers:** an exact resource characterisation; polynomial-time feasibility via
Edmonds separation, against NP-hard membership at `k = m`; explicit separating certificates;
constructive circuit schedules; a sharp product-versus-pair-entanglement boundary; and at
most `O(m²)` branches by Carathéodory.

**What it does not deliver.** It clears K4. It does not by itself carry a major
quantum-computing paper, because the combinatorics are classical matching theory. The paper
needs a complexity transition — width 2 reduces to maximum-weight matching, width 3
conjecturally NP-hard, candidate reduction Partition into Triangles with utility 1 on edges
and 0 on non-edges, so that utility `3q` on `3q` vertices is achievable iff the vertices
partition into `q` triangles. Unproved. Requires comparison against the size-constrained
clique-partitioning polytope literature.

---

## 5. Baselines

### 5.1 The sparse pure state is mandatory

A labelled schedule is never needed merely to prove feasibility. Applying Carathéodory to
the augmented first-and-second-moment vectors `(s, offdiag(ss^T))` bounds the support by

```
m + C(m,2) + 1  =  m(m+1)/2 + 1 ,
```

which is **79** at `m = 12`. (The first-pass review's 134 came from
`2(C(m,2)+1)` on the symmetric construction: valid, but not the tightest elementary bound.)
The resulting single pure state uses one setting.

Gate counts depend on the sparse-preparation algorithm and topology, so "`O(m³)` gates"
should be stated as a generic polynomial upper bound, never as measured or optimal.

v0.1 lists only `generic_full_state_preparation_then_transpile`. That guarantees CovQ an
artificial advantage. It must be supplemented by four realisation-matched baselines built
from the *same* decomposition:

1. sparse Carathéodory pure state;
2. randomised cat schedule;
3. one coherent sparse data state;
4. an explicit flag–data purification.

This isolates the real trade: **settings ↔ per-shot gates ↔ ancillas/coherent control ↔
noise.**

### 5.2 Measured (prototype, seed 2026, all-to-all, error `< 3e-16` throughout)

| m | atoms | program | settings | CX | 2q depth | anc |
|---|---|---|---|---|---|---|
| 4 | 4 | labelled global cat | 4 | 12 | 3 | 0 |
| 4 | 4 | Carathéodory sparse state | 1 | 184 | 168 | 1 |
| 4 | 4 | flagged cat (literal) | 1 | 334 | 310 | 3 |
| 4 | 4 | dense state prep *(weak baseline)* | 1 | 256 | 232 | 1 |
| 5 | 3 | labelled global cat | 3 | 12 | 4 | 0 |
| 5 | 3 | Carathéodory sparse state | 1 | 328 | 292 | 2 |
| 5 | 3 | flagged cat (literal) | 1 | 314 | 290 | 3 |
| 5 | 3 | dense state prep *(weak baseline)* | 1 | 1000 | 876 | 2 |
| 6 | 4 | labelled global cat | 4 | 20 | 5 | 0 |
| 6 | 4 | Carathéodory sparse state | 1 | 712 | 620 | 3 |
| 6 | 4 | flagged cat (literal) | 1 | 490 | 450 | 3 |
| 6 | 4 | dense state prep *(weak baseline)* | 1 | 3472 | 2996 | 3 |

Sparsity is worth 3–5× against the dense baseline at `m = 5, 6`. **The apparent
labelled-schedule win is not admissible as stated:** the prototype's sparse preparation is a
generic prefix-tree construction with v-chain lowering, not a state-of-the-art routine, so
every single-state number above is an upper bound on that baseline's true cost. Until a
competitive sparse-preparation routine is implemented, no per-shot gate-count advantage over
the single-state realisation may be claimed.

### 5.3 Correction — a retracted measurement

The first pass reported that a "relaxed" coherent-flag construction realises the same QFIM
roughly 12× cheaper than the charter's literal controlled-cat normal form. **That comparison
is invalid.** The cheap construction buys its gate count by moving the information into
flag–data coherence — it is the §2 loophole. Dephasing its flag sends its QFIM to zero,
while the literal form is unaffected:

| program | cond. width | flag qubits | Schmidt | CX | tr F coh | tr F deph | coherence load-bearing |
|---|---|---|---|---|---|---|---|
| `sign_purification` | 1 | 4 | 8 | 188 | 4.000 | 0.000 | **yes** |
| `flagged_cat_relaxed` | 1 | 3 | 8 | 50 | 4.000 | 0.000 | **yes** |
| `flagged_cat_literal` | 4 | 2 | 4 | 334 | 4.000 | 4.000 | no |

Programs on opposite sides of the dephasing test are not comparable on gate count. A cheaper
coherent program that fails dephasing has not compiled the same object.

---

## 6. Gates

Three v0.1 gates cannot fail. Re-specified:

**C2 → C2a / C2b / C2c.**

```
C2a  F_CQ            = Σ_r p_r F_r
C2b  F_coh           = Σ_r p_r F_r + Cov_p(μ_r)      ← must include non-zero-mean branches
C2c  F_dephased flag = Σ_r p_r F_r
```

C2b's non-zero-mean test branches may sit outside the main benchmark sector; without them
the identity is untested. C2c is additionally the **discriminator** that decides whether a
given program is labelled-equivalent (dephasing costs nothing, width is a resource) or
genuinely coherent (dephasing destroys the QFIM, width is meaningless and flag dimension,
Schmidt rank and joint readout must be priced instead).

**C3.** "Exact solver agrees with exhaustive enumeration" is hollow when the exact solver
*is* the enumeration. Re-specified as agreement between three independent code paths —
full-vertex LP, column generation with exact pricing, hypermetric certificate search — on a
suite containing feasible, infeasible **and** boundary instances.

**C4 — instantiate it.** v0.1 gives no ground-truth instance. Use

```
F_out = [[ 1.0, −0.4, −0.4],
         [−0.4,  1.0, −0.4],
         [−0.4, −0.4,  1.0]]
```

Spectrum `{0.2, 1.4, 1.4}`: positive definite, unit diagonal. But
`F₁₂ + F₁₃ + F₂₃ = −1.2 < −1` violates a correlation-polytope facet. Triangle inequalities
stop sufficing past `m = 3`, so `m ≥ 5` instances need hypermetric ones.

The certificate family has a genuinely quantum reading worth using in the paper: for integer
`b` with `Σ_i b_i` odd, every realisable `F` obeys

```
b^T F b = Var(Σ_i b_i P_i) ≥ 1 ,
```

because `Σ_i b_i P_i` has integer spectrum of the parity of `Σ_i b_i`. A violated inequality
is an **uncertainty-relation obstruction**, not merely a separating hyperplane. The instance
above is `b = (1,1,1)`, `b^T F b = 0.6`, margin `0.4`.

**C10 → conditioning-aware stability.** `‖AᵀΔF A‖ ≤ ‖A‖²‖ΔF‖` is submultiplicativity: it
holds for every matrix, so it is not a test. The operational object is `(AᵀFA)^{-1}`, or a
pseudoinverse on the identifiable quotient. With `B_⋆ = AᵀF_⋆A ⪰ γI` on its declared support
and `‖Aᵀ(F − F_⋆)A‖ ≤ ε < γ`,

```
‖B^{-1} − B_⋆^{-1}‖ ≤ ε / (γ(γ − ε)) .
```

If the kernels differ, **no** stable inverse guarantee follows. Report `rank F`, `λ⁺_min(F)`,
`κ⁺(F)`, kernel principal angles, and `‖F⁺ − F_⋆⁺‖`. This also resolves the awkwardness that
the polytope's extreme points `ss^T` are rank one: feasible extremes, usually terrible
full-field estimation targets.

**C6 — split the width metric.** v0.1 runs physical-Z and Clifford-conjugated generators
under one width number. Report `k_logical`, `k_physical`, and the Clifford two-qubit cost
separately. Minimum-physical-entanglement language is valid for the physical-Z arm only. The
prototype measures the gap in **both** directions: a Z-frame width-3 branch that is
physically a *product* state (so the frame is a compiler degree of freedom, not just a
caveat), and Z-frame width-`m` branches that stay width-`m`.

**C11 — the missing lemma.** Nothing in v0.1 establishes that the multiparameter bound is
attainable. It is:

```
Im ⟨∂_i ψ | ∂_j ψ⟩ = −(1/4) Im ⟨ψ| P_i P_j |ψ⟩ = 0
```

identically, for commuting Hermitian generators. State it with its qualifications:
attainability is *local*, at the operating point; the optimal POVM may depend on that point;
singular targets must be restricted to their identifiable support; labelled programs permit
branch-conditioned measurements while coherent programs require a joint flag–data
measurement; and hardware runs must report an explicit attainable *classical* Fisher matrix
unless the measurement is actually compiled. Without this, "certified downstream performance"
is unsupported.

**Enumeration caps.** Split `vertex_generation_max_m`, `exact_membership_max_m`, and
`minimum_support_MILP_max_m`. At `m = 12` the quotiented sign-vertex set has only `2¹¹ =
2048` columns against 66 rows: vertex generation is not the limiting operation. Minimum-
support MILP and certificate generation may still be.

---

## 7. Corrections to the first-pass review

1. **The relaxed ℓ1-rank hardness result is a cone result.** For a convex decomposition,
   `p ≥ 0` and `Σ p_r = 1` give `‖p‖₁ = 1` identically; Caprara et al. state explicitly that
   relaxed rank is therefore meaningless on the convex hull. It was also attributed to a
   protocol arm that does not exist: v0.1 has Frank–Wolfe sparse approximation, not ℓ1-rank
   minimisation. Frank–Wolfe support growth and exact ℓ1-rank minimisation are different
   objects and must not be conflated.
2. **K1 was reported fired; it is half-fired.** See §3.1.
3. **`F_ii ≤ 1` was proposed as the fix for §1; it is not safe.** Sub-unit diagonals do
   unfreeze the means, `|μ_i| = √(1 − F_ii)`, but feasibility is then a *joint* moment
   condition: with `M = F + μμᵀ`, one needs
   `(μ, M) ∈ conv{(s, ss^T) : s ∈ {±1}^m}`. Checking `M ∈ Q_m` is not enough, because the
   realising distribution must also reproduce the prescribed first moment. With only `F`
   given, this becomes an existential nonlinear projection over mean signs and joint
   moments, and the Bernoulli-compatibility / prescribed-moment literature has not been
   audited. "Genuinely uncharacterised" was not a safe claim.
4. **The 12× coherent-flag speedup is retracted.** §5.3.

### Adopted instead of (3): the information floor

Keep the realised probe sector unit-diagonal, `F_Π ∈ Q^lab_{m,k}`, and change the *target*
semantics:

```
min_Π  Cost(Π)      subject to      F_Π ⪰ G_⋆ ,
```

with `G_⋆ ⪰ 0` and possibly subunit-diagonal. Loewner dominance is now nontrivial, because
`G_⋆` need not have unit diagonal, and §1's L4 collapse does not apply. It also answers the
motivation question a referee asks first. Nobody prescribes an arbitrary covariance for its
own sake; a user declares the *minimum information downstream estimation requires*, and the
compiler returns the cheapest physical QFIM dominating it. For the width-2 hardware-native
backend this is a matching-polytope-constrained semidefinite problem with explicit
certificates. Exact matching survives as a verification mode and benchmark family.

---

## 8. Verdict table

| v0.1 item | Verdict | Proper role |
|---|---|---|
| Theorem A: full feasibility | correct but occupied | background proposition |
| Signed-cat primitive | correct but elementary | compiler primitive |
| Labelled additivity | correct | semantics lemma / unit test |
| Coherent-flag equivalence | correct only under an extra restriction | reformulate; separate resource model |
| `k = 2` bounded width | strong candidate | **main theorem** |
| Full-width hardness | occupied | complexity context |
| Downstream norm bound | elementary | background; replace with inverse stability |
| Task 0C v0.1 | premature | supersede before coding |

---

## 9. Numerical evidence

From `prototype/`, reproducible with the protocol seeds.

**The `k = 2` theorem, against brute force.** Edmonds separation versus full enumeration of
the extreme points of `Q^prog_{m,2}`, on 160 random unit-diagonal symmetric matrices
(`m ∈ {3,4,5,6}`, deliberately including infeasible ones): **0 disagreements**, and every
accepted instance reconstructs to `< 1e-8`.

**The odd-set facets are active, not redundant.** For a triangle carrying common
off-diagonal `c`, the degree constraints alone permit `c ≤ 1/2`, while the blossom
constraint permits only `c ≤ 1/3`. Both the Edmonds test and brute force switch at exactly
`1/3`, and the whole family stays inside `Q_m`:

| `c` | in `Q_m` | width-2 (Edmonds) | width-2 (brute force) | violated facet |
|---|---|---|---|---|
| 0.300 | yes | yes | yes | — |
| 0.3333 | yes | yes | yes | — (tight) |
| 0.340 | yes | no | no | odd set `{0,1,2}`, margin 0.02 |
| 0.400 | yes | no | no | odd set `{0,1,2}`, margin 0.20 |
| 0.450 | yes | no | no | odd set `{0,1,2}`, margin 0.35 |

So there is a family of perfectly feasible QFIMs that provably require width ≥ 3 for a
blossom reason. That is the sharpest available demonstration that the characterisation has
content beyond degree counting.

**The path threshold.** `F_{i,i+1} = c` has interior degree sums `2c`, so width-2 holds iff
`|c| ≤ 1/2`, while positive semidefiniteness survives to `1/(2cos(π/(m+1)))`. Reproduced
exactly at `m ∈ {4,5,6}`: feasible at `c ∈ {0.45, 0.49, 0.50}`, infeasible at
`c ∈ {0.51, 0.55}`. The window between the two thresholds is where a target is a perfectly
good QFIM that provably needs a width-3 branch.

**Hardware-native membership discriminates.** The `c = 0.4` path target is width-2 feasible
on `line` and `all_to_all` (3 branches) and infeasible on `square_grid`, rejected with an
`off_graph_support` certificate naming edge `(2,3)`.

**Convexity gap.** For a signed-cat schedule the unlabelled mixed state attains the branch
average exactly (trace gap `< 3e-15`), because `Z_i|C_s⟩ = s_i|C⁻_s⟩` maps the support into
the kernel and `(q_k−q_l)²/(q_k+q_l) = q_k+q_l` iff `q_k q_l = 0`. So the "forbidden
shortcut" costs nothing *for the compiler's own output family* — a finding about the
program-semantics section, not a licence. The positive control, `(|C_s⟩, |C⁻_s⟩)` at equal
weight, has covariance `ss^T` and QFI exactly zero: relative gap `1.0`.

**Gates.** `C1 C2 C2c C3 C4 C5 C6 C7 C8 C9 C11` PASS, `C10 C12` MEASURED, 0 FAIL, at
`m ∈ {3,4,5}`.

---

## 10. Retained from v0.1 without change

Novelty-before-code ordering. The `forbidden_shortcut` declaration. Resources derived from
emitted circuits. Pilot-then-freeze thresholds with no directional success criteria set in
advance. Spacetime as one deletable late benchmark. The stop rule: write the report, abort,
change no tolerance.
