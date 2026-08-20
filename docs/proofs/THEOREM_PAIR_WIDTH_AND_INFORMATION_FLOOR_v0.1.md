# Pair-width geometry and the information-floor compiler

Task 0B.5 step 2, parts 2A–2F. Every numerical claim is checked in
`prototype/tests/test_covq.py` (85 tests).

Throughout: `U_θ = exp[−(i/2) Σ_i θ_i P_i]` with `P_i = Z_i`; per-shot QFIM
`F_ij = ⟨Z_iZ_j⟩ − ⟨Z_i⟩⟨Z_j⟩`; realised probes zero-mean, so `diag F = 1`.
`Q^lab_{H,2}` is the set of per-shot QFIMs of labelled schedules whose every branch is a
pure product state over blocks of size ≤ 2, each 2-block sitting on an edge of `H = (V,E)`.

---

## 2A — Necessity

**Lemma (matching support).** Every width-2 branch `F_r` has `diag F_r = 1`, off-diagonal
support contained in a matching `M_r ⊆ E`, and entries in `[−1,1]` there.

*Proof.* `ψ_r = ⊗_B ψ_B` with `|B| ≤ 2`. For `i, j` in different blocks the observables act
on independent factors, so `Cov(Z_i,Z_j) = 0`. The blocks of size 2 are vertex-disjoint by
definition of a partition and lie on edges of `H` by hypothesis, so they form a matching.
`F_{r,ii} = 1 − ⟨Z_i⟩² = 1` in the zero-mean sector, and Cauchy–Schwarz then gives
`|F_{r,ij}| ≤ √(F_{r,ii} F_{r,jj}) = 1`. ∎

**Lemma (down-closedness).** If `x ∈ MATCH(H)` and `0 ≤ z ≤ x` componentwise, then
`z ∈ MATCH(H)`.

*Proof.* Edmonds' description consists of `x ≥ 0`, `Σ_{e∋i} x_e ≤ 1`, and
`Σ_{e⊆S} x_e ≤ (|S|−1)/2` for odd `S`. Every one of these has non-negative coefficients on
the left and a constant on the right, so it is preserved under componentwise decrease. ∎

**Proposition (necessity).** `F ∈ Q^lab_{H,2}` implies `diag F = 1`, `F_ij = 0` for
`ij ∉ E`, and `(|F_ij|)_{ij∈E} ∈ MATCH(H)`.

*Proof.* `F = Σ_r p_r F_r`. Off-graph entries vanish in every branch, hence in the average.
By the triangle inequality,

```
|F_ij|  ≤  Σ_r p_r |F_{r,ij}| · 1[ij ∈ M_r]  ≤  Σ_r p_r χ_{M_r}(ij) ,
```

and `Σ_r p_r χ_{M_r}` is a convex combination of matching incidence vectors, so it lies in
`MATCH(H)`. Down-closedness carries `(|F_ij|)` into `MATCH(H)`. ∎

## 2B — Magnitude sufficiency

**Proposition.** If `diag F = 1`, `F_ij = 0` off `E`, and `x := (|F_ij|)_{ij∈E} ∈ MATCH(H)`,
then `F ∈ Q^lab_{H,2}`, and a schedule realising it can be written down.

*Proof.* By Edmonds, `x = Σ_M λ_M χ_M` with `λ ≥ 0`, `Σ_M λ_M = 1`, over matchings of `H`
(the empty matching absorbs any slack). Fix a matching `M`. Because `M` is a matching, each
vertex lies in at most one edge, so the assignment

```
s_i = +1 for the lower-indexed endpoint of its edge,
s_j = sgn(F_ij) for the higher-indexed endpoint,
s_i = +1 for unmatched i
```

is consistent — **no vertex receives two conflicting demands.** Realise branch `M` by a
two-qubit signed cat `(|s_is_j⟩ + |−s_i,−s_j⟩)/√2` on each edge of `M` and `|+⟩` on every
unmatched qubit. Its QFIM is `1` on the diagonal, `sgn(F_ij)` on `ij ∈ M`, and `0` elsewhere.
Averaging,

```
Σ_M λ_M (branch_M)_ij  =  sgn(F_ij) Σ_M λ_M χ_M(ij)  =  sgn(F_ij) · x_ij  =  F_ij ,
```

and the diagonal averages to `1`. ∎

**Remark (why the signs are free, and why this is special to `k = 2`).** The sufficiency
proof works because a matching's edges are vertex-disjoint, so edge signs may be chosen
independently. At width 3 this fails: inside a 3-block the three signs obey
`(s_is_j)(s_js_k)(s_is_k) = s_i²s_j²s_k² = +1`, a parity constraint linking the three
off-diagonals. Magnitudes therefore *cannot* determine feasibility at `k ≥ 3`, and the clean
absolute-value characterisation is a genuine feature of the pair level rather than the first
case of a pattern.

## 2C — Hardware-graph form

Combining 2A and 2B:

> **Theorem 1.** `F ∈ Q^lab_{H,2}` **iff** `diag F = 1`, `F_ij = 0` for `ij ∉ E`, and
> `(|F_ij|)_{ij∈E} ∈ MATCH(H)`.

Taking `H = K_m` recovers the complete-graph form
`Σ_{j≠i}|F_ij| ≤ 1` and `Σ_{{i,j}⊆S}|F_ij| ≤ (|S|−1)/2` for odd `S`.

Membership is decidable in polynomial time by Edmonds separation, in contrast to NP-hard
membership in `Q_m` at unbounded width. *Verified against brute-force enumeration of the
extreme points of `Q^prog_{m,2}` on 160 random instances at `m ∈ {3,4,5,6}`: zero
disagreements.*

## 2D — Constructive emission

Each branch is a disjoint union of two-qubit cats, so per branch:

- two-qubit gate count `= |M| ≤ ⌊m/2⌋`;
- **two-qubit depth exactly 1** — every edge of a matching executes in parallel by
  construction, which is the whole point of compiling to a matching;
- no ancillas, no routing (every edge is native).

**Support bound.** The decomposition lives in the `|E|`-dimensional edge space with one
normalisation, so Carathéodory gives at most `|E| + 1` branches. Column generation typically
returns far fewer.

**Cost identity.** For non-negative native edge costs `c_e`, the expected per-shot bill of the
schedule is

```
Σ_r p_r Σ_{e ∈ M_r} c_e  =  Σ_e c_e Σ_r p_r χ_{M_r}(e)  =  Σ_e c_e |F_e| .
```

The objective in 2E is therefore not a proxy for cost. It *is* the expected native Bell-pair
cost per shot. *Verified: residual < 1e-9 on every emitted program.*

## 2E — The information-floor compiler

The compiler's contract is a downstream one. Given `A ∈ R^{m×d}` naming the parameter
combinations that matter and `G_req ⪰ 0` the information they must clear:

> **Theorem 2.** For non-negative edge costs, the least-cost pair-entangled labelled program
> meeting the contract `A^T F_Π A ⪰ G_req` on `H` is the value of the convex program
>
> ```
> min   Σ_{e∈E} c_e t_e
> s.t.  A^T F A ⪰ G_req
>       diag F = 1
>       F_ij = 0                    (ij ∉ E)
>       −t_ij ≤ F_ij ≤ t_ij         (ij ∈ E)
>       t ∈ MATCH(H)
> ```
>
> and any optimal `F` is realised by the schedule of 2B–2D at exactly that cost.

*Convexity.* The Loewner constraint is linear in `F` composed with the PSD cone;
`MATCH(H)` is a polytope; the remaining constraints are linear. ∎

*The `t = |F|` reduction.* Feasibility forces `t_e ≥ |F_e|`. Replacing `t` by `|F|` keeps
`t ∈ MATCH(H)` by down-closedness (2A) and cannot increase `Σ c_e t_e` because `c ≥ 0`. So
some optimum has `t_e = |F_e|`, which is what makes the objective the cost identity of 2D.
Negative costs would break this, and the implementation rejects them.

*Solution method.* The only non-polyhedral constraint separates by an eigenvector. If `v` is
a minimal eigenvector of `A^T F A − G_req` at the current iterate, then with `a_i` the `i`-th
row of `A`,

```
Σ_e F_e · 2 (a_i·v)(a_j·v)  ≥  v^T G_req v − v^T A^T A v
```

is valid for the whole feasible set. Adding such cuts to the LP over the matching polytope
gives a semi-infinite LP solved to prescribed precision.

**Claim discipline.** This is *convex, certificate-producing optimisation to a prescribed
precision*. It is **not** a strongly polynomial algorithm, and the paper must not say
otherwise: the matching half inherits strong combinatorics, but the contract half carries a
semidefinite constraint.

## 2F — Certificates

The compiler returns exactly one of three things, each independently checkable:

1. **Feasible** — `F`, its matching decomposition, the emitted circuits, and the minimum
   eigenvalue of `A^T F A − G_req` as the satisfaction margin.
2. **Lower bound** — every LP relaxation encountered along the way uses only cuts valid for
   the true feasible set, so its optimum lower-bounds the true cost. Reported alongside the
   achieved cost as an optimality gap.
3. **Infeasible** — an LP that is infeasible under valid cuts. Since each cut and each
   structural row holds for every genuine pair-width program, LP infeasibility is a proof
   that the requested floor is unreachable at width two on `H`, and Farkas supplies the
   explicit multipliers over the finite row set.

With a finite cut budget the method may also return `iteration_limit`, which is reported as
such and never as either of the three above.

---

## Corollaries

### C1 — Recovering the `k`-producibility witness

```
max { 1ᵀ F 1 : F ∈ Q^prog_{m,k} }  =  ⌊m/k⌋ k² + r² ,     r = m − ⌊m/k⌋k
```

*Verified exactly for every `(m,k)` with `m ≤ 7`, `k ≤ 3` (15 cases); at `k = 2` the matching
polytope gives it in closed form as `m + 2⌊m/2⌋`, checked to `m = 12`.*

The established entanglement-depth criterion is the support function of the reachable body in
**one** direction. Theorem 1 gives the complete facet description of that body at `k = 2`.
Witness theory asks whether a measured scalar rules out small depth; CovQ asks whether an
entire requested contract is achievable at small width, and emits the program.

### C2 — The odd/even collective-mode ceiling

For `u_S = 1_S/√k` on a subset `S` of size `k`, width-2 compilation caps `u_S^T F u_S` at

```
k even :  2
k odd  :  2 − 1/k
```

*Proof.* `u_S^T F u_S = 1 + (2/k) Σ_{i<j∈S} F_ij ≤ 1 + (2/k) Σ_{i<j∈S} |F_ij|`. For even `k`
the largest matching inside `S` has `k/2` edges, giving `1 + (2/k)(k/2) = 2`. For odd `k` the
blossom inequality on `S` caps the sum at `(k−1)/2`, giving `1 + (k−1)/k`. Both are attained
by a maximum matching on `S` with all signs positive. ∎

*Verified for `k = 2..8`, exact to 1e-4 at every size.*

This is Edmonds' odd-set constraint stated as physics: **an odd-sized collective mode always
leaves one qubit unmatched and therefore provably cannot reach the pair-entanglement ceiling
of 2, by a deficit of exactly `1/k`.** It is the cleanest operational meaning available for
branch width, and it is a blossom effect — degree constraints alone would permit 2.

### C3 — The common-mode compilation curve

For `A = u = 1/√m` and `G_req = [γ]`, the contract is `Σ_{i<j} F_ij ≥ m(γ−1)/2`, so at unit
edge costs

```
cost*(γ)  =  m (γ − 1) / 2 ,        1 ≤ γ ≤ 1 + 2⌊m/2⌋/m ,
```

a closed-form compilation curve. A product probe reaches only `γ = 1`; the upper limit is the
`k = 2` support-function value of C1. *Verified for `m ∈ {4,5,6}` at several `γ`: compiled
cost equals the prediction to 1e-8 with a zero optimality gap, the emitted circuits meet the
floor, and `γ` above the limit is rejected with an infeasibility certificate.*

### C4 — Multi-direction contracts bind jointly

With `d ≥ 2` the Loewner constraint is genuinely active (7 cutting planes on the benchmark
below, versus 1 when `d = 1`). Three thresholds at `m = 6`, `G_req = γI`:

| contract | `γ*` |
|---|---|
| one 3-mode on `{0,1,2}` | `1.666667` = `5/3`, exactly C2 at `k = 3` |
| two **disjoint** 3-modes | `1.666667` — no interaction |
| two 3-modes **sharing one qubit** | `1.330412` |

The overlapping case is not any per-mode cap: the shared qubit's degree constraint and the
off-diagonal of `A^T F A` bind together. This is the benchmark where the compiler must
balance directions rather than concentrate everything in one GHZ-like mode.

---

## The Paper 1 application mode

Paper 1 separates the physical forward map from the programmable probe information,
`F_θ = Jᵀ W_ρ J`. CovQ's application mode is then exactly

```
F_Π = W_ρ ,     A = J ,     Jᵀ F_Π J ⪰ G_{θ,req} ,
```

i.e. *given a required information floor on the identifiable spacetime parameters, compile the
cheapest quantum probe program whose link-level QFIM meets it.* This is more natural than
extracting one exact link-level `W_star` and demanding CovQ reproduce it, because many `W_ρ`
meet the same metric-level requirement and choosing among them is precisely the freedom the
compiler exists to exploit. The paper stays quantum-computing-only because `A` is generic;
`J` is one late benchmark.

## What remains open

- **Width 3 hardness.** Expected to be a short reduction from Partition into Triangles, hence
  a corollary. The sign-parity remark in 2B suggests where the difficulty enters. Attribution
  check against the bounded-cluster partitioning literature still outstanding.
- **Polynomial-time separation.** Odd-set constraints are enumerated (`m ≤ 16`); Padberg–Rao
  is what the complexity claim rests on and is not implemented.
- **Convergence rate** of the cutting-plane scheme is not analysed. Only finite-precision
  certificates are claimed.
