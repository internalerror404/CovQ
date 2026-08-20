# Manuscript v0.2 — independent numerical audit

Covers the 75-page `CovQ_Paper_v0.2.pdf` (18 sections) and the earlier 39-page draft.

Every checkable claim in `CovQ_Certified_QFIM_Compilation_Manuscript_v0.2.pdf` that the
prototype can reach, run against `prototype/` (120 tests). Nothing here is a re-derivation
from the manuscript's own algebra; each is an independent computation.

## Verified

| Claim | Check | Result |
|---|---|---|
| **Thm 5.3** hardware-native pair-width characterisation | Edmonds separation vs brute-force enumeration of the extreme points of `Q^prog_{m,2}`, 160 random instances, `m ∈ {3,4,5,6}`, feasible and infeasible | 0 disagreements |
| **Cor 5.4** complete-graph form | blossom family at common off-diagonal `c`: switches at exactly `1/3` while degree alone permits `1/2`, whole family inside `Q_m` | exact |
| **Cor 5.5** bipartite ⇒ odd sets redundant | degree-only test vs full Edmonds on `line`, `square_grid`, `ring(even)`, 225 instances | 0 disagreements |
| **Thm 7.1** width-3 NP-hardness reduction | `max L_G` over `Q^lab_{6,3}` on 6 graphs, 3 with and 3 without triangle partitions | `= 3q` iff partition exists, all 6 |
| **Eq (65)** branch separation identity | brute force over all matchings × sign patterns vs `tr Y − c₀ + MWM(2\|Y_ij\| − c_ij)`, 12 random PSD `Y` | agreement to 1e-15 |
| **Thm 8.2** strong duality | column generation over branches with the Eq-(65) oracle, `m ∈ {4,5,6}`, all-to-all and line | primal = dual to 1e-8 in every case |
| **Eq (63)** `OPT ≤ c₀λ_max(G⋆)` | same runs | holds; **and is frequently tight — see below** |
| **§8.1** dominance collapse | `F ⪰ F⋆` with equal unit diagonals | forces `F = F⋆`, confirmed |
| **Prop 6.12** Gershgorin conditioning | `1−ρ ≤ λ_min ≤ λ_max ≤ 1+ρ` and `κ₂ ≤ (1+ρ)/(1−ρ)` over 122 pair-feasible targets, `m ∈ {3..6}` | 0 violations |
| **Ex 6.13** three-edge schedule | independent compilation of `m=3`, off-diagonals `(0.4, 0.3, −0.2)` | **reproduced exactly**: same four branches at `p = 0.4, 0.3, 0.2, 0.1`, third edge `\|Φ⁻⟩`, expected pair activations `0.900000`, emitted circuits realise `F` to 4e-16 |
| **Ex 6.14** odd-set obstruction | degree loads `0.8 ≤ 1` yet `\|F₁₂\|+\|F₁₃\|+\|F₂₃\| = 1.2 > 1` | confirmed; smallest instance forcing blossom inequalities |
| **Frozen controls** `F(±)_ij = ±0.4` | three-way classification: PSD / global / pair | both match their frozen expectations |

The Eq-(65) identity is the load-bearing one: Theorem 8.2's oracle-polynomial claim rests
entirely on it, and it is exact.

## New: an answer to K9, in closed form

Running Problem 8.1 with `c₀ = c_e = 1` produces **cost equal to the Eq-(63) product-probe
bound in every case, with no entangled branch ever pricing in**. That is kill gate K9 firing on
the manuscript's own default cost model, and it needs to be in the paper before a referee finds
it.

It is not fatal, and the boundary is sharp. A perfect-matching branch has `uᵀBu = 2` on the
common mode, so it needs `γ/2` uses at `c₀ + (m/2)c_e` each, against `γ` uses at `c₀` for the
product probe. Hence:

> **Pair entanglement is cost-optimal on a common-mode floor iff `c_e < 2c₀/m`.**

Measured by bisection, exact to 1e-3 in all 8 cases tried:

| `m` | `c₀` | measured `c_e*` | `2c₀/m` |
|---|---|---|---|
| 4 | 1.0 | 0.500000 | 0.5000 |
| 6 | 1.0 | 0.333333 | 0.3333 |
| 8 | 1.0 | 0.250000 | 0.2500 |
| 10 | 1.0 | 0.200000 | 0.2000 |
| 4 | 2.0 | 1.000000 | 1.0000 |
| 6 | 2.0 | 0.666667 | 0.6667 |
| 8 | 2.0 | 0.500000 | 0.5000 |
| 10 | 2.0 | 0.400000 | 0.4000 |

And the payoff is bounded:

> **As `c_e → 0` the cost ratio against the product probe tends to exactly `1/2`, for every
> `m`.** Verified at `m = 4, 6, 8`: ratio `0.50000`.

That constant is not a coincidence — it is the width-two collective-mode ceiling
(`uᵀFu ≤ 2` for even `k`) showing up as a cost ratio. The two results are the same fact seen
from the primal and the dual side.

**Two consequences for §8.**

1. The compiler is non-trivial exactly in the cheap-entangler regime, and the paper should say
   so with the threshold rather than leave Eq (63) as a one-sided bound. A reader who plugs in
   `c_e = c₀` gets no advantage at all.
2. The threshold scales as `1/m`: **entanglers must get cheaper in proportion to system size**
   for pair entanglement to keep paying on a common-mode contract. That is an unflattering
   scaling and belongs in Limitations, not buried.

## Where the manuscript supersedes the repository's earlier formulation

The repository's `information_floor_compile` used a normalised per-shot QFIM plus a downstream
map `A`, with the contract `AᵀFA ⪰ G_req`. The manuscript's **shot-scaled** Problem 8.1 is
better and the repository now implements both:

| | repository v0.3 (`information_floor_compile`) | manuscript §8 (`shot_scaled_floor_compile`) |
|---|---|---|
| escapes the Loewner collapse by | introducing `A` | scaling by shots — no `A` needed |
| variables | edge correlations `F_e`, magnitudes `t_e` | shot allocations `t_{M,σ}` per branch |
| dual | over edges, PSD `Y` + congestion prices `μ, ν` | over branches, PSD `Y` only |
| separation | Loewner eigenvector cuts | **max-weight matching**, Eq (65) |
| claim | certificates to prescribed precision | **oracle-polynomial**, Thm 8.2 |

The manuscript's is the better primary formulation: it removes an entire modelling object
(`A`), and it earns a strictly stronger complexity claim. The repository's `A`-based version is
retained as the *downstream-projection mode* — it answers "which per-shot correlations at
minimum entangler cost", which is a different and still useful question — but it should not be
the paper's main statement.

The repository's earlier hedge that this is "not strongly polynomial" was correct as written
but weaker than necessary. Theorem 8.2's oracle-polynomial claim is right, and the Eq-(65)
check above is the reason to believe it.

## The two shot-scaled modes are not equally solvable

The manuscript gives the shot-scaled compiler in two forms. Both are implemented; they agree
to machine precision wherever both converge, but only one is usable with a scipy-only stack.

| | branch form (`shot_scaled_floor_compile`) | homogenised edge form (`shot_scaled_edge_compile`) |
|---|---|---|
| variables | shot allocations per `(M, σ)` | `t`, edge `y`, activations `a` |
| constraints | exponentially many, priced by Eq. (65) | polynomially many linear rows + one PSD |
| **solves at `c_e = c₀`** | yes | yes, agreeing to 1e-15 |
| **solves at `c_e = 0.1` or `0.02`** | **yes** | **no — hits the cut limit every time** |

The edge form's PSD constraint becomes active on a face with a non-trivial null space as soon
as entangled edges activate, and eigenvector cutting planes stall there. Adding every violated
eigendirection per round did not help; nor did bounding all variables by the valid Eq. (63)
bound. Closing it properly needs an interior-point SDP backend, deliberately not a dependency
here.

This is a point *in favour* of the manuscript's Theorem 8.2 route: the matching separation
oracle keeps the master low-dimensional, and that is exactly why it converges where the compact
formulation does not. The repository recommends the branch form and documents the other as
limited rather than silently returning a stalled iterate.

## Not checkable here

- §11 measurement compilation (attainable classical Fisher matrices) — the prototype computes
  QFIMs, not compiled POVMs. This remains the largest unverified block.
- §12 QFIM-IR / compiler architecture and §14–15 governance: **`MANUSCRIPT_INTEGRATION.md` was
  never received.** Gates C1–C16, the canonical JSON record schema, the status vocabulary, the
  manifest and DONE-record rules, and the freeze-before-run conditions are all specified there
  and cannot be implemented against a document the repository does not have. The prototype's
  gate set stops at C12 and uses its own status vocabulary; the two must be reconciled before
  any manuscript table is populated.
- Integer shot rounding and the matrix-concentration argument of Remark 8.3.
- Anything requiring QUEST, which remains unrun; no K3 answer is available.
- Noise-aware sections; the prototype is ideal-statevector only.
