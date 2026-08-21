# J7 — proof audit

**This is a self-audit and does not discharge J7.** The protocol asks for an *independent*
line-by-line review; that requires a second reader. What follows is the structured input
for that reader: for each theorem, the load-bearing step, what is verified numerically,
and — where they exist — the gaps I could not close from inside.

Verdict summary: **no unsupported first claim found**, four items flagged for the
independent reader, three of them about hypotheses rather than arguments.

---

## T1 — exact pair-width characterisation

`F ∈ Q^lab_{H,2}` iff `diag F = 1`, `F_ij = 0` off `E`, and `(|F_ij|) ∈ MATCH(H)`.

**Necessity.** A width-two branch partitions `V` into blocks of size ≤ 2, so its
off-diagonal support is a matching and its entries are `±1`. A labelled schedule is a
convex combination, so `(|F_ij|)` is a convex combination of matching incidence vectors,
i.e. in `MATCH(H)`. Sound.

**Sufficiency — the load-bearing step is the sign assignment, not the decomposition.**
Membership in `MATCH(H)` gives a convex decomposition of *magnitudes* into matchings. The
signs must then be realisable simultaneously. They are, for a reason specific to `k = 2`:
blocks are disjoint, so each block's sign pattern is free, and a 2-block realises
`s_is_j = ±1` either way. Every matching in the decomposition containing edge `e` is
assigned `sign(F_e)` — the same value each time, so no conflict can arise.

**Where this breaks at `k ≥ 3`, and why that is the real content.** Inside a 3-block,
`(s_is_j)(s_js_k)(s_is_k) = +1` identically, so magnitudes no longer determine
feasibility. The characterisation is not "matching polytope" boilerplate; it is the
observation that `k = 2` is exactly the width at which the sign-parity obstruction is
vacuous. The audit should keep that sentence.

*Verified:* Edmonds separation against full enumeration of the extreme points of
`Q^prog_{m,2}`, 160 instances at `m ∈ {3,4,5,6}` deliberately including infeasible ones,
zero disagreements. Blossom facets shown active (switch at exactly `1/3`, not the
degree-only `1/2`).

## T2 — cost identity

`E[N_2q] = Σ_e |F_e|`, via `t_e = |F_e|` and down-monotonicity of `MATCH(H)`.

Sound, and independently confirmed: the *measured* expected two-qubit count from emitted
schedules equals `Σ_e |F_e|` to `4.4e-16` across all eleven K3 instances — a check that
does not reuse the test written for the identity.

## T3 — compiler duality

`max ⟨Y, G_req − AᵀA⟩ − Σμ_i − Σν_S(|S|−1)/2` subject to the edge-wise bound and `Y ⪰ 0`.

**Flag 1 — constraint qualification is asserted, not proved.** Strong duality is observed
(gap `< 1e-15` at every noise level tested) but the argument requires Slater, i.e. a
strictly feasible point. The manuscript should state the qualification explicitly and say
what happens when `G_req` sits on the boundary of the reachable set. Numerically the gap
closes; that is evidence, not a proof of the hypothesis.

## T4 — block-local noisy pricing

`F^𝒩_{C,M,σ} = F^𝒩_idle + Σ_{e∈M} ΔF^𝒩_{e,σ_e}`, so pricing is a maximum-weight matching.

**Load-bearing step: three separate factorisations must hold together** — the state over
blocks, the channel over blocks, and the readout over blocks. Any one failing breaks it.
This is why the excluded list (crosstalk, correlated branch noise, route collisions,
coherent inter-block errors) is a hypothesis and not a caveat.

`max_σ` commutes with `max_M` because each edge's sign enters only its own term; disjoint.

**The pairless branch is a separate column.** It activates no pair, never waits through a
pair-preparation layer, and so carries no idle dephasing — a different reference from the
one the edge weights are measured against. Every non-empty branch has two-qubit depth
exactly one, hence the same idle wait, which is what makes the *idling* product template
the right reference for them.

*Verified:* additivity exactly `0.0`; oracle against exhaustive signed-matching
enumeration, 12 cases, worst gap `4.4e-16`. Using the idle-free reference broke it in
15 of 18 diagnostic cases, errors to `0.71`.

## T5 — block-local readout attainability

`p(x|θ) = 2^{−k}[1 + P_x cos(φ_B − A_B)]`, and `F_C^{(B)} = F_Q^{(B)} = s_Bs_Bᵀ` off the
singular set. Proved for arbitrary `k`; `k = 2,3,4` are implementation validation only.

**Flag 2 — the singular-set statement must stay first-order.** On
`φ_B − A_B ∈ πZ` the score vanishes and the model is *nonregular*; probabilities turn on
quadratically, so second-order distinguishability of phase magnitude may survive. The
claim is that the readout has zero score information and does not attain the signed local
QFIM. Any stronger phrasing is unsupported.

**A corollary the audit should not lose:** because `θ` enters a block only through the
scalar `φ_B = s_Bᵀθ_B`, the block Fisher matrix is *exactly* rank one in direction `s_B`
for **every** analyzer angle and every block-local channel — including asymmetric readout
confusion, where the outcome law is no longer a function of parity alone. An analyzer
angle can scale a block's CFI; it can never rotate it. This is what licenses the
deployable template being a scalar multiple of the matched-quadrature template.

## T6 — pinching

`Cov_{D_𝒢(ρ)}(P) = Cov_ρ(P)` and `F_Q(D_𝒢(ρ)) = 0`.

Two lines, both tight: `P_i` and `P_iP_j` are diagonal in the joint eigenbasis so no
generator moment moves; `[D_𝒢(ρ), P_i] = 0` makes the family constant in `θ`. Verified on
random commuting Pauli families: covariance preserved to `4.9e-15`, `|F_Q| ≤ 4.7e-28`.
This is the sharpest single result in the paper and the reviewer should check it first —
it is also the cheapest to check.

## T7 — width-three hardness

**Flag 3 — verified by instance, not by proof, on this side.** Partition Into Triangles is
NP-complete and the manuscript's §7.2 reduction is stated there; the prototype confirms
the *behaviour* on six graphs. I have not audited the reduction line by line, and the
independent reader must. The residual risk is not that it is false but that it is folklore
in the bounded-cluster partitioning literature — a novelty question, not a correctness one.

## T8 — estimator

Parity counts are sufficient; the parity model's Fisher matrix equals the schedule QFIM
exactly (`0.0`); the branch-conditioned MLE is consistent and efficient.

**Flag 4 — efficiency is empirical, and the reference band is a heuristic.** Whitened
efficiency lands in `[0.992, 1.039]` against a Marchenko–Pastur band of `[0.946, 1.056]`.
MP is a finite-sample *reference*, not a confidence interval, and its Wishart hypotheses
are not established here. The load-bearing claims are the exact Fisher equality and the
`1/N` scaling; the band should be presented as calibration, not inference.

The likelihood is genuinely multimodal — an arbitrarily seeded MLE inflates the variance
ratio to 10–2400×. The pilot therefore does two jobs: it phases the analyzer *and* selects
the correct periodic chart. The manuscript should say "identifies the correct local
likelihood chart, modulo the declared phase periodicity and only on the identifiable
support" and must not claim global identification of arbitrary `θ ∈ R^m`.

---

## Items the independent reader must not skip

1. **T3 Slater / constraint qualification** — stated, not proved.
2. **T7 width-three reduction** — audited only by instance here.
3. **T5 nonregularity wording** — first-order only, everywhere it appears.
4. **T8 MP band** — calibration reference, not a confidence interval.

## Citations

Not verifiable from this environment: `arxiv.org` is blocked by the network egress proxy,
so titles, authors and identifiers were confirmed only through search result metadata.
**Every citation must be re-verified against the actual records before submission**, and
`arXiv:2605.02367` in particular — the QUEST reimplementation is built from its published
method description, and the reimplementation claim depends on that description being
quoted accurately.
