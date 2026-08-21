# Block-local readout: attainability, singularity, and information death

Proved for arbitrary block width `k`. The `k = 2, 3, 4` runs in
`prototype/tests/test_covq.py` are independent implementation validation, **not** the
proof. Gate aliases: C10a, C10b, C10c, C15a–c (see `docs/audits/`).

## Setup

Commuting Hermitian generators `P_1..P_m`, encoding `U_θ = exp[−(i/2) Σ_i θ_i P_i]`,
per-shot QFIM `F_ij = ⟨P_iP_j⟩ − ⟨P_i⟩⟨P_j⟩`. Let `B` be a cat block of size `k` with
sign vector `s_B ∈ {±1}^k` and **effective block phase**

  `φ_B(θ) = s_Bᵀ θ_B`.

The encoded signed-cat state is

  `|C_{s_B}(θ)⟩ = ( e^{−iφ_B/2}|s_B⟩ + e^{+iφ_B/2}|−s_B⟩ ) / √2`.

Measure each `i ∈ B` in the local equatorial basis
`|x_i; α_i⟩ = (|0⟩ + x_i e^{iα_i}|1⟩)/√2`, `x_i ∈ {±1}`. Write the measured parity and
the **aggregate analyzer phase**

  `P_x = Π_{i∈B} x_i`,  `A_B = s_Bᵀ α_B`,  `δ = φ_B − A_B`.

## Proposition 1 (local equatorial probability law)

  `p(x | θ) = 2^{−k} [ 1 + P_x cos(φ_B − A_B) ]`.

*Proof.* `⟨x;α|s_B⟩ = 2^{−k/2} Π_i (…)`; collecting the two overlaps gives
`⟨x;α|C⟩ = 2^{−(k+1)/2}( e^{−iφ_B/2} + P_x e^{−iA_B} e^{+iφ_B/2} )` up to a global phase,
and squaring gives the claim. ∎

Only the parity carries parameter dependence; the remaining `k − 1` outcome bits are
ancillary uniform randomness. Verified to `2.5e-16` for `k = 2, 3, 4`.

## Theorem 2 (ancilla-free block-local attainability)

Off the singular set `δ ∈ πZ`,

  `F_C^{(B)} = s_B s_Bᵀ = F_Q^{(B)}`.

*Proof.* `∂_{θ_i} p = −s_i P_x sin(δ) 2^{−k}`. Splitting the `2^k` outcomes into the two
parity classes,

  `Σ_x (∂_φ p)²/p = (sin²δ / 2)[ 1/(1+cos δ) + 1/(1−cos δ) ] = sin²δ / sin²δ = 1`,

so the scalar parity information is exactly 1, and the chain rule `∂φ_B/∂θ_i = s_i`
lifts it to `s_B s_Bᵀ`. That equals the block QFIM. ∎

**Cost.** All but one qubit may be measured in `X`; a single pivot carries the analyzer
angle. Hence `N^readout_2q = 0` and `N^readout_anc = 0` — optimal in both. No joint Bell
measurement, no entangling readout, no ancilla. This supersedes the manuscript's
singleton / positive-Bell / negative-Bell case analysis, which is a special case.

## Corollary 3 (retained-label schedule attainability)

If branch `r` factors into disjoint signed-cat blocks, `F^{(r)} = ⊕_{B∈π_r} s_{r,B}s_{r,B}ᵀ`,
then the tensor product of block readouts has `F_C^{(r)} = F^{(r)}`, and retaining the label

  `F_C^schedule = Σ_r p_r F_C^{(r)} = Σ_r p_r F^{(r)} = F_Π`. ∎

## Proposition 4 (singular hyperplanes)

The first-order failure set of a fixed analyzer choice is
`{ θ : φ_B(θ) − A_B ∈ πZ for some block B of some branch }`. There
`∂_{θ_i} p(x|θ) = 0` for every outcome, so the score vanishes and the local Fisher matrix
is zero.

**Wording, deliberately careful.** The model is *nonregular* there, not
information-free in every sense: probabilities turn on quadratically off the hyperplane,
so nonregular second-order distinguishability of the phase *magnitude* may survive. The
claim is the first-order one — the readout has zero score information and does not attain
the signed local QFIM. That is what the compiler depends on.

For the all-X choice `A_B = 0` the set is `s_Bᵀθ_B ∈ πZ`, which contains `θ = 0` for every
block, and the entire uniform ray `θ = c·1` for any block with `s = (1,−1)`.

## Proposition 5 (matched quadrature, and why it must be solved not searched)

Turning the pivot angle alone, `⟨P_B⟩(α) = c cos α + d sin α` with `c`, `d` the parity
expectations at `α = 0` and `α = π/2`. Hence two evaluations determine the fringe, and

  `α* = atan2(−c, d)`  gives  `η_ro := |sin δ| = 1`.

Equivalently `s_Bᵀα_B = s_Bᵀθ_ref,B ± π/2 (mod π)`.

A grid argmax cannot substitute. On a pure block *every* non-degenerate angle ties at
Fisher information 1, so the tie-break returns an arbitrary detuning — free at `v = 1`,
and costly the instant visibility drops. This was a real defect in the first
implementation, invisible to every pure-state test and caught only by the noisy column.

## Proposition 6 (noisy detuned analyzer)

For a dephased cat block of fringe visibility `v`,

  `F_C(v, δ) = [ v² sin²δ / (1 − v² cos²δ) ] · s_B s_Bᵀ`,

so quadrature gives `F_C = v² s_Bs_Bᵀ = F_Q`. Verified against numerics to `3.8e-15` over
120 `(k, p, φ, A)` points.

Consequences: at `v = 1` every non-degenerate analyzer phase attains; at `v < 1` the
analyzer phase matters continuously and quadrature remains optimal; at `v = 0` the
information is zero **for every readout**, not for this one. There is therefore no
measurement gap to explain in the dephasing table.

## Proposition 7 (covariance-preserving information death)

Let `D_𝒢` be the pinching onto the joint eigenbasis of the commuting family. For every `ρ`,

  `Cov_{D_𝒢(ρ)}(P) = Cov_ρ(P)`   and   `F_Q(D_𝒢(ρ)) = 0`.

*Proof.* Pinching preserves all diagonal blocks in that eigenbasis, and `P_i`, `P_iP_j` are
diagonal there, so every first and second generator moment is unchanged. And
`[D_𝒢(ρ), P_i] = 0` gives `U_θ D_𝒢(ρ) U_θ† = D_𝒢(ρ)`: the state family is constant, so its
QFIM vanishes identically. ∎

This is strictly stronger than "mixed-state QFI is not covariance". Generator-basis
covariance can be **exactly** unchanged while all parameter information is destroyed.
Verified for random commuting Pauli families on 2–4 qubits: covariance preserved to
`4.9e-15`, `|F_Q| ≤ 4.7e-28`.

For the Bell block under independent phase flips of probability `p` per qubit,
`v = (1−2p)²` and `tr F_Q = 2(1−2p)⁴`; at `p = 0.2` this is `2(0.6)⁴ = 0.2592`.

## What is *not* claimed

Attainability is block-local for the **entire** signed-cat hierarchy. The special status
of pair width comes from exact matching geometry and certified optimization, not from
readout compatibility. Nothing here narrows to `k = 2`.
