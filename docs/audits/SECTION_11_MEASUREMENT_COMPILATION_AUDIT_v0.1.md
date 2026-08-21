# Section 11 audit — measurement, noise, and what "certified" survives

Scope: §11 of `CovQ_Paper_v0.2.pdf` (75 pp), pages 36–38. Evidence:
`prototype/src/covq/measurement.py`, gates M1–M3 in `prototype/tests/test_covq.py`.
Everything below is reproduced by `pytest -k "m1_ or m2_ or m3_"`.

§11 is a **scoping** section, not a table of numbers. §11.5 states plainly that v0.2
"does not yet solve minimum-cost measurement synthesis", and that the paper therefore
uses "realizes the QFIM" in the state-family sense. That is the honest framing and this
audit does not dispute it. What it does is close the operational gap for exactly the
class the compiler emits, and then test the two warnings §11 raises.

Naming note: these gates are labelled **M1–M3** rather than continuing the C-series.
The manuscript declares gates C1–C16; the prototype implements C1–C12; the contract that
fixes the numbering (`MANUSCRIPT_INTEGRATION.md`) has not been received. Renumbering is
deferred rather than guessed.

---

## M1 — local attainability is not just true, it is exact and elementary

**Claim (§11.1).** At a fixed operating point, on the identifiable tangent support, an
appropriate local measurement attains the SLD QFIM.

**What was checked.** For every branch of a compiled schedule, a **product of
single-qubit equatorial measurements** `M(α_q) = cos α_q X + sin α_q Y` attains the
branch QFIM. No joint measurement across a Bell block is needed, and no ancilla.

The reason is structural. A cat block `B` with sign vector `s` is
`(e^{-iφ/2}|s⟩ + e^{iφ/2}|−s⟩)/√2` with `φ_B = Σ_{i∈B} s_i θ_i`; the product readout
gives block-parity law `(1 + P cos(φ_B − A_B))/2^{|B|}` with `A_B = Σ_{i∈B} α_i`, whose
Fisher information is **exactly 1** for every `A_B` off the degenerate set
`sin(φ_B − A_B) = 0`. So only one angle per block is a real degree of freedom, the
optimum is attained on an interval rather than at an isolated point, and the whole
readout compiles blockwise in time linear in the number of blocks.

Because the label is retained, the schedule's total classical Fisher matrix is
`Σ_r p_r · CFI_r = Σ_r p_r F_r = F`. For this backend the state-family reading and the
operational reading of "realizes the QFIM" therefore **coincide**.

| check | result |
|---|---|
| `‖F_r − CFI_r‖₂` per branch, `m = 3..6`, θ = 0 and generic | `< 2.3e-15` |
| Loewner domination `CFI_r ⪯ F_r` | holds in every case |
| schedule-level `‖F − CFI‖₂`, task0 suite × 3 operating points | 48 cases, 0 failures, worst `2.2e-15` |
| block width `k = 2, 3, 4` | attains; worst gap `3.6e-15` |

Attainment is **not** special to pair width — it holds for any cat-block schedule. That
is worth knowing, because it means §11.1 is not what distinguishes `k = 2`; the matching
polytope is.

## M2 — the operating-point qualification is load-bearing

**Claim (§11.1(ii)).** The optimal measurement can depend on the operating point.

It does, and the failure of a *fixed* readout is not graceful. The unmatched all-X
readout attains `F` only off the union of hyperplanes

  `H = { θ : Σ_{i∈B} s_i θ_i ≡ 0 (mod π) for some block B of some branch }`.

`H` has measure zero, which makes the qualification sound minor. It is not, because `H`
contains the two operating points anyone would actually pick:

- **θ = 0.** Every block phase vanishes, and the fixed readout returns *identically zero*
  Fisher information — on every branch, for every instance tested. Not degraded: zero.
- **the uniform ray θ = c·1.** As soon as one branch carries a negatively signed edge,
  that block has `φ_B = c(s_i + s_j) = 0` for all `c`, and its contribution is lost at
  every point of the ray.

Worked fixture (`m = 4`, seed 11; branches `w=0.2` all-singletons and `w=0.8` blocks
`{0,2},{1,3}` with `s = (+,+,+,−)`):

| θ | compiled `‖F−CFI‖₂` | fixed-X `‖F−CFI‖₂` | tr F | tr CFI (fixed X) |
|---|---|---|---|---|
| `0` | 7.3e-16 | 1.800 | 4.000 | **0.000** |
| uniform `0.3` | 5.4e-16 | 1.600 | 4.000 | 2.400 |
| generic | 1.2e-15 | 0.000 | 4.000 | 4.000 |

The `2.400` is exactly `4 − 0.8×2`: the mixed-sign block `{1,3}` contributes nothing,
and nothing else is lost. The loss is fully accounted for by the sign pattern.

**Consequence for the compiler.** A CovQ output package that emits `(probe schedule,
parameter circuit)` without the matched readout has, at the natural operating point,
emitted a program with no information in it. §12.1's output tuple already lists
`Readout`; this is the argument for why that slot cannot be optional or defaulted.

## M3 — Eq (110) implemented; the covariance surrogate quantified

**Claim (§11.3).** For a noisy `ρ`, the SLD QFIM is not generally `4 Cov_ρ(G_i,G_j)`;
Eq (110) is the right object; the two must not be conflated by relabeling.

Eq (110) is implemented and validated two ways.

1. **Reduction.** On pure states, Eq (110) equals the covariance formula and equals
   `qfim_from_statevector` to `4.4e-16` over 20 states, `m = 1..4`.
2. **Closed form.** For the depolarised single-qubit phase probe with visibility
   `v = 1 − 4p/3`, Eq (110) returns exactly `v²` while `4 Cov` stays pinned at 1:

| `p` | `v` | Eq (110) | `v²` | Eq (109) | inflation |
|---|---|---|---|---|---|
| 0.00 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.00 |
| 0.05 | 0.933333 | 0.871111 | 0.871111 | 1.000000 | 1.15 |
| 0.10 | 0.866667 | 0.751111 | 0.751111 | 1.000000 | 1.33 |
| 0.25 | 0.666667 | 0.444444 | 0.444444 | 1.000000 | 2.25 |
| 0.50 | 0.333333 | 0.111111 | 0.111111 | 1.000000 | 9.00 |

The error is unbounded, not a correction. Over 80 random mixed states (`m = 2, 3`) the
surrogate never fell below the true QFIM, so it is at least a valid upper bound — but an
arbitrarily loose one.

### The sharp case: the Bell primitive is invisible to its own dominant noise

Every CovQ generator is diagonal in `Z`. Therefore `Z` dephasing **commutes with the
entire generator family and leaves the observable covariance matrix exactly invariant** —
including the edge correlator `F_e` the compiler was asked to hit — while the SLD QFIM
decays to zero. On a single Bell block at θ = (0.4, 0):

| `p` | tr CFI (compiled readout) | tr QFIM, Eq (110) | tr `4 Cov`, Eq (109) | `⟨Z₀Z₁⟩` |
|---|---|---|---|---|
| 0.00 | 2.000000 | 2.000000 | 2.000000 | 1.000000 |
| 0.10 | 0.415221 | 0.819200 | 2.000000 | 1.000000 |
| 0.20 | 0.106476 | 0.259200 | 2.000000 | 1.000000 |
| 0.35 | 0.006149 | 0.016200 | 2.000000 | 1.000000 |
| 0.50 | 0.000000 | **0.000000** | 2.000000 | 1.000000 |

At `p = 0.5` the program carries **no information whatsoever** and the covariance table
is bit-for-bit what it was at `p = 0`. The Loewner chain
`CFI(declared readout) ⪯ QFIM(110)` holds at every row.

This strengthens §11.3(a) rather than contradicting it. The paper already says a noisy
backend must report both matrices. The finding is *why* the "at least two" is not
belt-and-braces: for the Bell primitive specifically, report (a) alone and the report is
not merely imprecise, it is completely insensitive to the failure mode that matters. I
would state this in §11.3 explicitly, because a reader who skims will assume the
covariance table degrades along with the state, and here it provably does not.

---

## What this does and does not settle

**Settled.** For pure cat-block schedules, "certified" survives the move from state
family to emitted readout: the readout is a product of single-qubit measurements,
compilable in linear time, attaining `F` exactly, with the matched-readout requirement
now shown to be non-negotiable at the natural operating points.

**Not settled, and still the largest open block.**
- **Minimum-cost** measurement synthesis. M1 exhibits *an* attaining readout; it says
  nothing about the cheapest one under a declared measurement-setting cost. §11.5's
  frontier stands.
- The estimator layer. `(…, POVM, estimator)` in Eq (112) — attaining the CFI at a point
  is not the same as an estimator that reaches it at finite shots, and near `H` the
  outcome law is near-deterministic, so the finite-sample behaviour there needs its own
  treatment even though the CFI is exactly 1 up to the boundary.
- Noisy *compilation*: §11.4's objective (111) is not implemented. M3 supplies the
  measuring stick (Eq 110 plus the declared-readout CFI) but no compiler optimises it.
- Readout error and idle noise are not modelled; only depolarising and dephasing are.
