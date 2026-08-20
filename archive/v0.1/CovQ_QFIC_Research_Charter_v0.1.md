# CovQ / QFIC Research Charter v0.1

## Working title

**CovQ: Quantum Fisher-Information Compilation for Commuting Pauli Programs**  
*Feasibility Certificates, Minimum-Entanglement Program Synthesis, and Hardware-Aware Backends*

## Paper identity

This is a quantum-computing paper. Its primary object is a compiler:

\[
(\mathcal G,F_\star,H,\epsilon,\mathcal C)
\longmapsto
\Pi,
\]

where:

- \(\mathcal G=(P_1,\ldots,P_m)\) is a commuting Pauli generator program;
- \(F_\star\) is a target pure-state QFIM or observable covariance specification;
- \(H\) is a hardware coupling graph;
- \(\epsilon\) is a target error;
- \(\mathcal C\) is a declared resource/noise cost model;
- \(\Pi\) is a labeled or coherently flagged circuit program whose realized QFIM approximates \(F_\star\).

Spacetime tomography is one late benchmark only. The paper must remain coherent if all spacetime material is deleted.

## Provisional thesis

Traditional state preparation specifies a complete statevector or target unitary. CovQ specifies only a quantum information geometry. The compiler may choose any state, ensemble of flagged states, entanglement partition, and hardware mapping that realizes the requested QFIM within tolerance.

The strongest intended contribution is not generic expectation matching. It is the combination of:

1. exact QFIM feasibility for commuting Pauli programs;
2. correct labeled/coherent program semantics for additive branch QFI;
3. minimum branch-entanglement-width compilation;
4. sparse approximate decomposition with certificates;
5. generator-aware Clifford normalization;
6. hardware- and noise-aware choice among all QFIM-equivalent programs.

## Novelty boundary

### Already occupied

- Generic expectation-value-constrained pure-state synthesis (QUEST, arXiv:2605.02367).
- State-preparation optimization exploiting partial specification or “don’t cares” (arXiv:2409.01418).
- General full-state preparation and unitary synthesis.
- Cat-, graph-, and stabilizer-state circuit synthesis.
- Classical Bernoulli correlation/cut-polytope characterization (arXiv:1706.06182).
- QFIM-based entanglement witnessing and entanglement-dimensionality criteria (arXiv:2501.14595).
- Generic approximate Carathéodory sparsification (arXiv:1512.08602).

### Provisional white space to verify

- A QFIM-specific compiler intermediate representation for commuting Pauli programs.
- Exact feasibility/infeasibility certificates connected to emitted quantum programs.
- Minimum branch-entanglement-width synthesis of a target QFIM.
- Hardware-aware optimization over the full QFIM-equivalence class rather than one prescribed state.
- Labeled versus coherently flagged backend semantics and resource tradeoffs.
- Target-specific support/resource bounds stronger than generic approximate Carathéodory.

No “first” claim is allowed until the novelty kill-shot review is complete.

## Formal model

Let \(P_i\) be independent commuting Pauli observables with \(P_i^2=I\), and define

\[
U_{\boldsymbol\theta}
=
\exp\!\left[-\frac{i}{2}\sum_{i=1}^m\theta_iP_i\right].
\]

For a pure probe \(|\psi\rangle\),

\[
F_{ij}(|\psi\rangle)
=
\operatorname{Cov}_{\psi}(P_i,P_j).
\]

The initial exact theory should focus on the zero-mean unit-diagonal sector

\[
\langle P_i\rangle=0,
\qquad
F_{ii}=1.
\]

After simultaneous Clifford diagonalization, study \(P_i=Z_i\). For

\[
|\psi\rangle=\sum_{z\in\{\pm1\}^m}c_z|z\rangle,
\quad p_z=|c_z|^2,
\]

zero means give

\[
F=\mathbb E_p[zz^{\mathsf T}].
\]

Hence the exact target set is the symmetric Bernoulli correlation polytope

\[
\mathcal Q_m
=
\operatorname{conv}\{ss^{\mathsf T}:s\in\{\pm1\}^m\}.
\]

For a sign vector \(s\), the signed cat primitive

\[
|C_s\rangle
=
\frac{|s\rangle+|-s\rangle}{\sqrt2}
\]

realizes the extreme QFIM \(ss^{\mathsf T}\).

## Program semantics

A compiler output is a flagged program

\[
\Pi=\{(p_r,U_r,r)\}_{r=1}^q.
\]

Two valid backends are:

1. **Labeled schedule:** sample branch \(r\), retain the label in the classical record, and run \(U_r\).
2. **Coherently flagged program:** prepare
   \[
   \sum_r\sqrt{p_r}|r\rangle_{\mathrm{flag}}|\psi_r\rangle_{\mathrm{data}},
   \]
   retain the flag through readout, and let parameter generators act only on data.

For zero-mean branches, the program QFIM is

\[
F_\Pi=\sum_rp_rF_r.
\]

An unlabeled mixed state is not interchangeable with this program and must not be used silently.

## Compiler resource hierarchy

For a branch partition \(\pi_r=\{B\}\), define branch width

\[
w_r=\max_{B\in\pi_r}|B|,
\qquad
w(\Pi)=\max_rw_r.
\]

The compiler resource is **branch entanglement width**, not an unqualified claim about the entanglement depth of one mixed state.

Define

\[
\mathcal Q_{m,k}^{\mathrm{prog}}
=
\operatorname{conv}
\left\{
\bigoplus_{B\in\pi}F_B:
\max_{B\in\pi}|B|\le k,
\;F_B\in\mathcal Q_{|B|}
\right\}.
\]

Then

\[
k_\epsilon(F_\star)
=
\min\{k:\operatorname{dist}(F_\star,\mathcal Q_{m,k}^{\mathrm{prog}})\le\epsilon\}.
\]

This is the central compiler resource target.

## Target theorem program

### Theorem A — Exact zero-mean QFIM feasibility

For independent commuting Pauli generators, a unit-diagonal zero-mean target is realizable iff it belongs to the corresponding sign-correlation polytope after Clifford normalization.

Deliverables:

- exact statement;
- generator-dependence for nontrivial commuting Pauli relations;
- infeasibility certificates;
- pure-state and flagged-program realizations.

### Theorem B — Flagged program normal form

Every feasible target admits a finite signed-cat schedule; every such schedule has exactly the convex-combination QFIM when branch labels are retained. Establish the coherent-flag equivalent under the stated zero-mean hypothesis.

### Theorem C — Bounded-width feasibility hierarchy

Characterize \(\mathcal Q_{m,k}^{\mathrm{prog}}\) for the declared instruction set and show that the minimum feasible \(k\) is an exact compiler resource.

### Theorem D — Sparse approximation and downstream stability

If \(\|F_\Pi-F_\star\|_{\mathrm{op}}\le\epsilon\), then for every downstream linear map \(A\),

\[
\|A^{\mathsf T}(F_\Pi-F_\star)A\|_{\mathrm{op}}
\le
\epsilon\|A\|_{\mathrm{op}}^2.
\]

Seek support bounds depending on structure, stable rank, block pattern, or target spectrum—not merely generic Carathéodory.

### Theorem E — Complexity boundary

Document exact membership/decomposition hardness and identify tractable regimes:

- fixed \(m\);
- bounded branch width \(k\);
- chordal/block targets;
- low support;
- approximation with certified error.

## Kill criteria

The paper stops or pivots if any of the following occurs:

1. A prior method already accepts commuting Pauli generators plus a target QFIM/covariance and emits equivalent hardware circuits with feasibility certificates.
2. The exact theory reduces entirely to citing the correlation polytope and preparing ordinary GHZ states.
3. No new theorem survives beyond generic expectation targeting or approximate Carathéodory.
4. CovQ does not materially improve gate count, depth, settings, compile time, or target error over QUEST/generic moment matching.
5. Freedom to choose among QFIM-equivalent states yields no advantage over compiling one fixed full state.
6. The bounded-width hierarchy produces no useful or provable resource tradeoff.
7. Only spacetime-derived targets work.
8. Only a VQA result survives.

## Paper contribution checklist

A viable paper needs at least:

- one exact feasibility/normal-form theorem;
- one nontrivial compiler-resource theorem or bound;
- one exact or certified algorithm;
- one scalable approximation algorithm;
- hardware-aware emitted circuits and measured resource accounting;
- generic, structured, infeasible, and application benchmark families;
- a direct comparison to generic expectation targeting;
- a negative result or regime boundary where product/small-width programs are preferable.

## Working abstract v0.1

Quantum compilers traditionally receive a target unitary, statevector, or circuit and optimize a fixed semantic object. We study a weaker and more flexible specification: a target quantum Fisher information matrix for a commuting Pauli parameter program. The compiler may choose any quantum state or flagged circuit schedule whose information geometry matches the target, enabling optimization over an equivalence class of states rather than one prescribed wavefunction. For the zero-mean unit-diagonal sector, we connect exact feasibility to a sign-correlation polytope and give signed-cat program semantics for its extreme points. We formulate minimum branch-entanglement-width and minimum-setting synthesis, develop exact and approximate compilation algorithms, and map the resulting programs to restricted hardware graphs under gate, depth, routing, and noise costs. The evaluation compares structure-aware QFIM compilation with generic expectation-value targeting, variational moment matching, and full-state preparation across feasible, infeasible, structured, Clifford-transformed, and application-derived targets. The central question is whether observable-level information geometry can be compiled more cheaply than a complete quantum state while retaining certified downstream performance.

## Paper structure

1. Introduction
2. Related work and novelty boundary
3. Commuting Pauli programs and QFIM-IR
4. Exact feasibility and flagged normal forms
5. Minimum branch-entanglement-width synthesis
6. Exact, sparse, and approximate algorithms
7. Generator-aware Clifford normalization
8. Hardware-aware circuit backend
9. Noise-aware compilation
10. Experimental methodology
11. Results
12. Application benchmark from Lorentzian tomography
13. Limitations and open problems
14. Conclusion

## Scope exclusions for v1

- noncommuting generator QFIMs;
- arbitrary mixed-state QFI compilation;
- universal quantum speedup claims;
- QAOA as the primary solver;
- fault-tolerant QEC synthesis;
- spacetime as the main paper object.
