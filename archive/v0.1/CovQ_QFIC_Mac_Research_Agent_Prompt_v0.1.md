# CovQ / QFIC — First Mac Research-Agent Assignment

You are beginning a novelty-gated quantum-computing research program. Do not start by training a variational circuit or writing a large benchmark harness. The first objective is to determine whether the proposed compiler abstraction is genuinely new and whether the core exact mathematics is correct.

## Working paper

**CovQ: Quantum Fisher-Information Compilation for Commuting Pauli Programs**  
*Feasibility Certificates, Minimum-Entanglement Program Synthesis, and Hardware-Aware Backends*

## Central problem

Input:

- independent commuting Pauli generators \(\mathcal G=(P_1,\dots,P_m)\);
- a target zero-mean unit-diagonal QFIM \(F_\star\);
- hardware graph \(H\);
- approximation tolerance \(\epsilon\);
- resource/noise cost model.

Output:

- a labeled randomized circuit schedule or coherent accessible-flag program \(\Pi\);
- a certificate that \(F_\Pi\) matches or approximates \(F_\star\);
- emitted and transpiled circuits;
- exact resource accounting.

## Task 0A — novelty kill-shot review

Read the full papers, not only abstracts, for at least:

1. Mahapatra and Kadiri, *Quantum State Engineering Under Multiple Expectation-Value Constraints* (QUEST), arXiv:2605.02367.
2. Wang, Tan, and Cong, *Quantum State Preparation Circuit Optimization Exploiting Don't Cares*, arXiv:2409.01418.
3. Huber and Maric, *Bernoulli Correlations and Cut Polytopes*, arXiv:1706.06182.
4. Caprara et al., *Hardness of Some Optimization Problems over Correlation Polyhedra*, arXiv:2605.02896.
5. Mirrokni et al., *Tight Bounds for Approximate Carathéodory and Beyond*, arXiv:1512.08602.
6. Du et al., *Quantifying Entanglement Dimensionality from the Quantum Fisher Information Matrix*, arXiv:2501.14595.
7. van den Berg and Temme, *Circuit Optimization of Hamiltonian Simulation by Simultaneous Diagonalization of Pauli Clusters*, arXiv:2003.13599.
8. Current cat-, graph-, and stabilizer-state compilation papers from 2024–2026.

Search explicitly for:

- target QFIM to circuit synthesis;
- observable-covariance circuit compilation;
- inverse QFIM/state-synthesis problems;
- minimum entanglement depth for prescribed QFIM;
- commuting-Pauli covariance feasibility;
- circuit ensembles or flagged programs matching target moments;
- hardware-aware synthesis over an observable-equivalence class.

Produce `NOVELTY_MATRIX.md` with columns:

- prior work;
- exact input object;
- exact output object;
- theorem/algorithm;
- hardware awareness;
- feasibility certificate;
- target-state versus target-observable semantics;
- overlap with CovQ;
- remaining differentiation;
- kill risk.

Hard stop: if a prior work already accepts commuting Pauli generators plus a target QFIM/covariance and emits hardware circuits with equivalent guarantees, report the collision before coding.

## Task 0B — line-by-line mathematical audit

Prove or refute the following candidate statements.

### A. Exact zero-mean feasibility

After Clifford normalization to independent \(Z_i\), a zero-mean unit-diagonal QFIM is feasible iff

\[
F\in\operatorname{conv}\{ss^{\mathsf T}:s\in\{\pm1\}^m\}.
\]

Clearly separate what is inherited from classical Bernoulli correlation-polytope theory from any new quantum program statement.

### B. Signed-cat primitive

For

\[
|C_s\rangle=(|s\rangle+|-s\rangle)/\sqrt2,
\]

verify exactly that \(F=ss^{\mathsf T}\) for generators \(Z_i/2\).

### C. Flagged schedule semantics

For zero-mean branches, prove that a retained classical branch label yields total experimental QFIM

\[
F_\Pi=\sum_rp_rF_r.
\]

Then analyze the coherent flag state

\[
|\Psi\rangle=\sum_r\sqrt{p_r}|r\rangle|\psi_r\rangle
\]

when generators act only on data. State exactly when its QFIM equals the branch average and what changes when branch means are nonzero.

### D. Bounded branch width

Define the exact feasible set for labeled programs whose every pure-state branch factors into blocks of size at most \(k\). Do not conflate this resource with the entanglement depth of an unlabeled mixed state.

### E. Clifford generator mapping

For independent commuting Pauli generators, state the simultaneous-Clifford normalization theorem needed by the compiler and identify the circuit cost of the map.

Deliver `THEOREM_AUDIT.md` with:

- theorem statement;
- assumptions;
- proof;
- counterexamples/edge cases;
- novelty status;
- whether it is paper-grade, background, or false.

## Task 0C — minimal exact prototype

Only after Tasks 0A and 0B pass:

1. Implement exact vertex enumeration for \(m\leq8\), quotienting \(s\sim-s\).
2. Implement LP feasibility and decomposition of \(F_\star\).
3. Emit signed-cat circuits in Qiskit or a framework-neutral circuit IR.
4. Compute QFIM directly from statevectors and compare with the target.
5. Implement infeasible-target rejection.
6. Implement labeled schedule semantics and a coherent-flag statevector test.

Pinned correctness gates:

- signed-cat QFIM relative error < 1e-12;
- exact feasible target reconstruction error < 1e-10;
- exhaustive enumeration and LP agree on all small instances;
- known infeasible targets are rejected;
- no tolerance modification at runtime.

Outputs:

- `src/` implementation;
- `tests/`;
- `results/task0_exact.json`;
- `STATUS.md` separating DONE/HOLLOW/FAILS/ABSENT;
- one concise report.

## Prohibited shortcuts

- Do not call a generic VQA result a compiler contribution.
- Do not equate unlabeled mixed-state QFI with average branch QFI.
- Do not claim the correlation polytope as new.
- Do not claim minimum entanglement depth unless the resource and state/program semantics are exact.
- Do not use spacetime as the only benchmark.
- Do not add QAOA unless a separate speedup question is justified.
- Stop and report if a theorem fails; do not patch the statement silently.
