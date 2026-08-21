# CovQ Repository Handoff v0.3

## Controlling decision

This handoff supersedes the v0.2 repository plan. The paper is no longer organized around exact target-QFIM reconstruction with an information-floor mode attached later. The controlling object is now a **Fisher-information contract compiler**:

\[
(\mathcal G,H,A,G_{\rm req},\mathcal C)
\longmapsto
(\Pi,\text{certificate},\text{circuits},\text{readout}),
\]

with primary optimization

\[
\min_{\Pi\in\mathfrak P^{\rm lab}_{H,2}}
\operatorname{Cost}_{\mathcal C}(\Pi)
\quad\text{subject to}\quad
A^{\mathsf T}F_\Pi A\succeq G_{\rm req}.
\]

The exact pair-width QFIM characterization is the mathematical engine that makes this compiler constructive and certifiable. It must not be presented in code, documentation, or reports as a second independent project.

The intended end-to-end object is:

```text
downstream Fisher floor
    -> exact pair-width QFIM body on H
    -> minimum-cost convex optimization
    -> matching decomposition
    -> signed native Bell circuits
    -> primal/dual/infeasibility certificate
```

Exact matching `F_Pi = F_target` remains implemented, but only as a secondary compiler mode, theorem canary, nearest-projection control, and interface to state-preparation baselines.

## Paper identity

**CovQ: Certified Fisher-Information Contract Compilation for Commuting Pauli Programs**  
*Exact Pair-Width Geometry and Minimum-Cost Hardware-Native Bell Schedules*

The repository must preserve the following claim boundary:

- Classical optimal design, information-matrix constraints, semidefinite sensor selection, the matching polytope, maximum-weight matching, and bounded-cluster combinatorics are inherited tools.
- CovQ contributes their identification with an exact **achievable QFIM body** for retained-label pair-width quantum programs, sign-complete circuit realization, an information-floor compiler interface, and certificate-bearing hardware program emission.
- No strongly polynomial claim is permitted for the complete information-floor problem. The pair-width geometry has polynomial separation; the compiler contains a semidefinite constraint and should be described as convex and solvable to prescribed precision with verifiable certificates.
- The width-three reduction is retained to explain why width two is privileged. Attribute the classical source problem and do not sell the reduction as a new complexity-theory result.

## Scientific disposition

| Item | Status | Repository role |
|---|---|---|
| Downstream information-floor compiler | Main problem | Primary CLI, API, experiments, and paper results |
| Hardware-native pair-width characterization | Main structural theorem | Exact feasible body used by compiler |
| `F,t` matching-constrained SDP | Main optimization theorem | Primary finite formulation |
| Signed-matching branch master | Main optimization theorem | Column-generation form and emitted schedule |
| Dual lower bound and matching pricing | Main certificate theorem | Optimality certificate |
| PSD support-function witness | Main certificate theorem | Infeasibility certificate |
| Common-mode closed-form cost curve | Main analytic benchmark | C13 and primary figure |
| Multidirectional information floors | Main empirical benchmark | Tests genuine compiler freedom |
| Exact target-QFIM reconstruction | Secondary mode | Correctness, projection, and baseline control |
| Full sign-correlation polytope | Background | Small-instance control only |
| Signed-cat and Bell primitives | Background | Circuit primitive and unit test |
| Coherent flags | Separate resource model | Semantics appendix/control, not pair-width backend |
| Width-two/width-three boundary | Contextual theorem | Explains tractable regime |
| Lorentzian light-ray application | Late case study | `A=J`, removable without weakening core paper |

## Primary mathematical contract

### Realized QFIM body

For physical generators `Z_i/2` and direct native pair operations on a graph `H=(V,E)`, the primary backend realizes

\[
\mathcal Q^{\rm lab}_{H,2}
=
\left\{
F:
\operatorname{diag}F=\mathbf 1,
\ F_{ij}=0\; (ij\notin E),
\ (|F_{ij}|)_{ij\in E}\in\operatorname{MATCH}(H)
\right\}.
\]

A branch is a signed matching:

- every active edge carries `Phi+` for sign `+1` or `Psi+` for sign `-1`;
- unmatched vertices carry `|+>`;
- all active Bell blocks in one branch are disjoint and may be prepared in parallel;
- the retained branch label selects the corresponding measurement/estimator metadata.

### Primary finite formulation

For nonnegative native-edge costs `c_e`, solve

\[
\begin{aligned}
\min_{F,t}\quad & \sum_{e\in E}c_et_e\\
\text{s.t.}\quad & A^{\mathsf T}FA\succeq G_{\rm req},\\
&\operatorname{diag}F=\mathbf 1,\\
&F_{ij}=0\quad(ij\notin E),\\
&-t_{ij}\le F_{ij}\le t_{ij}\quad(ij\in E),\\
&t\in\operatorname{MATCH}(H).
\end{aligned}
\]

Because `MATCH(H)` is downward closed, a nonnegative-cost optimum may be represented with `t_e = abs(F_e)`. The objective is then exactly the expected native Bell-pair preparation cost per shot for the direct pair-block instruction set.

### Branch-master form

Let `b=(M,sigma)` index a signed matching branch, `B_b` be its unit-diagonal QFIM, and `c(b)=sum_{e in M} c_e`. Solve

\[
\begin{aligned}
\min_{p\ge0}\quad & \sum_b p_bc(b)\\
\text{s.t.}\quad & \sum_bp_bA^{\mathsf T}B_bA\succeq G_{\rm req},\\
&\sum_bp_b=1.
\end{aligned}
\]

The implementation must cross-check the `F,t` and branch-master formulations on every graph small enough for explicit signed-matching enumeration.

### Dual and pricing

The branch-master dual is

\[
\begin{aligned}
\max_{Y\succeq0,\beta}\quad & \langle G_{\rm req},Y\rangle-\beta\\
\text{s.t.}\quad &
\langle AYA^{\mathsf T},B_b\rangle-\beta\le c(b)
\quad\text{for every signed matching branch }b.
\end{aligned}
\]

For `Q=A Y A^T`, the most violated branch is obtained from

\[
\operatorname{tr}Q-\beta+
\max_{M\in\operatorname{Match}(H)}
\sum_{e\in M}(2|Q_e|-c_e).
\]

Thus the pricing/separation oracle is maximum-weight matching. The code must return:

- the matching;
- selected signs `sign(Q_e)`;
- raw violation;
- tolerance-adjusted violation;
- graph/cost hash;
- deterministic tie-breaking metadata.

### Infeasibility witness

For feasibility independent of cost, define

\[
h_H(Q)=\operatorname{tr}Q+
2\max_{M\in\operatorname{Match}(H)}\sum_{e\in M}|Q_e|.
\]

A PSD matrix `Y` certifies that a requested floor is impossible when

\[
\langle G_{\rm req},Y\rangle
>
h_H(A Y A^{\mathsf T}).
\]

Every repository-level `INFEASIBLE` conclusion must include a serialized witness `Y`, both sides of the inequality, the maximum-weight matching used to evaluate the support function, and a positive certified margin greater than the frozen tolerance.

## Proof and implementation order

Do not start broad benchmarking until steps 2A-2F are separately closed.

### 2A - Necessity

Prove and test that every retained-label width-two branch yields a signed matching:

- two-qubit blocks contribute only one native edge each;
- disjointness gives a matching;
- arbitrary nonextremal two-qubit correlations lie below the matching incidence vector;
- convex aggregation places `abs(F_edges)` in `MATCH(H)`.

Required outputs:

```text
docs/proofs/2A_NECESSITY.md
tests/test_pair_necessity.py
results/proof_checks/2A_necessity.json
```

### 2B - Sign-complete sufficiency

Prove and test that edge signs create no additional compatibility constraints:

- decompose `abs(F_edges)` into matchings;
- apply the sign of `F_e` independently on every active edge;
- emit `Phi+` or `Psi+` accordingly;
- reconstruct every matrix entry, not only a scalar witness.

Required outputs:

```text
docs/proofs/2B_SIGN_SUFFICIENCY.md
tests/test_sign_sufficiency.py
results/proof_checks/2B_sign_sufficiency.json
```

### 2C - Hardware graph theorem

Prove the theorem directly on general `H`, not only `K_m`:

- nonedges must vanish in direct mode;
- bipartite graphs use degree inequalities only;
- nonbipartite graphs require odd-set inequalities;
- routing is a separately named backend and may not be hidden in theorem tests.

Required outputs:

```text
docs/proofs/2C_HARDWARE_GRAPH.md
tests/test_hardware_matching_body.py
results/proof_checks/2C_hardware_graph.json
```

### 2D - Constructive emission and cost

Implement:

- matching-polytope decomposition;
- support reduction to at most `|E|+1` branches when requested;
- signed Bell branch IR;
- emitted circuits;
- branch-conditioned measurement metadata;
- expected-cost identity from parsed circuits.

Required outputs:

```text
docs/proofs/2D_EMISSION_AND_COST.md
src/covq/decompose.py
src/covq/circuits.py
src/covq/measurements.py
tests/test_emission_cost.py
results/proof_checks/2D_emission_cost.json
```

### 2E - Information-floor convex program

Close, in this order:

1. implement and verify the `F,t` primal;
2. implement the explicit branch master on tiny graphs;
3. prove numerical objective equivalence;
4. verify `t=abs(F_edges)` at an optimum under nonnegative costs;
5. emit a schedule from the primal `F`;
6. verify the emitted schedule satisfies `A^T F A >= G_req`;
7. reproduce the analytic common-mode cost curve.

Required outputs:

```text
docs/proofs/2E_INFORMATION_FLOOR.md
src/covq/contracts.py
src/covq/floor_sdp.py
src/covq/branch_master.py
tests/test_floor_equivalence.py
tests/test_common_mode_curve.py
results/proof_checks/2E_information_floor.json
results/common_mode.json
```

### 2F - Certificates

Implement:

- dual lower bounds;
- matching pricing/separation;
- primal-dual residuals;
- relative gap;
- PSD support-function infeasibility witnesses;
- certificate verification independent of the solver object that produced them.

Required outputs:

```text
docs/proofs/2F_CERTIFICATES.md
src/covq/dual.py
src/covq/certificates.py
tests/test_floor_duality.py
tests/test_floor_infeasibility.py
results/proof_checks/2F_certificates.json
```

Only after 2E and 2F pass may the repository describe CovQ as a certified minimum-cost information-floor compiler.

## Code architecture

```text
src/covq/
  qfim.py                 # pure-state and program QFIM validation
  flags.py                # labeled, coherent, and dephased semantics
  contracts.py            # A, G_req, support, product slack, contract margin
  matching.py             # membership, separation, support, max-weight matching
  signed_branches.py      # signed matching QFIMs and deterministic branch IDs
  floor_sdp.py            # F,t primal
  branch_master.py        # explicit/column-generated branch formulation
  dual.py                 # Y,beta dual and pricing
  certificates.py         # standalone verification of feasible/cost/infeasible certs
  decompose.py            # matching decomposition and support reduction
  circuits.py             # Bell branch emission and routing-separated backend
  measurements.py         # local branch readout and classical-Fisher checks
  exact_target.py         # secondary exact-match/projection modes
  full_sector.py          # small-instance background controls only
  clifford.py             # separately reported logical-frame extension
  noise.py                # noisy covariance/QFI/readout Fisher separation
  schema.py               # canonical artifact validation
```

Do not overload one function with both direct and routed semantics. A result record must state:

```text
backend_mode = direct_native | routed | clifford_normalized
```

and only `direct_native` inherits the main pair-width theorem without additional qualifications.

## Certificate interface

The compiler returns exactly one top-level disposition:

### `FEASIBLE_OPTIMAL`

Contains:

- `F`, `t`;
- `A`, `G_req`, and all input hashes;
- downstream surplus spectrum;
- matching decomposition;
- emitted circuit and measurement hashes;
- parsed expected cost;
- dual `Y,beta`;
- primal/dual objectives;
- all residuals and relative gap;
- standalone verifier status.

### `FEASIBLE_BOUND_ONLY`

Used only when a feasible emitted program exists but optimality has not been certified. It must not be described as least cost.

Contains:

- feasible program certificate;
- current upper bound;
- current dual lower bound if available;
- unresolved gap and reason.

### `INFEASIBLE`

Contains:

- PSD witness `Y`;
- `lhs=<G_req,Y>`;
- `rhs=h_H(A Y A^T)`;
- positive margin;
- support-function maximizing matching and signs;
- standalone verifier status.

### `FAILED` or `ABSENT`

`FAILED` means an attempted run violated a gate or solver/certificate invariant. `ABSENT` means not run. Neither counts as success.

## Analytic unit tests

### Common mode

For `H=K_m`, unit edge costs,

\[
u=\mathbf1/\sqrt m,\qquad A=u,\qquad G_{\rm req}=[\gamma],
\]

verify

\[
\operatorname{OPT}_m(\gamma)=
\begin{cases}
0, & \gamma\le1,\\
\frac{m(\gamma-1)}2,
&1<\gamma\le1+\frac{2\lfloor m/2\rfloor}{m},\\
+\infty,&\gamma>1+\frac{2\lfloor m/2\rfloor}{m}.
\end{cases}
\]

Pinned examples:

| `m` | `gamma` | Expected disposition | Expected cost |
|---:|---:|---|---:|
| 6 | 1.0 | feasible, product control | 0.0 |
| 6 | 1.5 | feasible, nontrivial | 1.5 |
| 6 | 2.0 | feasible, pair-width maximum | 3.0 |
| 6 | 2.05 | infeasible | n/a |

### Multidirectional balance

Use a downstream map whose columns contain:

- normalized common mode `u`;
- at least one normalized contrast vector orthogonal to `u`;
- optionally a third localized or module-difference direction.

The registered `G_req` must make the product program infeasible and force allocation across more than one edge family. Record the contract margin per downstream eigenvector, not only the minimum eigenvalue.

### K3 exact controls

Retain all three fixed matrices:

1. pair-feasible signed target `(0.4,-0.3,0.2)` with cost `0.9`;
2. PSD but full-sector-infeasible all-negative `-0.4` target;
3. full-sector-feasible but pair-infeasible all-positive `+0.5` target.

These remain secondary controls; do not use them as the main experimental story.

## Correctness gates

All statuses are `PASS`, `FAIL`, or `ABSENT`; `ABSENT` never counts as `PASS`.

| Gate | Frozen requirement |
|---|---|
| C1 | Signed Bell/cat QFIM relative error `<1e-12` |
| C2a | Retained-label additivity below `1e-12` |
| C2b | Full coherent-flag covariance formula, including nonzero branch means |
| C2c | Dephasing removes the coherent covariance contribution |
| C3 | Exact matching-body equivalence across inequalities, columns, enumeration, and circuits |
| C4 | Three K3 sector-separation canaries behave exactly as specified |
| C5 | Exact-target reconstruction error `<1e-10` |
| C6 | Branch width, support bound when requested, and parsed expected-cost identity pass |
| C7 | Clifford-frame QFIM integrity passes; logical and physical width remain separate |
| C8 | Direct-backend hardware legality has zero illegal entanglers |
| C9 | All resource claims are parsed from emitted circuits |
| C10 | Compiled readout, downstream pullback, support/kernel, and inverse-stability controls pass |
| C11 | `F,t` primal feasibility and emitted schedule certificate pass |
| C12 | Formulation equivalence, dual feasibility, pricing, and relative gap pass |
| C13 | Full analytic common-mode curve is reproduced, including infeasible points |
| C14 | Every infeasible floor has a PSD witness, and every primary noncontrol floor defeats `F=I` |

A failed gate must write a failure artifact, invalidate all dependent results, preserve the frozen configuration, and stop the dependent campaign.

## Primary benchmark families

### IF1 - Common mode

This is the first scientific figure and the strongest analytic canary. Plot:

- requested `gamma`;
- exact cost curve;
- solver objective;
- product baseline;
- pair-width ceiling;
- dual lower bound;
- number of emitted branches.

### IF2 - Multidirectional contracts

Use downstream dimension two and three. Contracts should require the optimizer to balance common, contrast, and localized directions. Report the full downstream surplus spectrum and active edge pattern.

### IF3 - Hidden feasible contracts

Generate a feasible source schedule, choose `A`, and set

```text
G_req = A^T F_source A - S
```

for a registered PSD slack `S`. Hide the source schedule from the compiler. The task is not to recover the source schedule; it is to find the least-cost schedule satisfying the same floor.

### IF4 - Boundary and infeasible floors

Include:

- support-tangent feasible floors;
- floors just inside/outside a support hyperplane;
- hardware-nonedge requirements;
- odd-set bottlenecks;
- cost-induced alternative optima.

Every infeasible case needs a PSD support witness.

### IF5 - Structured multidirectional floors

Use block, banded, clustered, signed-community, graph-local, and low-effective-rank maps/contracts. Product-trivial cases must be filtered before execution.

### IF6 - Heterogeneous native costs

Hold the contract fixed while changing costs:

- unit active-pair count;
- calibrated native two-qubit duration;
- calibrated error proxy;
- intermodule penalty.

This tests whether contract freedom changes the emitted program in an operationally meaningful way.

### IF7 - Paper 1 application

Set

```text
F_Pi = W_rho
A = J
constraint = J^T F_Pi J >= G_theta_req
```

Declare the identifiable support of `J`. Never claim that changing `F_Pi` lifts `ker(J)`. Different feasible link-level QFIMs meeting the same metric-level contract are precisely the compiler freedom being evaluated.

## Baselines

### Primary information-floor baselines

1. product program `F=I`;
2. fixed maximum-matching mixture;
3. fixed hardware-local Bell policy;
4. variational direct contract-loss optimizer;
5. relaxed information-design lower bound that ignores QFIM achievability, clearly labeled as a lower bound;
6. a registered two-stage `choose exact target -> compile target` pipeline.

### Secondary exact-target baselines

1. QUEST or a faithful generic expectation-targeting implementation;
2. sparse Caratheodory pure-state preparation;
3. the same decomposition as a retained-label schedule;
4. coherent sparse data-state preparation;
5. explicit flag-data purification;
6. variational moment matching;
7. fixed global cat;
8. fixed block cats.

QUEST is a baseline for exact first/second-moment targeting, not an information-floor optimizer. Do not force it into the primary role without a separately specified outer optimization over targets.

## Experiment order

1. **Phase 0: theorem and semantics** - C1-C10 and fixed canaries.
2. **Phase 1: primary compiler** - 2E, C11-C13, common-mode curve, hidden feasible floors.
3. **Phase 2: certificates and multidirectional contracts** - 2F, C14, balance and infeasibility campaigns.
4. **Phase 3: secondary exact mode and baselines** - exact targets, projections, full-sector controls, QUEST/state-preparation comparisons.
5. **Phase 4: hardware and noise** - direct native first; routing and Clifford frames separately labeled.
6. **Phase 5: Paper 1 application** - `J^T F J` contract on a declared identifiable quotient.

No hardware/noise headline may be produced before the primary ideal compiler is certified.

## Canonical artifacts

```text
results/correctness_gates.json
results/information_floor_primary.json
results/common_mode.json
results/multidirectional.json
results/infeasible_floors.json
results/exact_target_controls.json
results/baselines.json
results/hardware.json
results/noise.json
results/application.json
results/traces/
results/programs/
results/circuits/
results/certificates/
```

Every information-floor record must include:

- schema version;
- source commit and dirty-state flag;
- config SHA-256;
- environment and solver versions;
- deterministic seed/instance ID;
- graph, `A`, `G_req`, and cost hashes;
- product-contract slack;
- primal `F,t` and all residuals;
- downstream surplus eigenvalues;
- matching decomposition and branch probabilities;
- emitted circuit and readout hashes;
- parsed resources and expected cost;
- dual `Y,beta`, lower bound, pricing residual, and gap;
- disposition and certificate type;
- standalone certificate-verifier result.

## Registration and claim firewall

- Correctness tolerances are frozen in `configs/experiment_v0.3.yaml`.
- Pilot runs may estimate runtime, memory, and effect sizes only.
- Main dimensions, graph families, seeds, thresholds, and primary metrics are frozen after the pilot and before the main run.
- Failed, absent, corrective, and post-hoc arms remain visible.
- A solver status alone is not a certificate.
- A feasible program without a verified dual bound is not called least cost.
- An infeasibility status without a PSD support witness is not called certified infeasible.
- A contract already satisfied by `F=I` is not a primary positive result unless explicitly registered as a product control.
- A formula-derived resource count is not an empirical hardware count until reconciled with emitted circuits.
- Noisy covariance is not mixed-state QFI.
- Routing does not silently enlarge the direct hardware graph.
- Exact target matching is not restored as the paper's sole motivation.

## Immediate merge target

The first repository merge under v0.3 should contain:

1. the v0.3 config and result schema;
2. proof/check artifacts for 2A-2D;
3. C1-C10;
4. the three K3 controls and coherent-flag loophole;
5. skeleton implementations for `contracts.py`, `floor_sdp.py`, `branch_master.py`, `dual.py`, and `certificates.py` with tests marked `ABSENT`, not falsely passing.

The second merge should close 2E, C11-C13, and `common_mode.json`.

The third merge should close 2F, C14, `multidirectional.json`, and `infeasible_floors.json`.

Only after those merges should the larger hardware, noise, QUEST, and Paper 1 campaigns begin.
