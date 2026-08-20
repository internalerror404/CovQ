# CovQ

Pair-entangled compilation of commuting-Pauli information requirements.

**Read first:** [`STATUS.md`](STATUS.md) — what is background, what is paper-grade, what is
superseded, and what is not claimed.

## Layout

```
STATUS.md                                       disposition of every item
docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md
docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md
docs/audits/PRIOR_ART_SWEEP_TASK0B5_v0.1.md
docs/proofs/THEOREM_PAIR_WIDTH_AND_INFORMATION_FLOOR_v0.1.md   2A-2F
charter/CovQ_Research_Charter_v0.3.md           current charter
charter/CovQ_Experiment_Protocol_v0.3.yaml      current protocol
charter/CovQ_Paper_Skeleton_v0.3.tex            current skeleton
archive/v0.1/, archive/v0.2/                    superseded documents, unmodified
prototype/                                      PARTIAL; the audit's numerical evidence
```

## Where the project stands

v0.1 pivoted. The declared v0.1 target sector — zero-mean, unit-diagonal — is the sector in
which four independent degeneracies collapse the distinctive machinery, and its headline
theorems are either inherited (Pitowsky, Huber–Marić), already proved elsewhere
(Caprara et al.), or one line of submultiplicativity. The specification is also a special
case of QUEST's input language, so the "we specify something weaker than a state" framing is
dead.

What CovQ compiles is a **Fisher-information contract**: given `A` naming the parameter
combinations that matter, a required floor `G_req`, and a hardware graph `H`, return the
least-cost program with `AᵀF_ΠA ⪰ G_req` — or a certificate that none exists at the requested
width.

The enabling theorem makes that tractable:

> **A labelled schedule of entanglement width two on `H` realises `F` exactly when
> `diag F = 1`, `F` is supported on `E(H)`, and the matrix of off-diagonal magnitudes lies in
> the matching polytope of `H`.**

Constructive — it emits parallel Bell pairs, every branch of two-qubit depth exactly 1 —
polynomial-time decidable by Edmonds separation against NP-hard membership at unbounded width,
and with blossom certificates when a contract provably needs more than pair entanglement. It
turns the contract into a convex program whose objective *is* the expected per-shot Bell-pair
bill.

Two consequences are analytic and make good unit tests:
`cost*(γ) = m(γ−1)/2` for a common-mode contract, and a collective mode on `k` qubits capped at
`2` for even `k` but `2 − 1/k` for odd `k` — Edmonds' blossom inequality restated as physics.

## Next

Steps 1, 2A–2F and 5 of Task 0B.5 are done. Remaining: the width-3 attribution check (step 3),
and closing the prototype gaps (noise, a competitive sparse-preparation routine, QUEST).

The width-3 hardness is expected to be a short corollary of Partition into Triangles rather
than a theorem; its role is to explain why pair width is a privileged tractable regime.

## Reproducing the prototype numbers

```bash
cd prototype
python3 -m pip install -r requirements.txt      # numpy, scipy, pytest
PYTHONPATH=src python3 -m pytest tests -q
PYTHONPATH=src python3 scripts/run_evidence.py --out results/prototype_evidence.json
```

Runs on a stock macOS Python 3.11+ with numpy and scipy. No Qiskit, no GPU, no network.
