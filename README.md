# CovQ

Pair-entangled compilation of commuting-Pauli information requirements.

**Read first:** [`STATUS.md`](STATUS.md) — what is background, what is paper-grade, what is
superseded, and what is not claimed.

## Layout

```
STATUS.md                                       disposition of every item
docs/decisions/DECISION_001_PIVOT_TO_BOUNDED_WIDTH.md
docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md
charter/CovQ_Research_Charter_v0.2.md           current charter
charter/CovQ_Experiment_Protocol_v0.2.yaml      current protocol
charter/CovQ_Paper_Skeleton_v0.2.tex            current skeleton
archive/v0.1/                                   the superseded v0.1 documents, unmodified
prototype/                                      PARTIAL; the audit's numerical evidence
```

## Where the project stands

v0.1 pivoted. The declared v0.1 target sector — zero-mean, unit-diagonal — is the sector in
which four independent degeneracies collapse the distinctive machinery, and its headline
theorems are either inherited (Pitowsky, Huber–Marić), already proved elsewhere
(Caprara et al.), or one line of submultiplicativity. The specification is also a special
case of QUEST's input language, so the "we specify something weaker than a state" framing is
dead.

What survives, and is now the spine:

> **A target QFIM is realisable by a labelled schedule of entanglement width two on a
> hardware graph `H` exactly when its matrix of off-diagonal magnitudes, supported on `E(H)`,
> lies in the matching polytope of `H`.**

Constructive (it emits parallel Bell pairs on physical edges), polynomial-time decidable by
Edmonds separation against NP-hard membership at unbounded width, and with explicit blossom
certificates when a target provably needs more than pair entanglement.

Compiler semantics move from exact matching to an **information floor**,
`min Cost(Π) s.t. F_Π ⪰ G⋆`, which un-degenerates Loewner dominance without leaving the
audited sector.

## Next

Task 0B.5, in the order given in `DECISION_001`: prior-art sweep → formal proof of the `k=2`
theorem for `K_m` and `H` → prove or abandon the width-2/width-3 boundary → write up the
coherent-flag resource semantics → rewrite charter/protocol/abstract *(drafted)* → then the
prototype.

The paper needs the width-3 hardness result. The `k=2` theorem clears the project's own K4
gate but does not by itself carry a major paper, because the combinatorics underneath it are
classical matching theory.

## Reproducing the prototype numbers

```bash
cd prototype
python3 -m pip install -r requirements.txt      # numpy, scipy, pytest
PYTHONPATH=src python3 -m pytest tests -q
PYTHONPATH=src python3 scripts/run_evidence.py --out results/prototype_evidence.json
```

Runs on a stock macOS Python 3.11+ with numpy and scipy. No Qiskit, no GPU, no network.
