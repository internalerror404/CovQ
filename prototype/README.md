# prototype — PARTIAL

**This is not Task 0C v0.1.** Task 0C v0.1 is `SUPERSEDED` (see `../STATUS.md`). This tree
exists for one reason: to supply the numerical evidence quoted in
`../docs/audits/TASK0AB_NOVELTY_THEOREM_AUDIT_v0.2.md`.

It was written before the pivot decision arrived, against the v0.1 protocol. It is retained
rather than discarded because it produced three results that changed the documents:

1. **The `k = 2` theorem checks out against brute force.** Edmonds separation versus full
   enumeration of the extreme points of `Q^prog_{m,2}`, on 160 random unit-diagonal symmetric
   matrices at `m ∈ {3,4,5,6}` deliberately including infeasible ones: zero disagreements.
2. **The blossom facets are active, not redundant.** For a triangle carrying common
   off-diagonal `c`, degree constraints permit `c ≤ 1/2` but the odd-set constraint permits
   only `c ≤ 1/3`, and both the Edmonds test and brute force switch at exactly `1/3` while
   the whole family stays inside `Q_m`. Without this the characterisation could have been a
   restatement of degree counting.
3. **It retracts a claim.** The first-pass review reported the "relaxed" coherent-flag
   construction as ~12× cheaper than the literal controlled-cat form. It is not a like-for-
   like comparison: the cheap form moves the information into flag–data coherence and its
   QFIM vanishes under flag dephasing. Gate `C2c` now measures exactly this.

It does **not** close step 6 of Task 0B.5. Missing against that list: noise stages beyond
the ideal statevector, the information-floor SDP, a competitive sparse-preparation routine,
and QUEST.

## Run

```bash
python3 -m pip install -r requirements.txt
PYTHONPATH=src python3 -m pytest tests -q                       # ~2 s
PYTHONPATH=src python3 scripts/run_evidence.py --out results/prototype_evidence.json
```

Stock macOS Python 3.11+, numpy and scipy only. No Qiskit, no GPU, no network. The gate suite
stops at the first failure and writes the report without touching any tolerance.

## Modules

| module | contents |
|---|---|
| `paulis.py` | symplectic Pauli algebra; simultaneous Clifford diagonalisation. Signs are read off *numerically* rather than propagated through a phase rule, which removes a whole class of silent sign bugs. |
| `qfim.py` | per-shot QFIM from statevectors; SLD saturability residual; conditioning and pseudoinverse stability. |
| `polytope.py` | `Q_m` vertex enumeration, exact LP, Carathéodory reduction, minimum-support MILP, column generation, Frank–Wolfe **with away steps**, and the hypermetric / uncertainty-relation certificates. |
| `width.py` | **the spine.** Matching-polytope membership and Edmonds separation for `Q^lab_{m,2}` and the hardware-native `Q^lab_{H,2}`; constructive decomposition by column generation over an exact subset-DP max-weight-matching oracle; brute-force width-`k` enumeration as the independent cross-check. |
| `circuits.py` | framework-neutral IR, lowering to `{1q, CX}`, routing, resource accounting, QASM export. |
| `sim.py` | dense statevector simulator; product-partition and finest-partition tests used by the width gates. |
| `programs.py` | signed cats, labelled schedules, three flagged constructions including the width loophole, and the flag resource report. |
| `baselines.py` | the realisation-matched baselines, with the deliberately-weak ones labelled as such. |
| `instances.py` | benchmark families. Ground truth is `True`, `False`, or **`None`** — the third case is deliberate, so the gates see instances whose answer was not baked in. |
| `gates.py` | C1–C12, with C2/C3/C10 re-specified so they can fail, and C10/C12 reported as MEASURED rather than PASS. |

## Conventions

- QFIMs are **per shot**, everywhere.
- Resource counts come from **emitted circuits** after lowering and routing, never a formula.
- Qubit ordering is little-endian.
- Caps are explicit and raise rather than silently degrade.

## Known limitations that bound what may be claimed

- The sparse state preparation is a generic prefix-tree construction with v-chain lowering,
  not a competitive routine. **Every single-state gate count it produces is an upper bound**,
  so no per-shot gate-count advantage over the single-state realisation is admissible from
  these numbers.
- `unlabelled_mixed_qfi` builds a dense density matrix and is capped at `m = 8`.
- Odd-set separation is by enumeration, capped at `m = 20`. The polynomial-time route
  (Padberg–Rao) is what the complexity claim rests on and is not implemented.
- Exact pricing for `Q_m` column generation enumerates `2^(m-1)` sign vectors. That is the
  known NP-hardness showing up, not an implementation shortcut.
- No noise model. No QUEST.
