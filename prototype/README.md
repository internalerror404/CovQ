# prototype — the CovQ reference implementation

Supplies every number quoted in the manuscript and in `../docs/audits/`. Pure
Python + NumPy + SciPy; no quantum SDK, no network, no GPU. `SCIENTIFIC_SCOPE = FROZEN`
(see `../REGISTRATION.md`).

## Running it

```bash
python3 -m pip install -r requirements.txt
PYTHONPATH=src python3 -m pytest tests -q
PYTHONPATH=src python3 scripts/run_evidence.py --out results/prototype_evidence.json
PYTHONPATH=src python3 scripts/emit_measurement_records.py
PYTHONPATH=src python3 scripts/emit_k3_baseline.py
PYTHONPATH=src python3 scripts/emit_noise_records.py
PYTHONPATH=src python3 scripts/emit_campaign_records.py
```

Records land in `../results/` carrying source commit, dirty flag, tolerances, content
hashes, and a `status` of `DONE` or `FAILS`. **No manuscript number is transcribed by
hand.**

## Modules

| module | contents |
|---|---|
| `paulis` | commuting Pauli families, simultaneous diagonalisation |
| `circuits`, `sim` | circuit IR, lowering to CX, routing, statevector simulation |
| `qfim` | pure-state QFIM, identifiability, CRB conditioning |
| `polytope`, `width` | `Q_m`, `MATCH(H)`, Edmonds separation, signed schedules |
| `programs` | labelled schedules, coherent flags, unlabelled mixtures |
| `floor` | the information-floor compiler, conic dual, Farkas certificates |
| `measurement` | block-local readout, mixed-state QFIM, SLD solve, pinching |
| `estimator` | parity likelihood, branch-conditioned MLE, adaptive recentering |
| `noise` | block-local noise, noisy pricing oracle, Eq (111), deployable policy |
| `quest` | the QUEST baseline (arXiv:2605.02367), reimplemented from its method |
| `baselines`, `instances`, `gates` | comparators, target families, correctness gates |

## What the tests are for

They are written so that a **false theorem** fails them, not merely a broken
implementation. Where a claim is checkable against brute force it is: Edmonds separation
against full enumeration, the noisy pricing oracle against exhaustive signed matchings,
the spectral mixed-state QFIM against an independent Lyapunov solve, greedy support
selection against exhaustive enumeration.

Two tests exist to record **failures**, and must keep failing what they fail:

- `test_n9_greedy_support_selection_is_not_exact` — greedy is 7.24 % over optimum in the
  multi-setting regime.
- the zero-noise row of the N10 campaign, where QUEST is the more shot-efficient arm.

## Known limitations, stated rather than discovered

- `shot_scaled_edge_compile` stalls once entanglement activates; the branch formulation is
  the working path. It needs an interior-point SDP backend.
- Sparse-state preparation is a generic dense routine (232 CX at `m = 4`). Every count it
  produces is an **upper bound**, so no per-shot gate-count advantage over it is admissible.
- Noise coverage is block-local only: dephasing, depolarizing, two-qubit gate noise, idle,
  readout confusion. Crosstalk, correlated branch noise, route collisions and coherent
  inter-block errors are out of scope by construction.
- The QUEST arm has no compiled readout, so its column is a QFI **upper bound**, not an
  attainable number.
