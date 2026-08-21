# CovQ manuscript v0.11 - QUEST variant-coverage correction

This revision adds the missing evidentiary distinction between the QUEST terminal-insertion variant (tE) and best-position variant (bE). It does not change any CovQ theorem, compiler result, noise result, or numerical record.

## What changed

- **Abstract:** replaces the unqualified `2.50-32.65x` QUEST statement with two confidence tiers:
  - `2.50-8.89x` on the eight instances completed by both tE and bE, which agree exactly;
  - `17.78x`, `25.56x`, and `32.65x` on `path-m4`, `path-m5`, and `banded-m5`, where bE exceeded the registered compute budget.
- **Methodology (Section 13.4):** records bE coverage as `8/11`, states that no per-instance minimum is taken, and names the three tE-only rows.
- **Table 3:** marks every row as either `tE=bE` or `tE only*`; the caption and continuation explain that the tE-only values upper-bound the separation from potentially shorter bE circuits.
- **Figure 3(e):** adds stars to the three tE-only bars and explains the coverage distinction in the caption.
- **Limitations (Section 16.4):** explicitly treats `32.65x` as a tE comparison and an upper bound on the separation from bE, not a fully cross-validated maximum.
- **Conclusion:** carries the same distinction into the final claims.

## Final claim hierarchy

1. Fully cross-validated exact-target ratio range: **2.50-8.89x**.
2. Additional tE-only deep-target ratios: **17.78x, 25.56x, 32.65x**.
3. The tE-only values are measured comparisons against tE and upper bounds on the separation from bE.
4. The unconditional structural result remains two-qubit preparation depth one for every nonempty CovQ branch.

## Unchanged operational result

The deployable noise crossover remains supported by deterministic binomial enumeration and Monte Carlo pilot evaluation:

- CovQ beats the QUEST QFI ceiling at every tested point from `2%` onward;
- narrowest margin: `3.09%`;
- conservative margin-to-uncertainty ratio: `512x`;
- exposure ratio at `20%`: `3.20x`.

## Delivery

The rendered PDF remains 51 pages and retains its existing metadata, bookmarks, links, and page numbering. The file hash is recorded in `DELIVERY_SHA256.txt`.
