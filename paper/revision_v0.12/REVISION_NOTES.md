# CovQ v0.12 revision notes

This pass addresses the genuine points in the two external reviews while leaving the remaining bE runs intentionally unexecuted.

## Scientific corrections

- States the universal pair-width limitation `0 <= F <= 2I` and the resulting twofold ideal-exposure ceiling relative to product probes.
- Distinguishes the fixed-exposure Bernstein proposition from the primary shot-scaled compiler and adds the shot-scaled branch-allocation corollary.
- States the flagship campaign's ledger explicitly: `c0 = 1`, `ce = 0`; Table 5 is an exposure comparison, not full execution cost.
- Adds cost curves over `ce/c0`, including the pair-price points at which exposure savings cease to be cost savings.
- Adds a deliberately multidirectional contract; its noiseless product/CovQ exposure ratio is 1.199 rather than 1.823 for the registered contract.
- Adds refined per-CX grid points and an independent amplitude-damping model. The crossing brackets are `(0.0125, 0.015]` for depolarization and `(0.015, 0.0175]` for amplitude damping on the registered instance.
- Reframes deterministic enumeration versus Monte Carlo as numerical implementation validation only, not physical-model uncertainty.
- Separates the specification-matched randomized-schedule comparator from the single-state QUEST comparison. The matched Pauli-rotation construction emits exactly two CX per active pair versus one for CovQ.
- Keeps the fully cross-validated QUEST range at `2.50-8.89x`; the three tE-only ratios remain marked as upper bounds on separation from bE. No new bE runs were made.

## Production corrections

- Rebuilt the abstract, contents, introduction, related work, Sections 9-17, and Appendix F in a consistent layout without floating tcolorboxes.
- Removed the duplicated full Section 8; Section 8 is now a pointer and the background appears once in Appendix G.
- Restored Section 5.8's definition, proof, and k-producibility attribution.
- Made setup and switching costs consistent between Eq. (63) and the campaign ledger.
- Repaired blank/misordered contents entries and severed sentences.
- Numbered Tables 1-8 uniquely; the notation table is unnumbered.
- Replaced the unavailable companion-paper reference with the public CovQ repository citation.
- Added the public repository URL and explicit DOI-pending language.
- Shortened the downstream-map application to one illustrative paragraph.

## Scope

The alternative noise model demonstrates model sensitivity within the block-local framework; it does not establish a hardware threshold. Real-device execution, correlated noise, routing, and globally optimal fixed-charge synthesis remain outside scope.

Repository: https://github.com/internalerror404/CovQ
