"""Apply the QUEST-variant coverage correction to the CovQ v0.10 PDF.

Usage:
    python patch_quest_variant_coverage.py \
        CovQ_Protected_Manuscript_v0.10.pdf \
        CovQ_Protected_Manuscript_v0.11.pdf
"""
from __future__ import annotations

import sys
import fitz

REG = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"
BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf"
WHITE = (1, 1, 1)
BLACK = (0, 0, 0)


def redact_many(page: fitz.Page, rects: list[fitz.Rect]) -> None:
    for rect in rects:
        page.add_redact_annot(rect, fill=WHITE)
    page.apply_redactions(
        images=fitz.PDF_REDACT_IMAGE_NONE,
        graphics=fitz.PDF_REDACT_LINE_ART_NONE,
        text=fitz.PDF_REDACT_TEXT_REMOVE,
    )


def put(page: fitz.Page, rect: fitz.Rect, text: str, size: float,
        *, align: int = 3, lineheight: float = 1.14,
        fontfile: str = REG) -> None:
    rem = page.insert_textbox(
        rect,
        text,
        fontname="djvbold" if fontfile == BOLD else "djvreg",
        fontfile=fontfile,
        fontsize=size,
        color=BLACK,
        align=align,
        lineheight=lineheight,
        overlay=True,
    )
    if rem < -0.01:
        raise RuntimeError(f"Page {page.number + 1}: replacement did not fit ({rem:.3f})")


def main(src: str, out: str) -> None:
    doc = fitz.open(src)

    page = doc[0]
    redact_many(page, [fitz.Rect(68, 488, 542, 614)])
    put(page, fitz.Rect(72, 492, 539, 611),
        "Across the eight exact targets completed by both published QUEST variants, tE and bE agreed "
        "exactly and QUEST emitted 2.50-8.89 times as many CX gates as CovQ's expected count. "
        "Terminal-insertion tE completed three additional deep targets with ratios 17.78, 25.56, and "
        "32.65, but best-position bE exceeded the registered compute budget there; those values are "
        "therefore tE comparisons and upper bounds on the separation from potentially shorter bE "
        "circuits. CovQ retains two-qubit depth one and pays additional settings. Under one depolarizing "
        "event per emitted CX, QUEST's optimistic QFI ceiling is more exposure-efficient at the tested "
        "0% and 1% points, whereas CovQ's attained deployable CFI is more efficient at every tested "
        "point from 2% onward and by 3.20 times at 20%. The narrowest crossover margin is more than "
        "512 times the worst disagreement between deterministic binomial enumeration and an independent "
        "Monte Carlo evaluation. A fixed-setting greedy extension is 7.24% suboptimal. The result is a "
        "certified resource phase diagram, not a claim of universal entanglement superiority.",
        7.85, lineheight=1.17)

    page = doc[34]
    redact_many(page, [fitz.Rect(60, 82, 552, 211)])
    put(page, fitz.Rect(63, 85, 549, 208),
        "We implement the published QUEST-tE and QUEST-bE procedures [1]. Each iteration inserts one "
        "Pauli rotation, then jointly reoptimizes all accumulated angles by L-BFGS. tE inserts "
        "terminally; bE searches insertion position, and the adjoint gradient is pinned against "
        "finite differences. bE completed 8 of 11 targets and agreed exactly with tE on each. For "
        "path-m4, path-m5, and banded-m5, bE exceeded the registered compute budget; Table 3 reports "
        "tE values and marks those rows. We never take a per-instance minimum across variants. The "
        "start state is |+>^m and is included in the emitted circuit. QUEST is used only for exact "
        "moment targeting. The earlier terminal-greedy routine remains a named ablation.",
        8.75, lineheight=1.16)

    page = doc[36]
    status_words = [w for w in page.get_text("words")
                    if w[4] == "converged" and 515 < w[1] < 690]
    if len(status_words) != 11:
        raise RuntimeError(f"Expected 11 Table 3 status cells; found {len(status_words)}")
    rects = [fitz.Rect(59, 430, 553, 487), fitz.Rect(60, 697, 552, 737)]
    rects.extend(fitz.Rect(w[0] - 2, w[1] - 1, w[2] + 2, w[3] + 1)
                 for w in status_words)
    redact_many(page, rects)
    put(page, fitz.Rect(63, 434, 550, 484),
        "Table 3: Exact-target realization using the independently implemented published QUEST "
        "procedures. 'Rot.' is the total Pauli-rotation ansatz length; '2q rot.' counts logical "
        "weight-two Pauli rotations. CX and CX depth are parsed after common lowering. The ratio is "
        "QUEST-tE emitted CX divided by CovQ expected CX. QUEST-bE completed 8 of 11 cases and agreed "
        "exactly with tE there; * marks tE-only rows.",
        7.45, lineheight=1.14)
    targets = ["matching-m3", "path-m3", "toeplitz-m3", "matching-m4", "path-m4",
               "toeplitz-m4", "matching-m5", "path-m5", "toeplitz-m5", "star-m4",
               "banded-m5"]
    te_only = {"path-m4", "path-m5", "banded-m5"}
    for target, word in zip(targets, sorted(status_words, key=lambda w: w[1])):
        put(page, fitz.Rect(484, word[1] - 1.5, 545, word[3] + 1.5),
            "tE only*" if target in te_only else "tE=bE",
            6.75, align=1, lineheight=1.0)
    put(page, fitz.Rect(63, 701, 550, 735),
        "QUEST-bE completed 8 of 11 instances and agreed exactly with tE on every completed case; "
        "the fully cross-validated emitted-CX ratio range is 2.50-8.89. Best-position insertion "
        "exceeded the registered compute budget on the three deepest targets, marked *.",
        8.05, lineheight=1.13)

    page = doc[37]
    redact_many(page, [fitz.Rect(59, 57, 553, 141)])
    put(page, fitz.Rect(62, 60, 550, 138),
        "QUEST-tE ratios on the three unconfirmed rows are 17.78, 25.56, and 32.65. Because these "
        "are the deepest cases - precisely where bE has the most opportunity to shorten the circuit - "
        "32.65 is a tE comparison and an upper bound on separation from bE, not a fully cross-validated "
        "maximum. Matching and Toeplitz targets show modest advantages and expose the settings tradeoff. "
        "The unconditional structural difference is depth: every nonempty CovQ branch has pair-preparation "
        "depth one; registered QUEST circuits have depth 2-36.",
        9.0, lineheight=1.11)

    page = doc[38]
    for x, y in [(317.8, 441.2), (351.8, 429.1), (385.7, 421.0)]:
        page.insert_text((x, y), "*", fontname="djvbold", fontfile=BOLD,
                         fontsize=5.8, color=BLACK, overlay=True)
    redact_many(page, [fitz.Rect(59, 633, 553, 673)])
    put(page, fitz.Rect(62, 636, 550, 671),
        "ceiling. (d) The noise-aware compiler uses three entangled settings through 0.15 and switches "
        "to the pairless column at 0.20. (e) Exact-target CX ratios: 8 bE-completed cases agree with "
        "tE and span 2.50-8.89; starred tE-only bars extend to 32.65 and upper-bound separation from "
        "bE. (f) Four Pauli rotations, three weight-two rotations, and six emitted CX are compared "
        "with at most two CX in a CovQ branch.",
        7.45, lineheight=1.08)

    page = doc[42]
    redact_many(page, [fitz.Rect(60, 477, 552, 560)])
    put(page, fitz.Rect(63, 480, 550, 557),
        "QUEST receives noisy QFI as an optimistic upper bound because no attaining noisy multiparameter "
        "readout was compiled; CovQ uses an emitted readout. Exact-target variant coverage is uneven: bE "
        "completed 8 of 11 targets and matched tE exactly there, so 2.50-8.89 is fully cross-validated. On "
        "path-m4, path-m5, and banded-m5, bE exceeded the registered budget; tE ratios 17.78, 25.56, and "
        "32.65 are upper bounds on separation from potentially shorter bE circuits. The generic sparse "
        "circuit remains a 232-CX synthesis upper bound.",
        9.0, lineheight=1.11)

    page = doc[43]
    redact_many(page, [fitz.Rect(60, 652, 552, 731)])
    put(page, fitz.Rect(63, 655, 550, 728),
        "Block-local noise preserves matching combinatorics while changing edge values and the optimal "
        "program. In the registered contract, the compiler moved from three entangled settings to one "
        "product setting as noise increased. The operational comparison is a phase diagram, not universal "
        "dominance: QUEST's optimistic QFI ceiling is lower at 0% and 1%, while CovQ's attained deployable "
        "CFI is lower-cost at every tested point from 2% onward.",
        9.0, lineheight=1.11)

    page = doc[44]
    redact_many(page, [fitz.Rect(59, 57, 553, 198)])
    put(page, fitz.Rect(62, 60, 550, 195),
        "In exact-target resources, bE completed 8 of 11 targets and matched tE exactly, giving a fully "
        "cross-validated ratio range of 2.50-8.89. On path-m4, path-m5, and banded-m5, bE exceeded the "
        "budget; the tE ratios 17.78, 25.56, and 32.65 are therefore upper bounds on separation from bE. "
        "The noise crossing is instance- and model-specific: the 2% margin is 3.09%, the conservative "
        "two-channel margin-to-uncertainty ratio is 512x, and the 20% exposure ratio is 3.20. Extra settings "
        "have an explicit amortization boundary, and the fixed-setting greedy method remains a negative result.\n\n"
        "Together, these results identify a tractable, certifiable regime for information-first quantum "
        "compilation. Exact pair-width geometry, emitted readout, finite-campaign guarantees, compiler "
        "scaling, and noise-aware cost support an end-to-end compiler claim without universal entanglement "
        "superiority. Clifford frames, coherent flags, correlated noise, routing, and real-device execution "
        "remain follow-on directions.",
        9.0, lineheight=1.10)

    meta = doc.metadata
    meta.update({
        "title": "CovQ: Certified Quantum Experiment Design from Fisher-Information Requirements",
        "author": "Hina Dixit; Abhinav Chauhan",
        "subject": ("Independent Research, Saratoga, CA. Certified quantum experiment design with "
                    "variant-qualified QUEST evidence, compiler scaling, finite-campaign guarantees, "
                    "and operational resource phase diagrams."),
        "keywords": ("quantum experiment design; Fisher information; quantum Fisher information; "
                     "matching polytope; Bell schedules; compiler scaling; matrix Bernstein; "
                     "distributed quantum sensing; QUEST-tE; QUEST-bE"),
    })
    doc.set_metadata(meta)
    doc.save(out, garbage=4, deflate=True, clean=True)
    doc.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: patch_quest_variant_coverage.py INPUT.pdf OUTPUT.pdf")
    main(sys.argv[1], sys.argv[2])
