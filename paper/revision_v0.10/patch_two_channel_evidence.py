import fitz
import shutil

SRC = "/mnt/data/CovQ_Protected_Manuscript_v0.10.pdf"
TMP = "/mnt/data/CovQ_Protected_Manuscript_v0.10.tmp.pdf"
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"

TEXT = (
    "At per-CX depolarization q in {0, 0.01}, the optimistic QUEST QFI ceiling "
    "requires less exposure than CovQ's deployable emitted-readout CFI. This does "
    "not establish an operational QUEST win because no attaining noisy QUEST readout "
    "is compiled. At every tested point from q = 0.02 onward, CovQ's attained CFI "
    "beats the QUEST QFI ceiling; the first registered crossover therefore lies in "
    "(0.01, 0.02]. At q = 0.02, the ceiling-to-deployable ratio is 1.031. The "
    "deployable template was evaluated by exact enumeration of the full binomial "
    "pilot support and independently by Monte Carlo pilot draws. Their maximum "
    "relative disagreement over the registered grid is 6.0e-5, so the narrowest "
    "3.09% crossing margin is more than 512 times larger than this two-channel "
    "uncertainty. At q = 0.20, CovQ requires 2.799 versus 8.961, a ratio of 3.201. "
    "A CovQ exposure below the QUEST QFI ceiling is operationally decisive because "
    "an actual QUEST measurement cannot exceed that ceiling. The mechanism is "
    "transparent: QUEST emits six CX gates, while a CovQ branch activates at most "
    "two disjoint Bell-pair CX gates at depth one. The frozen pilot increases CovQ "
    "exposure by only 1.1--1.4% relative to the matched-analyzer diagnostic."
)

doc = fitz.open(SRC)
page = doc[38]
rect = fitz.Rect(58, 58, 555, 248)
page.add_redact_annot(rect, fill=(1, 1, 1))
page.apply_redactions(
    images=fitz.PDF_REDACT_IMAGE_NONE,
    graphics=fitz.PDF_REDACT_LINE_ART_NONE,
    text=fitz.PDF_REDACT_TEXT_REMOVE,
)
remaining = page.insert_textbox(
    rect + (-1, 3, -1, -1),
    TEXT,
    fontname="djvreg",
    fontfile=FONT,
    fontsize=9.05,
    color=(0, 0, 0),
    align=3,
    lineheight=1.22,
    overlay=True,
)
if remaining < 0:
    raise RuntimeError(f"replacement did not fit: {remaining}")
doc.save(TMP, garbage=4, deflate=True, clean=True)
doc.close()
shutil.move(TMP, SRC)
print(f"patched {SRC}; remaining vertical space = {remaining:.3f}")
