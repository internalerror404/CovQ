import fitz
import shutil

SRC = "/mnt/data/CovQ_Protected_Manuscript_v0.10.pdf"
TMP = "/mnt/data/CovQ_Protected_Manuscript_v0.10.tmp.pdf"
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"
TEXT = (
    "Distributed-sensing theory already separates local multiparameter estimation "
    "from global-function estimation: entanglement may offer little intrinsic gain "
    "for many independent local parameters, but it can enhance selected global "
    "linear functions; limited-data and graph-state analyses further qualify that "
    "advantage [34-37]. Equations (57)-(60) specialize that distinction to CovQ's "
    "pair-width cost model. The contribution here is not the sensing principle, "
    "but the exact achievable matrix body, minimum-cost hardware compilation, "
    "emitted readout, and certificates."
)

doc = fitz.open(SRC)
page = doc[21]
redaction = fitz.Rect(60, 116.5, 552, 164.5)
page.add_redact_annot(redaction, fill=(1, 1, 1))
page.apply_redactions(
    images=fitz.PDF_REDACT_IMAGE_NONE,
    graphics=fitz.PDF_REDACT_LINE_ART_NONE,
    text=fitz.PDF_REDACT_TEXT_REMOVE,
)
remaining = page.insert_textbox(
    fitz.Rect(63, 120, 549, 162.5),
    TEXT,
    fontname="djvreg",
    fontfile=FONT,
    fontsize=7.35,
    color=(0, 0, 0),
    align=3,
    lineheight=1.08,
    overlay=True,
)
if remaining < 0:
    raise RuntimeError(f"replacement did not fit: {remaining}")
doc.save(TMP, garbage=4, deflate=True, clean=True)
doc.close()
shutil.move(TMP, SRC)
print(f"patched {SRC}; remaining vertical space = {remaining:.3f}")
