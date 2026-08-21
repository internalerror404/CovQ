import fitz
import shutil

SRC = "/mnt/data/CovQ_Protected_Manuscript_v0.10.pdf"
TMP = "/mnt/data/CovQ_Protected_Manuscript_v0.10.tmp.pdf"
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf"
TEXT = (
    '[34] Timothy J. Proctor, Paul A. Knott, and Jacob A. Dunningham. "Multiparameter Estimation in Networked Quantum Sensors." Physical Review Letters 120, 080501 (2018). doi:10.1103/PhysRevLett.120.080501.\n'
    '[35] Zachary Eldredge, Michael Foss-Feig, Jonathan A. Gross, S. L. Rolston, and Alexey V. Gorshkov. "Optimal and Secure Measurement Protocols for Quantum Sensor Networks." Physical Review A 97, 042337 (2018). doi:10.1103/PhysRevA.97.042337.\n'
    '[36] Jesus Rubio and Jacob Dunningham. "Bayesian Multiparameter Quantum Metrology with Limited Data." Physical Review A 101, 032114 (2020). doi:10.1103/PhysRevA.101.032114.\n'
    '[37] Nathan Shettell and Damian Markham. "Graph States as a Resource for Quantum Metrology." Physical Review Letters 124, 110502 (2020). doi:10.1103/PhysRevLett.124.110502.\n'
    '[38] Joel A. Tropp. "User-Friendly Tail Bounds for Sums of Random Matrices." Foundations of Computational Mathematics 12, 389-434 (2012). doi:10.1007/s10208-011-9099-z.'
)

doc = fitz.open(SRC)
page = doc[50]
page.add_redact_annot(fitz.Rect(60, 584.5, 554, 715), fill=(1, 1, 1))
page.apply_redactions(
    images=fitz.PDF_REDACT_IMAGE_NONE,
    graphics=fitz.PDF_REDACT_LINE_ART_NONE,
    text=fitz.PDF_REDACT_TEXT_REMOVE,
)
remaining = page.insert_textbox(
    fitz.Rect(64, 589, 550, 710),
    TEXT,
    fontname="djvreg",
    fontfile=FONT,
    fontsize=7.1,
    color=(0, 0, 0),
    align=0,
    lineheight=1.13,
    overlay=True,
)
if remaining < 0:
    raise RuntimeError(f"references did not fit: {remaining}")
doc.save(TMP, garbage=4, deflate=True, clean=True)
doc.close()
shutil.move(TMP, SRC)
print(f"patched {SRC}; remaining vertical space = {remaining:.3f}")
