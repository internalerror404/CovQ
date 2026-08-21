import fitz

SRC='/mnt/data/CovQ_Protected_Manuscript_v0.9.pdf'
REP='/mnt/data/covq_replacement_pages.pdf'
OUT='/mnt/data/CovQ_Protected_Manuscript_v0.10.pdf'
REG='/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf'
BOLD='/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf'
WHITE=(1,1,1); BLACK=(0,0,0); BLUE=(0.0,0.0,0.55)

src=fitz.open(SRC); rep=fitz.open(REP); out=fitz.open()
out.insert_pdf(src,from_page=0,to_page=22,links=True,annots=True)
out.insert_pdf(rep,from_page=0,to_page=0,links=True,annots=True)
out.insert_pdf(src,from_page=24,to_page=31,links=True,annots=True)
out.insert_pdf(rep,from_page=1,to_page=1,links=True,annots=True)
out.insert_pdf(src,from_page=33,to_page=39,links=True,annots=True)
out.insert_pdf(rep,from_page=2,to_page=2,links=True,annots=True)
out.insert_pdf(src,from_page=40,to_page=45,links=True,annots=True)
out.insert_pdf(rep,from_page=3,to_page=4,links=True,annots=True)
out.insert_pdf(src,from_page=48,to_page=49,links=True,annots=True)


def whiteout(page,rect):
    page.add_redact_annot(rect,fill=WHITE)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,graphics=fitz.PDF_REDACT_LINE_ART_NONE,text=fitz.PDF_REDACT_TEXT_REMOVE)


def box(page,rect,text,size=8.5,path=REG,color=BLACK,align=0,lineheight=1.2):
    return page.insert_textbox(rect,text,fontname=('djvbold' if path==BOLD else 'djvreg'),fontfile=path,fontsize=size,color=color,align=align,lineheight=lineheight,overlay=True)


def line(page,x,y,text,size=10,path=REG,color=BLACK):
    page.insert_text((x,y),text,fontname=('djvbold' if path==BOLD else 'djvreg'),fontfile=path,fontsize=size,color=color,overlay=True)


def footer(page,n):
    whiteout(page,fitz.Rect(280,738,332,775))
    box(page,fitz.Rect(280,742,332,765),str(n),size=8.5,path=REG,align=1)

p=out[0]
whiteout(p,fitz.Rect(55,276,558,681))
line(p,273,292,'Abstract',size=10.2,path=BOLD)
abstract=(
"Whether entanglement is worth using is a regime question, not a universal yes-or-no claim. CovQ turns a downstream Fisher-information requirement into a minimum-cost quantum experiment and certifies the resulting resource phase diagram across exposure, pair activation, settings, and noise. Given commuting Pauli generators G, a hardware graph H, a downstream map A, a floor G_req >= 0, and a declared cost model, it emits an exposure-weighted retained-label schedule, native circuits, a two-quadrature pilot, ancilla-free readout, an estimator, and a primal-dual certificate satisfying the requested downstream Fisher floor.\n\n"
"The exact pair-width regime is completely characterized. For physical Z_i/2 generators on H=(V,E), a unit-diagonal matrix F is achievable by width-two branches exactly when its nonedge entries vanish and (|F_ij|) on E lies in Edmonds' matching polytope. Matching decompositions become depth-one schedules of parallel signed Bell pairs, maximum-weight matching separates the dual circuit constraints, and every pair-width matrix obeys 0 <= F <= 2I. Thus pair entanglement can redistribute information but cannot concentrate more than two Fisher units into any eigenmode; width-three utility optimization is already NP-hard.\n\n"
"Product equatorial measurements attain every emitted signed-cat block, parity is sufficient, and a branch-conditioned maximum-likelihood estimator is empirically efficient. A two-quadrature pilot both phases the analyzer and selects the correct periodic likelihood chart. Under block-local noise, nonempty branches remain a common depth-one idling template plus signed edge increments, so pricing remains matching. A compiler-only scaling study closed all 30 registered m=50-500 grid and heavy-hex instances; at m=500, median times were 2.39 s and 1.67 s with 40 and 35 generated columns.\n\n"
"Across eleven exact targets, CovQ used 2.5-32.7 times fewer expected emitted CX gates than the published QUEST procedures while retaining two-qubit depth one and paying additional settings. Under one depolarizing event per emitted CX, QUEST's optimistic QFI ceiling is more exposure-efficient at the tested 0% and 1% points, whereas CovQ's attained deployable CFI is more efficient from 2% onward and by 3.20 times at 20%. A fixed-setting greedy extension is 7.24% suboptimal. The result is a certified resource phase diagram, not a claim of universal entanglement superiority."
)
box(p,fitz.Rect(72,302,540,677),abstract,size=8.15,path=REG,align=3,lineheight=1.20)

p=out[3]
whiteout(p,fitz.Rect(55,347,557,714))
toc4=(
"14  Experimental results                                                                  36\n"
"    14.1  Structural theorem and solver checks  . . . . . . . . . . . . . . . . . . . . . .   36\n"
"    14.2  Local readout, parity sufficiency, and estimator efficiency  . . . . . . . . . .   36\n"
"    14.3  Deployable pilot recentering  . . . . . . . . . . . . . . . . . . . . . . . . .   36\n"
"    14.4  Exact-target comparison with QUEST  . . . . . . . . . . . . . . . . . . . . .   37\n"
"    14.5  Block-local noise preserves matching compilation  . . . . . . . . . . . . .   38\n"
"    14.6  Equal-accounting noisy information-floor campaign  . . . . . . . . . . . .   38\n"
"    14.7  Fixed-setting greedy selection fails  . . . . . . . . . . . . . . . . . . . .   40\n"
"    14.8  Setting-cost amortization  . . . . . . . . . . . . . . . . . . . . . . . . . .   40\n"
"    14.9  Compiler-only scale on grid and heavy-hex graphs  . . . . . . . . . . . .   41\n"
"    14.10 Reproducibility status  . . . . . . . . . . . . . . . . . . . . . . . . . . .   42\n\n"
"15  Illustrative application interface: Lorentzian light-ray tomography                 42\n"
"    15.1  The application is exactly a Fisher-information contract  . . . . . . . . .   42\n"
"    15.2  The geometric kernel remains absolute  . . . . . . . . . . . . . . . . . .   42\n"
"    15.3  Status of this interface  . . . . . . . . . . . . . . . . . . . . . . . . . . .   42\n\n"
"16  Limitations and open problems                                                        43\n"
"    16.1  Exact scope of the reachable-body theorem  . . . . . . . . . . . . . . . .   43\n"
"    16.2  Program width is instruction-set relative  . . . . . . . . . . . . . . . . .   43\n"
"    16.3  Block-local noise is a real restriction  . . . . . . . . . . . . . . . . . .   43\n"
"    16.4  One-state baselines are not fully operational  . . . . . . . . . . . . . .   43\n"
"    16.5  Fixed setting charges remain nonconvex  . . . . . . . . . . . . . . . . .   43\n"
"    16.6  Adaptive readout and periodic identifiability  . . . . . . . . . . . . . .   43\n"
"    16.7  Finite schedules and finite samples  . . . . . . . . . . . . . . . . . . .   44\n"
"    16.8  Simulation and repository scale  . . . . . . . . . . . . . . . . . . . . .   44\n"
"    16.9  Novelty remains falsifiable  . . . . . . . . . . . . . . . . . . . . . . .   44"
)
box(p,fitz.Rect(63,352,550,708),toc4,size=7.75,path=REG,color=BLUE,lineheight=1.20)

p=out[4]
whiteout(p,fitz.Rect(55,55,557,278))
toc5=(
"17  Conclusion                                                                              44\n\n"
"A   Notation and implementation conventions                                                45\n\n"
"B   Proof details for the coherent-flag formula                                            45\n\n"
"C   Proof details for matching necessity                                                    46\n\n"
"D   A small exact example                                                                   46\n\n"
"E   Repository result schema and provenance                                                 46\n\n"
"F   Reproducibility and data availability  . . . . . . . . . . . . . . . . . . . . . . .   47\n\n"
"G   Full unit-diagonal feasibility as background                                             48\n"
"    G.1  The symmetric-Bernoulli correlation sector  . . . . . . . . . . . . . . . . .   48\n"
"    G.2  Signed cats and the baseline they do not defeat  . . . . . . . . . . . . . .   49\n"
"    G.3  A positive-semidefinite matrix outside the full sector  . . . . . . . . . .   49"
)
box(p,fitz.Rect(63,62,550,270),toc5,size=8.35,path=REG,color=BLUE,lineheight=1.31)

p=out[35]
whiteout(p,fitz.Rect(55,512,557,589))
bootstrap_results=("The branch-conditioned parity MLE was consistent on both operating points and all three shot counts. At N=50,000, the whitened empirical-efficiency eigenvalues across the two operating points ranged from 0.9525 to 1.0579; all biases remained within four reported standard errors. A 100,000-draw Gaussian parametric bootstrap gives a simultaneous 95% envelope [0.893, 1.117] across the two operating points, containing all six observed eigenvalues. The machine-readable bootstrap artifact records the seed, Wishart construction, per-case intervals, and joint envelope.")
box(p,fitz.Rect(63,517,550,586),bootstrap_results,size=8.45,path=REG,align=3,lineheight=1.18)

p=out[41]
whiteout(p,fitz.Rect(55,55,557,186))
line(p,63,78,'14.10   Reproducibility status',size=10.5,path=BOLD)
body=("All scientific components required for the paper's stated scope are present: exact pair geometry, primal-dual certificates, operational readout and estimator, adaptive recentering, block-local noisy pricing, the corrected published-QUEST comparison, settings amortization, the N9 negative result, and the compiler-only m=50-500 scale campaign. The symmetric-accounting correction affected only QUEST-dependent resource and noise quantities; no CovQ theorem, certificate, readout, estimator, or internal compiler result changed. Numerical tables and figures are generated from canonical records, with proof and citation audits used as internal quality-control checks.")
box(p,fitz.Rect(63,88,550,181),body,size=8.55,path=REG,align=3,lineheight=1.22)

p=out[46]
whiteout(p,fitz.Rect(55,132,557,292))
records=(
"The principal records are:\n"
"- results/measurements/ideal_block_readout.json;\n"
"- results/measurements/schedule_attainability.json;\n"
"- results/measurements/estimator_efficiency.json;\n"
"- artifacts/estimator_efficiency/bootstrap_parametric_v0.1.json;\n"
"- artifacts/compiler_scale_v0.1/compiler_scale.json;\n"
"- results/noise/adaptive_readout_performance.json;\n"
"- results/noise/noise_aware_floor_compiler.json;\n"
"- results/noise/noisy_pricing_oracle.json;\n"
"- results/noise/quest_operational_comparison.json;\n"
"- artifacts/quest_reference_reimplementation_v1/;\n"
"- results/noise/fixed_setting_cost_solver.json;\n"
"- results/noise/setting_cost_phase_diagram.json; and\n"
"- results/baselines_k3_quest.json.\n\n"
"Every paper table and figure is generated from these artifacts; manual transcription is not used for numerical results."
)
box(p,fitz.Rect(63,137,550,288),records,size=7.65,path=REG,lineheight=1.17)

for page_no in range(42,48):
    footer(out[page_no-1],page_no)
for page_no in (50,51):
    footer(out[page_no-1],page_no)

meta=out.metadata
meta.update({'title':'CovQ: Certified Quantum Experiment Design from Fisher-Information Requirements','author':'Hina Dixit; Abhinav Chauhan','subject':'Independent Research, Saratoga, CA. Certified quantum experiment design, operational resource phase diagram, compiler scaling, finite-campaign guarantee, and reproducible comparator artifacts.','keywords':'quantum experiment design; Fisher information; quantum Fisher information; matching polytope; Bell schedules; compiler scaling; matrix Bernstein; distributed quantum sensing; QUEST'})
out.set_metadata(meta)
toc=[
 [1,'Abstract',1],
 [1,'1 Introduction',6],[2,'1.1 A certified resource phase diagram from information requirements',6],[2,'1.6 Scope and validation status',8],
 [1,'2 Related work and novelty boundary',8],[1,'3 Fisher-information contracts and program semantics',9],
 [1,'4 Exact pair-width geometry: the compiler engine',13],[1,'5 Certified information-floor compilation',17],
 [2,'5.7 Two exact common-mode curves',21],[2,'5.8 Multidirectional contracts and witness validation',21],
 [1,'6 Operational attainability and measurement compilation',22],[2,'6.5 Adaptive recentering and periodic identifiability',24],[2,'6.6 Estimator efficiency on the identifiable quotient',24],
 [1,'7 A width-two/width-three complexity boundary',25],[1,'8 Full unit-diagonal feasibility (moved to Appendix G)',26],
 [1,'9 Compiler algorithms',26],[1,'10 Generator-aware and hardware-aware backends',28],
 [1,'11 Noise-aware compilation and operational information',30],[1,'12 Contract robustness, conditioning, and quotients',31],[2,'12.4 Finite randomized campaigns',33],
 [1,'13 Experimental methodology',33],[1,'14 Experimental results',36],[2,'14.9 Compiler-only scale on grid and heavy-hex graphs',41],[2,'14.10 Reproducibility status',42],
 [1,'15 Illustrative application interface: Lorentzian light-ray tomography',42],[1,'16 Limitations and open problems',43],[1,'17 Conclusion',44],
 [1,'Appendix A Notation and implementation conventions',45],[1,'Appendix B Proof details for the coherent-flag formula',45],[1,'Appendix C Proof details for matching necessity',46],[1,'Appendix D A small exact example',46],[1,'Appendix E Repository result schema and provenance',46],[1,'Appendix F Reproducibility and data availability',47],
 [1,'Appendix G Full unit-diagonal feasibility as background',48],[2,'G.1 The symmetric-Bernoulli correlation sector',48],[2,'G.2 Signed cats and the baseline they do not defeat',49],[2,'G.3 A positive-semidefinite matrix outside the full sector',49],
 [1,'References',50]
]
out.set_toc(toc)
out.save(OUT,garbage=4,deflate=True,clean=True)
print(OUT,'pages',out.page_count)
