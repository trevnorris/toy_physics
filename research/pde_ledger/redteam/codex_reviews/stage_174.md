---
stage: 174
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage174_static_self_similarity_sympy_audit.py, moving_throat_pde_stage174_static_self_similarity_sympy_audit.txt, moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl, moving_throat_pde_stage174_static_self_similarity_mathematica_audit.txt, stage_174.md, stage_174.tex]
---
# Codex review — stage 174

## What the edit was supposed to do
F1 required fixing the Mathematica audit because it was a line-for-line SymPy port with the first-order differentials `b01Two`, `z01Two`, and `n01Two` hand-typed, making transcription errors invisible to the second engine. The required change was to derive those differentials from perturbed primitive expressions using `D[..., eps] /. eps -> 0`, while leaving the `expectZero` checks, assumptions, weighted-slope targets, theorem substitutions, and carry-forward block unchanged. It also required correcting the stale banner to `STAGE 174`.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl:67-76 | BdG differential from `b0Eps`, checked against weighted log-slope formula | Yes. A wrong sign, factor 2, or denominator perturbation would make `b01Two/b0Two - deltaBWeighted` nonzero; it is not `X - X`. | Yes, this is the BdG slope piece in the weighted static self-similarity failure. | PASS |
| 2 | moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl:86-95 | Conservative Maxwell/mixed differential from `z0Eps`, checked against weighted quotient log-slope formula | Yes. Quotient-rule sign or denominator errors are independently exposed by the perturbation derivative. | Yes, this is the conservative port slope piece in the stage output. | PASS |
| 3 | moving_throat_pde_stage174_static_self_similarity_mathematica_audit.wl:105-114 | Outgoing-transfer differential from `n0Eps`, checked against weighted squared-ratio log-slope formula | Yes. Dropping the squared-ratio factor 2 or the `delta` variation would leave a residual. | Yes, this is the outgoing port slope piece in the stage output. | PASS |

## Bottom line
PASS. I checked the directive, the stage card, both scripts, and both saved transcripts. The edited Mathematica checks now derive the three first-order differentials from perturbed primitive expressions rather than retyping the SymPy closed forms, so the attacks I tried mentally — sign flips, missing factor 2, wrong quotient denominator variation, or omitted perturbation terms — would not cancel tautologically. The Mathematica output corroborates all three edited checks with `= 0` and `PASS`, and the banner is now `STAGE 174`; the SymPy transcript shows the corresponding audit identities evaluating to `0`.
