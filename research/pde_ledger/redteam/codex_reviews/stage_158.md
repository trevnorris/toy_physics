---
stage: 158
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [stage_158.md, stage_158.tex, moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py, moving_throat_pde_stage158_linear_defect_transport_sympy_audit.txt, moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl, moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.txt]
---
# Codex review — stage 158

## What the edit was supposed to do
The directive held the original paper-misalignment issue for resolution, while requiring two script-side fixes. F2 required deleting the tautological `delta Ms law` assertion and preserving the nontrivial `delta Mq law`. F3 required adding composed boxed checks for the carry-forward \(\delta M_q\) and \(\delta\Pi\) identities in both SymPy and Mathematica, with transcripts showing the new zero/PASS lines.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py:48 | `delta Ms law` removed; `delta Mq law` remains | yes, remaining Mq residual fails for wrong sign/coefficient | yes, \(\delta M_q\) transport | PASS |
| 2 | moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl:38 | `delta Ms law` removed; `delta Mq law` remains | yes, same residual structure | yes, \(\delta M_q\) transport | PASS |
| 3 | moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py:72 | composed `delta Mq` boxed law | yes, wrong \(R_*=1/4\) or wrong \(\delta g\) sign leaves nonzero residual | yes, boxed \(\delta M_q\) carry-forward | PASS |
| 4 | moving_throat_pde_stage158_linear_defect_transport_sympy_audit.py:79 | composed `delta Pi` boxed law | yes, wrong \(dS\), \(S_*\), or \(\delta g\) coefficient leaves nonzero residual | yes, boxed \(\delta\Pi\) carry-forward | PASS |
| 5 | moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl:61 | composed `delta Mq` boxed law | yes, same coefficient/sign attacks fail | yes, boxed \(\delta M_q\) carry-forward | PASS |
| 6 | moving_throat_pde_stage158_linear_defect_transport_mathematica_audit.wl:65 | composed `delta Pi` boxed law | yes, same coefficient/sign attacks fail | yes, boxed \(\delta\Pi\) carry-forward | PASS |

## Bottom line
PASS. I checked the directive, card, scripts, and saved transcripts. The old tautological `delta Ms law` line is absent from both transcripts. The new composed checks are not `x == x`: changing the \(\delta g\) sign, the \(1/4\) canonical \(R_*\) coefficient, or the \(\delta S\)/\(S_*\) coefficients would produce nonzero residuals. The Mathematica file mirrors the SymPy algebra rather than providing an independent derivation, but the copied target matches the directive and the stage card’s linearized \(\delta M_q,\delta\Pi\) deliverables, and both saved outputs show the expected `= 0` / `PASS` lines.
