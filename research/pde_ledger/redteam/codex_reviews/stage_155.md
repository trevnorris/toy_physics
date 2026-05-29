---
stage: 155
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py, moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.txt, moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl, moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.txt, stage_155.md, stage_155.tex]
---
# Codex review — stage 155

## What the edit was supposed to do
The directive had two fixes. F1 required adding SymPy-side numeric assertions for the load-bearing fixed-point moments `g_fp`, `S_fp`, `R_fp`, and `Pi_fp`, matching constants already checked by the Mathematica audit, because the prior SymPy assertion only tested an internally consistent transport-law identity. F2 required correcting the printed banner in both audit scripts from `STAGE 138` to `STAGE 155`, matching the paper card title. The directive explicitly said not to alter the Stage 154 docstring reference.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:80 | SymPy banner says Stage 155 | Yes: transcript would retain wrong stage label | Yes: matches stage card identity | PASS |
| 2 | moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:94 | `g_fp` equals fixed numeric target | Yes: computed moment is compared to external literal | Yes: anchors the co-evolving fixed-point output | PASS |
| 3 | moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:95 | `S_fp` equals fixed numeric target | Yes: computed outlet moment is compared to external literal | Yes: supports the reported fixed-point/Pi ledger | PASS |
| 4 | moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:96 | `R_fp` equals fixed numeric target | Yes: computed residual is compared to external literal | Yes: directly checks the paper quote `R_fp approx 0.282714` | PASS |
| 5 | moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:97 | `Pi_fp` equals fixed numeric target | Yes: computed `Sigma0*(1 - R*S)` is compared to external literal | Yes: anchors the printed fixed-point output | PASS |
| 6 | moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl:26 | Mathematica banner says Stage 155 | Yes: transcript would retain wrong stage label | Yes: matches stage card identity | PASS |

## Bottom line
The new SymPy assertions are not `X - X` checks: each compares a value computed from the fixed-point solve against a hard-coded target, so changing the map, quadrature, moment definitions, or constants would make them fail. I tried to break them by tracing shared intermediates, by checking whether `Pi_fp` was merely compared to itself, by mapping the targets back to the paper card’s fixed-point claim, and by checking the saved output order; the SymPy transcript prints the target values and then reaches later post-assert output, while the Mathematica transcript shows matching `PASS` lines and the corrected Stage 155 banner. The Mathematica file is a close algorithmic twin rather than a separate derivation, so I do not treat dual-engine agreement as independent proof, but the edited checks satisfy the directive and non-tautologically anchor the numerical stage deliverables.
