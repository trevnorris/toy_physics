---
stage: 152
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py, moving_throat_pde_stage152_family1_actual_correction_sympy_audit.txt, moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.wl, moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.txt, stage_152.md, stage_152.tex]
---
# Codex review — stage 152

## What the edit was supposed to do
The directive identified that the SymPy audit computed the four load-bearing scale deliverables, `delta Pi_act`, `delta Tm_act`, `lambda_eff^(Pi)`, and `lambda_eff^(T)`, but only printed them. The required fix was to add four anchored `expect_close` assertions against the notes-quoted constants, mirroring the existing Mathematica `expectApprox` checks. The existing covariance consistency check was to remain untouched.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:117 | `deltaPi` equals `0.907084414842908` | Yes; computed covariance/g-prime result is compared to an independent literal target. | Yes; this is the finite `Pi` correction determining `Pi_corr`. | PASS |
| 2 | moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:118 | `deltaT` equals `0.271653979462338` | Yes; computed `AT*delta_g + BT*delta_S` is compared to an independent literal target. | Yes; this is the finite mouth-temperature correction determining `T_corr`. | PASS |
| 3 | moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:119 | `lam_Pi` equals `0.380487632771110` | Yes; computed affine interpolation parameter is compared to an independent literal target. | Yes; it quantifies the claimed moderate broadening scale for the `Pi` correction. | PASS |
| 4 | moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:120 | `lam_T` equals `0.378939241176339` | Yes; computed affine interpolation parameter is compared to an independent literal target. | Yes; it quantifies the claimed moderate broadening scale for the `T_m` correction. | PASS |

## Bottom line
PASS. I tried to break the new checks as tautologies by tracing whether each target was rebuilt from the same intermediate; it was not, since all four assertions compare computed quantities to hard-coded numerical constants. I checked paper alignment against the stage card’s finite mouth-profile correction claim: `deltaPi` and `deltaT` directly set the corrected point, while the two `lambda_eff` checks audit the stated moderate broadening toward the uniform family. The saved SymPy output shows all four new checks present with residuals below tolerance, and the Mathematica transcript independently corroborates the same four scales while deriving the canonical constants and interpolation baselines rather than merely hardcoding the SymPy setup.
