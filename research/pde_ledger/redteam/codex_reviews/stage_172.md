---
stage: 172
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py, moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.txt, moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl, moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.txt, stage_172.md, stage_172.tex]
---
# Codex review — stage 172

## What the edit was supposed to do
F1 required the Mathematica audit to stop extracting the four slopes by the same `Normal[Series[...]]/(eps*lam)` route used in SymPy, and instead derive them by implicit differentiation of the defining relations. The seven final residual assertions and printed forms were supposed to stay unchanged. F2 only corrected stale `STAGE 155` transcript banners to `STAGE 172`. The paper card’s output is the collapse identifying \(\mathfrak K_1=-D_0u_2^{(1)}\) and \(\mathfrak G_1=D_0P_1\), with the surrounding physical even-response reductions preserved.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:43 | `deltaU2` from implicit derivative of \(u_2D_0+D_2=0\) | Yes: sign or base-relation errors would change the printed slope and break the \(K_A\) residuals. | Yes, feeds \(\mathfrak K_1=-D_0u_2^{(1)}\). | PASS |
| 2 | moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:45 | `deltaP0` from implicit derivative of \(P_0D_0-N_0=0\) | Yes: numerator/sign errors would break `G_A - D0*delta_P0`. | Yes, directly verifies \(\mathfrak G_1=D_0P_1\). | PASS |
| 3 | moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:79 | canonical `deltaU2Star` by implicit derivative | Yes: wrong canonical branch data would break the hidden-even residual. | Yes, supports the physical even-response collapse. | PASS |
| 4 | moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:81 | `deltaU4Star` from implicit \(u_4D_0^2-(D_2^2-D_0D_4)=0\) | Yes: coefficient/sign changes leave a nonzero hidden-even residual. | Yes, checks the even-consistency reduction used by the stage. | PASS |
| 5 | moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:55 | seven residual assertions retained with independent slopes | Yes: they are not `X-X`; wrong slope, branch, or prefactor choices would leave symbolic residuals. | Yes, they cover the obstruction-pair and direct outlet physical-slope reductions. | PASS |
| 6 | moving_throat_pde_stage172_physical_slope_collapse_sympy_audit.py:31; moving_throat_pde_stage172_physical_slope_collapse_mathematica_audit.wl:26 | transcript banner says `STAGE 172` | Yes at transcript level; stale labels would show in saved output line 3. | No, print-only directive hygiene. | PASS |

## Bottom line
I tried to break the edit by treating the Mathematica slope derivations as copied quotient expansions, by checking whether the residuals were built as `X-X`, and by perturbing the relevant signs/coefficient roles mentally against the displayed formulas. Those attacks did not hold: the Mathematica slopes now come from `Solve[D[...]]` on implicit defining relations, the residuals compare those slopes against separately named obstruction/direct-outlet expressions, and both saved transcripts show the corrected `STAGE 172` banner plus the required zero/PASS lines.
