---
stage: 161
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py, moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.txt, moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl, moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.txt, stage_161.md, stage_161.tex]
---
# Codex review — stage 161

## What the edit was supposed to do
The directive required replacing two tautological checks. F2 had to linearize the exact slippage expression \(B_W\), not a hand-built first-order polynomial. F1 had to derive \(\epsilon_\gamma\) from \(\gamma_0=(1+r_c)(1+\epsilon_\gamma)/9\), linearize it, and use the branch relation \(\gamma_{0,*}=(1+r_c)/9\) to obtain \(d\epsilon_\gamma=d\ln\gamma_0-d\ln(1+r_c)\). F3 corrected the mislabeled Stage 144 banner to Stage 161 in both transcripts.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py:52 | linearized slippage law from exact `BW` | Yes; changing `BW` or the kappa/gamma coefficient changes `dBW` | Yes, checks the D/N slippage scalar feeding \(\Xi_\gamma-2\Xi_L\) | PASS |
| 2 | moving_throat_pde_stage161_dn_similarity_slippage_sympy_audit.py:82 | \(d\epsilon_\gamma=d\ln\gamma_0-d\ln(1+r_c)\) | Yes; changing the exact gamma coefficient or omitting the branch relation leaves a residual | Yes, verifies the odd similarity-defect part of the final slippage coordinate | PASS |
| 3 | moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:43 | linearized slippage law from exact `bW` | Yes; `dBW` is differentiated from `bW`, not copied from the target | Yes, same paper-side slippage deliverable | PASS |
| 4 | moving_throat_pde_stage161_dn_similarity_slippage_mathematica_audit.wl:73 | \(d\epsilon_\gamma=d\ln\gamma_0-d\ln(1+r_c)\) | Yes; the branch substitution and log-variation substitution are algebraically testable | Yes, same odd-defect deliverable | PASS |

## Bottom line
The edited checks are no longer `X-X` comparisons. I tried to break F2 by perturbing the exact slippage source mentally; because `dBW` is now differentiated from `BW`, wrong coefficients would propagate into the residual. I tried to break F1 by changing the gamma coefficient or removing the branch relation; the residual would not simplify to zero. The saved SymPy output shows the Stage 161 banner, `linearized slippage law = 0`, and `d eps_gamma = d ln gamma0 - d ln(1+r_c) = 0`; the Mathematica transcript shows the same checks with `PASS`.
