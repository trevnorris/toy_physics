---
stage: 107
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py, moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.txt, moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.wl, moving_throat_pde_stage107_general_dtn_deformation_mathematica_audit.txt, stage_107.md, stage_107.tex]
---
# Codex review — stage 107

## What the edit was supposed to do
F1 said the SymPy audit solved for `Sigma2` and `Sigma4` but only printed the resulting even-matching coefficients. The directive required two new `expect_zero` assertions locking those coefficients to the boxed closed forms already asserted by the Mathematica twin. The verifier expected new transcript lines `Sigma2 exact formula = 0` and `Sigma4 exact formula = 0` before the existing compensated `chi_Q` check.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py:60 | `Sigma2 exact formula` compares solved `Sigma2` to `-(3*S*beta**2 - 3*S + Sigma0)/9` | Yes. The RHS is not built from `sol[Sigma2]`; a sign, denominator, or DtN coefficient error leaves a nonzero residual. | Yes. It checks the compensated branch preserves the canonical `z^2` even coefficient. | PASS |
| 2 | moving_throat_pde_stage107_general_dtn_deformation_sympy_audit.py:64 | `Sigma4 exact formula` compares solved `Sigma4` to `-(3*S*beta**4 - 3*S + Sigma0)/27` | Yes. The RHS is independent of `sol[Sigma4]`; changing the target coefficient or solving equations breaks it. | Yes. It checks the compensated branch preserves the canonical `z^4` even coefficient. | PASS |

## Bottom line
PASS. I checked the new assertions against the directive, the stage card’s even-coefficient preservation requirement, and the saved transcripts. Independently deriving from `m2 = 1/9` gives `Sigma2 = S/3 - S beta^2/3 - Sigma0/9`; then using `m4 = 4/81` gives `Sigma4 = S/9 - S beta^4/9 - Sigma0/27`, matching the asserted forms. These are not `X - X` checks, and the SymPy output shows both new residuals as `0`; the Mathematica output also shows the corresponding checks passing, though I treated it as corroboration rather than the main proof because it follows the same algebraic structure.
