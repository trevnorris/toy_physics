---
stage: 165
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [stage_165.md, stage_165.tex, moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py, moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.txt, moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl, moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.txt]
---
# Codex review — stage 165

## What the edit was supposed to do
F1 required the SymPy audit to stop short-circuiting the headline deliverable `d ln L_W = d ln a` with a literal `1`, and instead compute the logarithmic derivative of the stated `L_W` law with respect to `a`. F2 required the tautological Stage 164 fixed-g/fixed-r channel assertions to be demoted to clearly labelled solver-consistency prints in both SymPy and Mathematica, because those expressions are just the solved equations substituted back into themselves.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:41-44 | Computes `a*d/da log(LW_law)` and asserts `d ln L_W - d ln a = 0` | Yes; any non-linear `a` dependence in `LW_law` would leave a nonzero residual | Yes; directly checks the stage output `\delta\ln L_W=\delta\ln a` at fixed `r_*` | PASS |
| 2 | moving_throat_pde_stage165_exact_branch_drifts_sympy_audit.py:81-89 | Former fixed-g/fixed-r channel assertions demoted to solver-consistency prints | N/A; no longer asserted | Yes; it avoids claiming an independent proof of `delta_perp = 0` | PASS |
| 3 | moving_throat_pde_stage165_exact_branch_drifts_mathematica_audit.wl:67-74 | Former fixed-g/fixed-r channel assertions demoted to solver-consistency prints | N/A; no longer asserted | Yes; mirrors the required downgrade without reporting tautological PASS lines | PASS |

## Bottom line
PASS. I tried to break the new SymPy `d ln L_W` check by looking for a reused intermediate or `X-X` construction; it is instead a derivative of the independently assigned `LW_law` compared to `1`, so wrong `a` scaling would fail. I also checked whether the former channel closures still report as assertions; both scripts now print them only as `(solver consistency)`, and the saved outputs contain no `PASS:` for those tautological channels. I read the directive, the stage card, both scripts, and both saved outputs; the SymPy output shows `d ln L_W - d ln a at fixed r_* = 0`, while the Mathematica output still corroborates the same headline check and the demoted channel prints.
