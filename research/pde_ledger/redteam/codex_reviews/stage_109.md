---
stage: 109
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py, moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.txt, moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl, moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.txt, stage_109.md, stage_109.tex]
---
# Codex review — stage 109

## What the edit was supposed to do
The directive required applying F1 and F2 only. F1 adds a SymPy closed-form check that the preservation solution is \(a_5=-5b/9-a_0/27\), before the existing substitution check. F2 rewrites the Mathematica derivation so it expands numerator and denominator separately, forms the ratio through a denominator inverse series, and solves directly from `chiSeries - 1 == 0`. F3 was explicitly held as an unresolved paper-misalignment issue and was not supposed to be fixed in this pass.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:54 | `a5 preservation closed-form` | yes; compares solve output to handwritten closed form | yes; branch preservation condition implied by \(\chi_Q=1\) | PASS |
| 2 | moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:55 | `preservation substitution` | no; substitutes the value solved from `coeff == 0` back into `coeff` | no independent exercise | FINDING |
| 3 | moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:51 | `linearized chi law` | yes; independent ratio expansion is compared to handwritten expected law | yes; checks quoted \(\Delta_Q=5b+a_0/3+9a_5\) | PASS |
| 4 | moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:57 | `overall scale cancels` | yes; would fail if the first-order defect retained `s` | partially; covers the pure-scale cancellation piece | PASS |
| 5 | moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:62 | `a5 preservation condition + 5 b/9 + a0/27` | yes; compares solve output to a closed-form target | yes; checks compensated branch selection | PASS |
| 6 | moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:63 | `preservation substitution` | no; substitutes the value solved from `chiSeries - 1 == 0` back into `chiSeries - 1` | no independent exercise | FINDING |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage109_linearized_branch_selection_sympy_audit.py:55`
- **What's wrong:** `expect_zero('preservation substitution', coeff.subs(a5, a5_sol))`. `a5_sol` is defined at line 51 by solving `coeff == 0` for `a5`, so substituting it back into the same `coeff` cannot falsify a wrong physical coefficient. The saved output line `preservation substitution = 0` only corroborates this circular check.
- **What it should be:** Substitute the independent closed-form target, e.g. `expect_zero('preservation with closed-form branch', sp.simplify(coeff.subs(a5, expected_a5_sol)))`, or remove this as evidence and rely on the closed-form check.

### R2 — tautological_check
- **Where:** `moving_throat_pde_stage109_linearized_branch_selection_mathematica_audit.wl:63`
- **What's wrong:** `expectZero["preservation substitution", (chiSeries - 1) /. a5 -> a5Pres];`. `a5Pres` is defined at line 60 by `Solve[chiSeries - 1 == 0, a5, Reals]`, so the check substitutes a solved root into the same equation used to define it. The saved output lines 14-15 show it passes, but the pass is not independent evidence.
- **What it should be:** Substitute the closed-form branch directly, e.g. `(chiSeries - 1) /. a5 -> (-5*b/9 - a0/27)`, or treat line 62 as the actual non-tautological preservation check.

## Bottom line
The new closed-form checks are corroborated by the saved outputs, and the Mathematica ratio derivation is no longer a direct line-by-line copy of the SymPy `coeff = (chiSeries - 1)/eps` path. The blocking problem is that both scripts still present preservation-substitution assertions that are solved-equation self-substitutions; those checks can stay green even when the only thing they prove is that the algebra system substituted its own root back into its own equation.
