---
stage: 119
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [moving_throat_pde_stage119_parent_balance_sympy_audit.py, moving_throat_pde_stage119_parent_balance_sympy_audit.txt, moving_throat_pde_stage119_parent_balance_mathematica_audit.wl, moving_throat_pde_stage119_parent_balance_mathematica_audit.txt, stage_119.md, stage_119.tex]
---
# Codex review — stage 119

## What the edit was supposed to do
F1 was supposed to stop the tube-length law from floating free of the parent-balance family by tying the Section III parameter `rc`/`rC` to the normalized hybridization parameter `rhat`/`rHat` via `r_c = rhat^2`. F2 was supposed to add assertions that the solved positive and negative branch traction scales `T_m` match the notes-given closed forms. The Mathematica fix was also instructed to strip `ConditionalExpression` wrappers before comparing those `T_m` formulas.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage119_parent_balance_sympy_audit.py:53 | `rc -> rhat**2` tube-length link | Only if the prior tube law fails; it assumes the link | Partially, but does not derive the link | FINDING |
| 2 | moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:60 | `rC -> rHat^2` tube-length link | Only if the prior tube law fails; it assumes the link | Partially, but does not derive the link | FINDING |
| 3 | moving_throat_pde_stage119_parent_balance_sympy_audit.py:75 | `T_m` plus branch closed form | Algebraically yes; domain-wise no | Only value form, not branch existence | FINDING |
| 4 | moving_throat_pde_stage119_parent_balance_sympy_audit.py:80 | `T_m` minus branch closed form | Algebraically yes; domain-wise no | Only value form, not branch existence | FINDING |
| 5 | moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:90 | `T_m` plus branch closed form after `stripCE` | Algebraically yes; domain-wise no | Only value form, not branch existence | FINDING |
| 6 | moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:95 | `T_m` minus branch closed form after `stripCE` | Algebraically yes; domain-wise no | Only value form, not branch existence | FINDING |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage119_parent_balance_sympy_audit.py:53`, `moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:60`
- **What's wrong:** `L_sel.subs(rc, rhat**2) - sp.pi*a*sp.sqrt((1+rhat**2)/3)/2` and `(lSel /. rC -> rHat^2) - (Pi*a*Sqrt[(1 + rHat^2)/3])/2` are just substitution instances of the immediately preceding tube-length identity. They do not derive or test `r_c = rhat^2`; they assume it. If the family connection were absent or wrong, these checks would still print `0`/`PASS`.
- **What it should be:** Independently derive the core ratio used as `r_c` from the Section I-II normalized variables and assert `rc_expr - rhat**2 == 0`, then use that derived expression in the tube-length law.

### R2 — symbol_assumption_error
- **Where:** `moving_throat_pde_stage119_parent_balance_sympy_audit.py:75`, `moving_throat_pde_stage119_parent_balance_sympy_audit.py:80`, `moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:88`
- **What's wrong:** `stripCE[expr_] := expr /. ConditionalExpression[v_, _] :> v;` discards the branch conditions Mathematica found. The saved output shows the plus branch requires `rHat > -1/Sqrt[3]` and the minus branch requires `rHat > 1/Sqrt[3]`. With only `rhat`/`rHat` real, the asserted formulas can be negative or invalid for positive `T_m` outside those ranges, e.g. plus at `rhat = -2` or minus at `rhat = 0`.
- **What it should be:** Keep or separately assert the branch conditions: plus branch `rhat > -1/sqrt(3)`, minus branch `rhat > 1/sqrt(3)`, or restrict the section assumptions to the intended physical branch before claiming positive `T_m`.

### R3 — transliteration
- **Where:** `moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:60`, `moving_throat_pde_stage119_parent_balance_mathematica_audit.wl:90`
- **What's wrong:** The Mathematica edits are a line-by-line copy of the SymPy targets: same assumed `rC -> rHat^2` substitution and same manually typed `T_m` branch formulas. This gives a second simplifier, not an independent derivation. A copied wrong target would pass in both engines.
- **What it should be:** The Mathematica audit should independently derive the family link and preserve/check the `ConditionalExpression` domains instead of erasing them to match the copied target.

## Bottom line
The saved outputs do corroborate the new lines as `0`/`PASS`, but the important F1 check is still tautological: I tried to break it by imagining the actual family relation changed or disappeared, and the assertion still passes because it simply substitutes the desired relation into an already-proved identity. The `T_m` checks also pass while ignoring the branch-existence conditions that Mathematica itself printed. The most important problem is that the edits verify algebraic restatements, not the family linkage and positive-branch domain needed for the stage.
