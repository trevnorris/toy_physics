---
stage: 157
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py, moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.txt, moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl, moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.txt, stage_157.md, stage_157.tex]
---
# Codex review — stage 157

## What the edit was supposed to do
The directive required removing the old tautological `-16 sigmaStar * dRfamily` delta-C assertion, because it only multiplied a previously checked zero. It also required the Mathematica audit to stop merely transliterating the SymPy algebra and to derive at least one load-bearing delta-C/even-preservation step independently, preferably through a 2x2 canonical-even solve. The saved outputs were expected to show the new solve-derived check and refreshed passing output.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:112 | repeat `solve(dE2=0,dE4=0)` and extract `deltaC` | Not after the immediately preceding full-solution check, except by solver inconsistency | Only duplicates canonical-even preservation; it does not test the tangent/family handoff | FINDING |
| 2 | moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:114 | `expect_zero("delta C from canonical-even Solve", deltaC_from_pair)` | No independent fail mode beyond line 109 | Adjacent to the paper checklist, but not a distinct verification of `delta_perp=0` | FINDING |
| 3 | moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:103 | hard-coded 2x2 solve for `dCsym,dKsym` | Yes for some coefficient errors, but not for copied wrong targets with unique zero solution | Tests the same canonical-even algebra as SymPy, not an independent derivation | FINDING |
| 4 | moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:105 | `expectZero["delta C from canonical-even Solve", deltaCIndep]` | No independent fail mode after the same solve | Corroborated by output, but not independent and not tied to the tangent family | FINDING |

## Findings (only if verdict = FINDINGS)
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage157_core_mouth_coevolution_status_sympy_audit.py:112`
- **What's wrong:** `sol_deltaC = sp.solve([sp.Eq(dE2, 0), sp.Eq(dE4, 0)], [deltaC, dkappa], dict=True)[0]`. Lines 107-110 already solved the same equations and asserted the whole solution is `{deltaC: 0, dkappa: 0}`. The edited assertion only repeats that solve and extracts one component, so after the preceding check passes it adds no new way to catch a wrong tangent-to-defect map or a wrong paper target.
- **What it should be:** Either make the solve check the single canonical-even preservation assertion, or derive `deltaC`/`delta_perp` through a separate tangent-family transport expression and compare that to the canonical-even result.

### R2 — transliteration
- **Where:** `moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:103`
- **What's wrong:** `solDeltaC = Solve[{dCsym - 9 sigmaStar dKsym == 0, 5 dCsym - 72 sigmaStar dKsym == 0}, {dCsym, dKsym}];`. This is not an independent second-engine derivation; it is the same hard-coded numerator system used by the SymPy-side solve. If both engines copy a wrong canonical-even target that still has unique zero solution, both pass.
- **What it should be:** Derive the canonical-even coefficients in Mathematica from an upstream Galerkin/canonical source, or otherwise compare against an independently constructed expression rather than restating the same 2x2 literals.

### R3 — symbol_assumption_error
- **Where:** `moving_throat_pde_stage157_core_mouth_coevolution_status_mathematica_audit.wl:93`
- **What's wrong:** `$Assumptions = Element[{sigmaStar, deltaC, dKappa}, Reals];`. The directive’s proposed independent solve explicitly required `sigmaStar > 0 && sigmaStar < 1`. Without that branch restriction, the full even-preservation claim is not valid at `sigmaStar = 0` for `dKappa`, and the formulas are singular at `sigmaStar = 1`.
- **What it should be:** Encode the physical branch domain, at minimum `0 < sigmaStar < 1`, or make the exceptional cases explicit.

## Bottom line
The saved outputs do show the new checks passing: SymPy prints `delta C from canonical-even Solve = 0`, and Mathematica prints the same with `PASS`. But the fix replaced the old known-zero multiplication with a repeated same-equation solve, and the Mathematica version remains a mirrored hard-coded target rather than an independent derivation. The most important problem is that these edits do not create an independent failure mode for the co-evolving tangent/even-preservation handoff claimed by the stage card.
