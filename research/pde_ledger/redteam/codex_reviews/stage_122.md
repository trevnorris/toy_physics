---
stage: 122
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [audit_directive_122.md, verify_directive_122.md, stage_122.tex, moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py, moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.txt, moving_throat_pde_stage122_mouth_source_compensation_test.md]
---
# Codex review — stage 122

## What the edit was supposed to do
The listed directive path `redteam/directives/stage_122.md` is absent, so I used the matching directive material in `redteam/tmp_prompts/audit_directive_122.md` plus the verifier note. The fix was supposed to resolve the notes/script `168π²` versus `100π²` mismatch in favor of the script frame, then add six SymPy assertions for the compensation quadratic, natural-branch defect, off-surface result, and traction-ratio identities. No Mathematica script was required.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:56 | `gminus` satisfies compensation quadratic | Yes, for algebraic sign/factor errors in the branch formula or quadratic. | Yes: notes §2 / compensation-family roots. | PASS |
| 2 | moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:57 | `gplus` satisfies compensation quadratic | Yes, for algebraic sign/factor errors in the branch formula or quadratic. | Yes: notes §2 / compensation-family roots. | PASS |
| 3 | moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:60 | natural defect closed form | Yes, if the boxed closed form or `rF` frame is wrong. | Yes: notes §3 defect formula. | PASS |
| 4 | moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:62 | natural branch off compensation surface | Weakly yes: it fails if `comp_def` simplifies to literal zero. | Yes: headline off-surface claim. | PASS |
| 5 | moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:64 | lower traction ratio identity | No: `T_ratio_minus` is defined as `1/gminus`. | No: it only checks reciprocal algebra, not an independently derived traction law. | FINDING |
| 6 | moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:65 | upper traction ratio identity | No: `T_ratio_plus` is defined as `1/gplus`. | No: it only checks reciprocal algebra, not an independently derived traction law. | FINDING |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage122_mouth_source_compensation_test_sympy_audit.py:43-44,64-65`
- **What's wrong:** `T_ratio_minus = sp.simplify(1/gminus)` and `T_ratio_plus  = sp.simplify(1/gplus)` are later asserted by `expect_zero("traction ratio (-) identity", gminus * T_ratio_minus - 1)` and `expect_zero("traction ratio (+) identity", gplus  * T_ratio_plus  - 1)`. These checks reduce to `gminus*(1/gminus)-1` and `gplus*(1/gplus)-1` by construction. If the Stage 221 proportionality between traction amplitude and mouth-coupling ratio were misquoted, missing a normalization factor, or had the wrong branch normalization, these assertions would still pass.
- **What it should be:** Derive `T_m^{(\pm)}/T_m^{nat}` from the imported traction relation or a separately encoded source formula, then compare that derived expression to `1/gminus` and `1/gplus`. At minimum, the checked quantity must not be defined as the reciprocal immediately before the assertion.

## Bottom line
The compensation-quadratic and defect checks are corroborated by the saved output lines `= 0`, and the nonzero natural-defect check is present with the exact nonzero expression in the transcript. The traction-ratio checks are not trustworthy: I tried to break them by imagining the physical traction law carrying a missing factor or wrong normalization, and the script would still pass because the asserted ratios are defined as `1/g±` before being checked. That makes the edit fail the non-tautology requirement despite the saved output showing `traction ratio ... identity = 0`.
