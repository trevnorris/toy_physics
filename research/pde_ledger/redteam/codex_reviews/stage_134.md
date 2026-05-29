---
stage: 134
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py, moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.txt, moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl, moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.txt, stage_134.md, stage_134.tex]
---
# Codex review — stage 134

## What the edit was supposed to do
F1 required replacing the SymPy print-only body with substantive checks for the shell limit, the specialized \(S_q\) branch, \(S_q(\Pi_*)\), and the canonical gain line. F2 required replacing the Mathematica tautological specialized-kernel comparison with numeric checks against pasted literal targets. F3 was a paper-misalignment hold, while F4 was to be neutralized by the F1/F2 independent numeric literals.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:48 | `S_shell - 1 == 0` | Yes; denominator/sign changes break it. | Yes, shell-channel limit. | PASS |
| 2 | moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:58 | `S_q(1/2), S_q(1), S_q(2)` vs literals | Yes; the directive’s old placeholder values or a bad kernel would fail. | Yes, spot-checks the closed \(S_q\) branch. | PASS |
| 3 | moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:76 | `S_q(Pi_star)` vs notes value | Yes; bad \(S_q\) or \(\Pi_*\) changes the residual. | Yes, checks the recorded numerical value. | PASS |
| 4 | moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:82 | gain-line intercept/slope | No for the intercept; slope only restates check #3 with a minus sign. | No, not independently. | FINDING |
| 5 | moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:43 | static shell `expectZero` | Yes; bad shell limit fails. | Yes, shell-channel limit. | PASS |
| 6 | moving_throat_pde_stage134_family1_mouth_fixedpoint_mathematica_audit.wl:57 | `S_q` numeric literals | Yes; bad specialized branch fails. | Yes, spot-checks the closed \(S_q\) branch. | PASS |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage134_family1_mouth_fixedpoint_sympy_audit.py:82`
- **What's wrong:** `intercept = sp.N(Pi_star, 30)` is compared to `sp.Float("1.50882951349316", 30)`, the same literal already used to define `Pi_star`. The slope check is `slope = -S_star` compared to `-0.658075937605428`, which only repeats the prior `S_star` check with a sign flip. This cannot detect a wrong gain-line derivation or an unsupported gain selection; it only confirms constants already inserted into the script.
- **What it should be:** Either remove this as a substantive check, or compare an independently derived/stored gain-line expression or outlet-consistent gain pair against the fixed-point residual. If outlet consistency belongs to Stage 135 as the current card says, Stage 134 should not claim this as an independent gain-line verification.

## Bottom line
The shell and \(S_q\) numeric checks are corroborated by the saved outputs and are not runtime tautologies. The failure is the SymPy “canonical gain line” assertion: it is built from the same \(\Pi_*\) and \(S_q(\Pi_*)\) intermediates it claims to verify, so the stage still lacks a non-tautological check of that gain-line deliverable.
