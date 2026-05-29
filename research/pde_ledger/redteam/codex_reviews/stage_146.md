---
stage: 146
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py, moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.txt, moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl, moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.txt, stage_146.md, stage_146.tex]
---
# Codex review — stage 146

## What the edit was supposed to do
F1 required replacing tautological affine-law checks with integral-form checks for the convex source moments, using a concrete positive normalized test profile and comparing the integrated moments against the affine law. F2 required anchoring the closed-form \(g(\Pi)\) and \(S_q(\Pi)\) formulas to direct integrals before the numeric \(\Pi_*\) work. F3 required correcting the Mathematica banner to the stage-146 first-order expansion title.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:42 | `g(Pi)` direct integral vs formula | Yes; direct integral is independent of `gPi`. | Yes | PASS |
| 2 | moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:44 | `S_q(Pi)` direct integral vs formula, numeric fallback | Yes, at sampled points only. | Partially | PASS |
| 3 | moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:50 | `S_q(Pi)` direct integral vs formula | Yes; Mathematica integrates independently. | Yes | PASS |
| 4 | moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:109 | `g_eps` affine law via integral form | Yes, but only with relaxed `1e-15` numeric tolerance. | Yes | FINDING |
| 5 | moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:110 | `S_eps` affine law via integral form | Yes, but only with relaxed `1e-15` numeric tolerance. | Yes | FINDING |
| 6 | moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:103 | `g_eps` affine law via integral form | Only for residuals larger than `10^-6`; raw residual is chopped in output. | Weakly | FINDING |
| 7 | moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:113 | `S_eps` affine law via integral form | Only for residuals larger than `10^-6`; raw residual is chopped in output. | Weakly | FINDING |
| 8 | moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:32 | stage banner literal | N/A | No | FINDING |

## Findings
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage146_positive_deformation_expansion_sympy_audit.py:116`
- **What's wrong:** `if abs(float(g_res)) > 1e-15:` relaxes the directive's numeric fallback tolerance from `1e-25` to `1e-15`. The saved output shows `g_eps` residuals around `10^-18`, so this would not satisfy the directive's own fallback standard even though the script prints `PASS`.
- **What it should be:** Use a symbolic zero check, or evaluate with enough precision that the numeric fallback can assert residuals below `1e-25`.

### R2 — insufficient_verification
- **Where:** `moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:106`
- **What's wrong:** `Print["g_eps affine law (integral form) at eps=1/10: ", fmt[Chop[gEpsSample1, 10^-6]]];` hides the raw residual, and lines 108 and 118 accept anything below `10^-6`. That is far weaker than the required `10^-25` tolerance and could pass a nonzero first-order error.
- **What it should be:** Print the raw residuals and require `Abs[...] < 10^-25`, or prove the residuals are exactly zero with `FullSimplify`.

### R3 — paper_misalignment
- **Where:** `moving_throat_pde_stage146_positive_deformation_expansion_mathematica_audit.wl:32`
- **What's wrong:** `banner["STAGE 146 — FINITE-CORRECTION EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];` does not match the directive or the paper card title. The saved Mathematica output corroborates the wrong banner text.
- **What it should be:** `banner["STAGE 146 — FIRST-ORDER EXPANSION FOR POSITIVE MOUTH-LAYER DEFORMATIONS"];`

## Bottom line
The main problem is that the edited affine-law checks were weakened during the numeric fallback: SymPy accepts `10^-15`, and Mathematica accepts and prints chopped zeros at `10^-6`. Those checks can miss nontrivial residuals, so the stage does not meet the directive's required non-tautological verification standard.
