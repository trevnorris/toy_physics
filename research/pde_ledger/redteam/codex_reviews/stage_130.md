---
stage: 130
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py, moving_throat_pde_stage130_mouth_bias_map_sympy_audit.txt, moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl, moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.txt, stage_130.md, stage_130.tex]
---
# Codex review — stage 130

## What the edit was supposed to do
The directive required two fixes. F2 adds a direct equality check that the integral-defined mouth-bias map equals the boxed closed form \(g_\Pi=2\Pi(2\Pi e^\Pi+\pi)/((4\Pi^2+\pi^2)(e^\Pi-1))\). F1 adds a strict monotonicity sweep for \(dg_\Pi/d\Pi>0\) at six positive sample points bracketing \(\Pi_*\), because the paper card claims exact \(g_\Pi\), monotonicity, and a unique Family-1 compensation point.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:16 | `gPi` equals boxed closed form | yes | yes, exact \(g_\Pi\) | PASS |
| 2 | moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:43 | `gPi` equals boxed closed form | yes | yes, exact \(g_\Pi\) | PASS |
| 3 | moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:40 | finite sweep of `dg/dPi > 0` | yes, but only at six points | only partially; not monotonicity on \(\Pi>0\) or uniqueness | FINDING |
| 4 | moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:63 | finite sweep of `dg/dpiM > 0` | yes, but only at six points | only partially; not monotonicity on \(\Pi>0\) or uniqueness | FINDING |

## Findings
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage130_mouth_bias_map_sympy_audit.py:40`; `moving_throat_pde_stage130_mouth_bias_map_mathematica_audit.wl:63`
- **What's wrong:** The added checks test `dgPi` only at `{1/10, 1/2, 1, 15088/10000, 3, 10}`. That can fail if one sampled derivative is non-positive, so it is not tautological, and the saved outputs do show positive numerical values / PASS lines. But it does not prove the paper-side monotonicity claim \(dg_\Pi/d\Pi>0\) for all \(\Pi>0\), nor the resulting uniqueness of \(\Pi_*\). A wrong implementation with a negative derivative between the sampled points, or for \(\Pi<0.1\) or \(\Pi>10\), would pass this edit.
- **What it should be:** Add an exact sign certificate for `D[gPi, Pi]` under `Pi > 0`, or otherwise reduce the derivative to a numerator/denominator form whose denominator is positive and whose numerator is proven positive on the whole positive domain. Then use endpoint values plus strict monotonicity to justify uniqueness of the compensation point.

## Bottom line
The boxed-form equality checks are substantive: I independently checked the integral structure against the closed form, and the Mathematica transcript explicitly reports zero/PASS while the SymPy script would stop before printing if the assertion failed. The monotonicity edits are the weak point. They satisfy the directive’s literal “multi-point sweep,” but they do not verify the paper card’s monotonicity or unique-point claim; the most important problem is that finite sampling is being treated as evidence for a global strict inequality.
