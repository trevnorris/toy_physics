---
stage: 131
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [stage_131.md, stage_131.tex, moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py, moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.txt, moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl, moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.txt]
---
# Codex review — stage 131

## What the edit was supposed to do
The directive required replacing print-only or tautological checks with assertions for the quoted compensation point, slope, parent threshold identity, and lower-branch discrimination. It also required fixing the Mathematica banner from stage 114 to stage 131. Finally, both engines were supposed to cite/check the imported \(g_-^{F1}\) lower-branch closed form against the numeric literal.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:17 | \(g_-^{F1}\) closed form vs literal | Yes | Supports imported lower branch | OK |
| 2 | moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:42 | \(\Pi_*\) vs quoted value | Yes | Supports compensation point | OK |
| 3 | moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:50 | \(g'(\Pi_*)\) vs quoted slope | Yes | Supporting notes anchor, not the threshold itself | OK |
| 4 | moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:59 | threshold-at-star equals expected form | No meaningful physics failure | No; duplicates the target residual | FINDING R1 |
| 5 | moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:68 | \(g(2\Pi_*)\neq g_-\) | Yes | No; does not test singular-branch exclusion | FINDING R2 |
| 6 | moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:37 | \(g_-^{F1}\) closed form vs literal | Yes | Supports imported lower branch | OK |
| 7 | moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:55 | \(\pi_*\) vs quoted value | Yes | Supports compensation point | OK |
| 8 | moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:59 | slope at \(\pi_*\) vs quoted value | Yes | Supporting notes anchor, not the threshold itself | OK |
| 9 | moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:69 | threshold-at-star identity residual is zero | No meaningful physics failure | No; duplicates the target residual | FINDING R1 |
| 10 | moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:77 | \(g(2\pi_*)\neq g_-\) | Yes | No; does not test singular-branch exclusion | FINDING R2 |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:57`, `moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:63`
- **What's wrong:** SymPy checks `assert sp.simplify(threshold_at_star - expected_form) == 0`, where `threshold_at_star` is just `threshold_residual.subs(Pi, Pi_star)` and `expected_form` is the same residual formula with `Pi_star` copied in. Mathematica repeats the same construction with `identityResidual = Chop[Simplify[thresholdAtStar - expectedForm], 10^-30]`. This verifies self-consistency of a hardcoded residual, not the paper claim \(T_m-q_*A'_0=\Pi_*\Theta_\sigma/L\).
- **What it should be:** Derive the parent bias mismatch from the independent mouth-bias/source equations, then compare that derived expression at the paper \(\Pi_*\) to \((T_m-q_*A'_0)-\Pi_*\Theta_\sigma/L\), and separately show the residual vanishes under the canonical threshold relation.

### R2 — insufficient_verification
- **Where:** `moving_throat_pde_stage131_parent_mouth_threshold_sympy_audit.py:66`, `moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:76`
- **What's wrong:** The check `gPi(2*Pi_*) - g_-` only proves an arbitrary off-star point is not another lower-branch root. It does not compare \(\Pi_*\) against the singular equal-normalized branch, and \(\Pi_*\) itself was computed from `gPi == g_minus`, so lower-branch membership is not independently tested.
- **What it should be:** Evaluate the paper-quoted \(\Pi_*\) against the lower-branch condition and against an independently imported/derived singular-branch condition, asserting lower residual small and singular residual separated from zero.

### R3 — transliteration
- **Where:** `moving_throat_pde_stage131_parent_mouth_threshold_mathematica_audit.wl:33`
- **What's wrong:** The Mathematica audit is a direct transliteration of the SymPy targets: same \(g_-\) closed form, same literals, same root solve, same threshold tautology, and same \(2\pi_*\) off-star check. The saved Mathematica output corroborates execution, but it is not an independent second-engine derivation capable of catching a wrong target copied into both scripts.
- **What it should be:** One engine should independently derive or import the parent bias map and singular-branch discriminator rather than mirroring the same asserted targets.

## Bottom line
The most important problem is the parent-threshold assertion: both scripts report PASS while comparing the threshold residual to a restatement of itself. The saved outputs show the checks ran, but they do not establish the paper’s canonical branch equation or the lower-vs-singular branch discrimination required by the stage card.
