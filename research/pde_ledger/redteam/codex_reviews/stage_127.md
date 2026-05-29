---
stage: 127
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage127_penetration_families_sympy_audit.py, moving_throat_pde_stage127_penetration_families_sympy_audit.txt, moving_throat_pde_stage127_penetration_families_mathematica_audit.wl, moving_throat_pde_stage127_penetration_families_mathematica_audit.txt, stage_127.md, stage_127.tex]
---
# Codex review — stage 127

## What the edit was supposed to do
The directive identified a paper/documentation stage-number mismatch as F1 and held that as a paper-misalignment item, with Stage 127 as the recommended canonical direction if resolved. The actionable mathematical fix was F2: the Mathematica audit was a line-by-line port of the SymPy audit and needed independent symbolic-integration checks for the slab and truncated-exponential mouth-bias closed forms. The new checks were required to derive `g_slab` and `g_exp` from their source integrals, print two new `PASS:` lines, and leave the existing compensation-root checks intact.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:36-41 | derive slab bias by `Integrate[(1/x) Cos[Pi z/2], {z,0,x}]` and compare to `2 Sin[Pi x/2]/(Pi x)` | Yes; a wrong factor, sign, limit, or normalization leaves a nonzero residual | Yes; anchors the slab mouth-bias map used for the penetration family reaching the lower branch | PASS |
| 2 | moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:44-49 | derive exponential bias by integrating the normalized truncated exponential source against the boundary-layer cosine and compare to closed form | Yes; the integral is built independently from the closed form, so denominator/sign/exponential mistakes would fail | Yes; anchors the truncated-exponential mouth-bias map used for the compensation point | PASS |
| 3 | moving_throat_pde_stage127_penetration_families_mathematica_audit.wl:61-62 | existing slab/exponential compensation roots still satisfy `g == gMinus` | Yes numerically, if `FindRoot` converged to a wrong value or the formulas no longer hit `gMinus` | Yes; directly checks the card’s “reach the lower branch at moderate penetration depths” claim | PASS |

## Bottom line
PASS. I tried to break the added checks by looking for same-intermediate `X - X` construction, copied closed-form reuse, unjustified assumptions, and output gaps. The Mathematica additions instead integrate the source densities directly and subtract separately defined closed forms under the justified penetration-depth domain `x > 0`; the saved output shows both new `PASS:` lines plus the original slab and exponential compensation-root `PASS:` lines. The remaining SymPy audit is not the independent engine for the new derivation, but the directive only required that independence to be added on the Mathematica side, which it was.
