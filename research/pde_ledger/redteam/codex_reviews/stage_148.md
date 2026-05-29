---
stage: 148
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [stage_148.md, stage_148.tex, moving_throat_pde_stage148_representative_positive_families_sympy_audit.py, moving_throat_pde_stage148_representative_positive_families_sympy_audit.txt, moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl, moving_throat_pde_stage148_representative_positive_families_mathematica_audit.txt]
---
# Codex review — stage 148

## What the edit was supposed to do
The directive required replacing the old hardcoded or same-source `xi_*` comparison with the independent Stage 228 closed form `(-37√3 - 5π² + 2√(4107 - 168π²))/(5(8 - π²))`, then asserting `1 - λ_{Π,0} = ξ_*`. It also required the Mathematica audit to stop transliterating the SymPy `aT,bT` path by deriving `dT` through `dSigma`, while preserving the same printed `dTU`, `dTD`, and `dTLam` numerics as SymPy. The paper card’s relevant deliverable is that the uniform and self-matched derivative families bracket the first-order correction and interpolation reproduces the earlier compensation fraction.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:87-92 | Numeric residual `(1-lambda_(Pi,0)) - xi_*` | Only weakly; wrong target and loose tolerance | No, it uses `4107 - 100*pi^2`, not the mandated Stage 228 target | FAIL |
| 2 | moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:84-86 | Exact `expectZero` for `(1-lambda_(Pi,0)) - xi_*` | No for the intended bridge; it is the same-source closed form of `gMinus` | No, it verifies the wrong closed form | FAIL |
| 3 | moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:43-73 | Reworked `dSigma`/`dT` derivation | Yes; it did fail against SymPy numerics | No, the first-order correction values diverge from the paper audit target | FAIL |

## Findings (only if verdict = FINDINGS)
### R1 — paper_misalignment
- **Where:** `moving_throat_pde_stage148_representative_positive_families_sympy_audit.py:87`; `moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:84`
- **What's wrong:** `sqrt(4107 - 100*sp.pi**2)` / `Sqrt[4107 - 100*Pi^2]` is not the directive’s required Stage 228 closed form, which uses `4107 - 168π²`. With `100π²`, the expression is just the algebraic closed form induced by the same `rF1`/`gMinus` used to define `lamPiZero`, so it does not test the independent compensation bridge. The SymPy output also shows a nonzero residual `7.82010941918532143843153844223e-17`, accepted only because the assertion was weakened to `1e-15` instead of the required `1e-25`.
- **What it should be:** Use `(-37*sqrt(3) - 5*pi^2 + 2*sqrt(4107 - 168*pi^2))/(5*(8 - pi^2))` in both engines, and require the stronger exact or high-precision residual check specified by the directive.

### R2 — insufficient_verification
- **Where:** `moving_throat_pde_stage148_representative_positive_families_mathematica_audit.wl:43-73`
- **What's wrong:** The reworked Mathematica `dTOfDeltas` path no longer reproduces the SymPy first-order correction values. The directive explicitly says the Mathematica `dTU`, `dTD`, and `dTLam` numerics must agree with SymPy to printed precision, but the saved outputs disagree: SymPy prints `dTU = 0.5087563022150839...`, `dTD = -0.1169438021518107...`, while Mathematica prints `dTU = 0.4976692633908835...`, `dTD = -0.1144451420715406...`. This means the independent engine is not corroborating the same correction.
- **What it should be:** The Mathematica `dSigma` route must include all first-order contributions needed to match the SymPy `AT/BT` correction, then the transcript should preserve the SymPy `dT` coefficients to printed precision.

## Bottom line
The edits do not pass adversarial review. The load-bearing `xi_*` check uses the wrong closed form and is still same-source with the neutral-point algebra, while the Mathematica rework breaks the first-order correction numerics it was supposed to independently corroborate.
