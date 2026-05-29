---
stage: 160
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.py, moving_throat_pde_stage160_bare_mixed_port_slippage_sympy_audit.txt, moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl, moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.txt, stage_160.md, stage_160.tex]
---
# Codex review — stage 160

## What the edit was supposed to do
F1 required the Mathematica audit to stop being a line-by-line transliteration of the SymPy audit. The load-bearing derivation block had to remove the perturbation parameter, rename the internal symbols, replace `Series`/`Coefficient` extraction with direct total-differential formulas at the canonical point, and retarget the existing residual and pure-scale checks. The directive also required leaving the SymPy file and downstream tangential print-only block untouched.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl:44-45 | `deltaGammaW - (1/3) deltaKappaW` equals transported core residual | yes: changing the canonical `gamma0Canon = kappa0Canon/3` relation or the relative differential coefficients leaves a nonzero residual | yes: this is the compensated slippage identity whose `deltaKappaW = 0` gate yields the card's quoted `delta gamma_W` formula | PASS |
| 2 | moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl:47-49 | pure-scale harmlessness under `dGamma0 -> dKappa0/3` | yes: a wrong numerator coefficient, e.g. not `1/3`, would leave a residual | yes: checks the carry-forward corollary that the remaining odd defect vanishes for pure-scale tangent motion | PASS |
| 3 | moving_throat_pde_stage160_bare_mixed_port_slippage_mathematica_audit.wl:31-45 | second-engine derivation uses direct differential path, not SymPy series extraction | yes: retained `eps`, `Series`, `Coefficient`, or copied names would violate F1 | yes: derives the same paper-side residual from the canonical quotient formulas by a structurally different route | PASS |

## Bottom line
PASS. I tried to break the edited residual by treating the Mathematica block as if it had copied the SymPy perturbation expansion, by perturbing the canonical `1:3` relation between `gamma0Canon` and `kappa0Canon`, and by checking whether the pure-scale assertion was just `X - X`; those attacks fail against the current source. The `.wl` block no longer uses `eps`, `Series`, or `Coefficient`, the saved Mathematica output shows both residual checks evaluating to `0` with `PASS`, and the unchanged SymPy transcript corroborates the same target identity.
