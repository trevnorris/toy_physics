---
stage: 151
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py, moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.txt, moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl, moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.txt, stage_151.md, stage_151.tex]
---
# Codex review — stage 151

## What the edit was supposed to do
The directive required replacing tautological covariance checks with an integration-based derivation of the moment shifts from the centered first-order source correction. It also required adding actual assertions for the bias retuning and traction retuning, rather than merely printing them. The Mathematica audit was supposed to be independent of the SymPy path by deriving `deltaSigma` from a `Series` expansion of the normalized full mouth source.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:80-81 | numeric canonical normalization | Yes, for the one numeric anchor | Adjacent only | OK |
| 2 | moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:104-105 | numeric centering of `delta_Sigma` | Yes, for the one numeric anchor | Only as a numeric sample | FINDING R1/R2 |
| 3 | moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:108-111 | numeric `delta_g_int`, `delta_S_int` versus covariances | Yes, for the one numeric anchor | Only as a numeric sample, not the exact stage claim | FINDING R1/R2 |
| 4 | moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:114-122 | numeric bias and traction retuning | Yes, for the one numeric anchor | Only as a numeric sample | FINDING R1/R2 |
| 5 | moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:53 | Series-derived centering | Yes | Yes | OK |
| 6 | moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:54 | `deltaSigma` agrees with centered hand form | Yes | Yes | OK |
| 7 | moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:56-59 | integrated moment shifts versus covariances | Yes | Yes | OK |
| 8 | moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:61-64 | bias and traction retuning assertions | Yes | Yes | OK |

## Findings
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:43-116`
- **What's wrong:** The directive required a SymPy symbolic derivation with `Pi_star`, `r1`, and `r2` symbolic. Instead the script uses `mpmath`, fixes `Pi_star`, `r1`, `r2`, `gprime`, `AT`, and `BT` to numeric anchors, and verifies the identities only by quadrature at that one point. A wrong formula that happens to work for this selected residual/source parameter set would pass, so this is not faithful to the paper card’s exact covariance-correction claim.
- **What it should be:** Restore the required symbolic `sp.symbols`, `sp.integrate`, and `expect_zero` checks for general positive `Pi_star` and symbolic real `r1`, `r2`, including the two retuning assertions.

### R2 — insufficient_verification
- **Where:** `moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.txt:11-15`
- **What's wrong:** The saved SymPy transcript does not show the required exact `= 0` / `PASS` residual checks. It shows tolerance-based approximate comparisons with nonzero diffs, for example centering diff `1.3758475920694696e-42`. That corroborates only that the numeric tolerance accepted the sample, not that the exact symbolic residuals vanished.
- **What it should be:** The SymPy transcript should contain exact zero residual lines from the symbolic `expect_zero` assertions specified in the directive.

## Bottom line
The Mathematica edit follows the requested independent Series path and its saved output corroborates the zero checks. The SymPy edit does not: it replaced the required symbolic SymPy audit with a numeric `mpmath` smoke test over one chosen residual and parameter set. The most important problem is that the stage is marked exact, but half of the two-engine audit now verifies only a single numerical instance.
