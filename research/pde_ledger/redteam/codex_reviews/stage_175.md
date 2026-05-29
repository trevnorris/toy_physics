---
stage: 175
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py, moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.txt, moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl, moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.txt, stage_175.md, stage_175.tex]
---
# Codex review — stage 175

## What the edit was supposed to do
F1 was supposed to remove the tautological `Sigma_N - (2 dln Lambda - dK)` check; the orchestrator ultimately chose the directive-approved fallback of deleting the redundant check and keeping only `Sigma_N - dln(Lambda^2/K)`. F2 was supposed to add an aggregate `Xi_load` check using normalized weights, so the no-go `Xi_load = -dK` is asserted rather than only printed. F3 was supposed to fix the stale stage banner and make the Mathematica differential block structurally independent, unless blocked; the final orchestrator note accepted the block and left Mathematica as a policy mirror.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:135 | `Sigma_N - dln(Lambda^2/K)` | yes; wrong `-kappa` sign/omission fails | yes | PASS |
| 2 | moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:98 | Mathematica mirror of `Sigma_N - dln(Lambda^2/K)` | locally yes, but not independently | yes, but copied choreography | FINDING |
| 3 | moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:162 | `Xi_load (all shapes frozen) + dK` | yes; missing normalization or wrong `Sigma_N_common` fails | yes | PASS |
| 4 | moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:118 | Mathematica `Xi_load` aggregate | yes | yes | PASS |

## Findings
### R1 — transliteration
- **Where:** `moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:26`
- **What's wrong:** `dlog[expr_] := FullSimplify[D[Log[FullSimplify[expr, Assumptions -> $Assumptions]], eps] /. eps -> 0, Assumptions -> $Assumptions];` The Sigma_N block at lines 95-98 still mirrors the SymPy `dlog(diff(log(...)))` route and expression choreography. The saved output corroborates that it passes, but it is not an independent second-engine derivation; a wrong target copied from the SymPy script would still be likely to pass in both.
- **What it should be:** Use the directive’s independent Mathematica route for the Sigma_N block, e.g. a `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps]` construction, and compare that independently derived slope against the same paper-side target.

## Bottom line
The SymPy F1 fallback and the new weighted `Xi_load` check are not `X - X` tautologies, and the saved transcripts show the expected `= 0` / `PASS` lines. The blocking issue is that the Mathematica audit remains a line-by-line mirror for the differential block after F3-step3 was waived, so the second engine does not provide independent protection against a copied wrong target.
