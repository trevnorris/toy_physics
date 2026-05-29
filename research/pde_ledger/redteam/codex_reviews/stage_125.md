---
stage: 125
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [audit_directive_125.md, verify_directive_125.md, stage_125.tex, moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py, moving_throat_pde_stage125_positive_source_theorem_sympy_audit.txt, moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl, moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.txt]
---
# Codex review — stage 125

## What the edit was supposed to do
The listed directive file `redteam/directives/stage_125.md` was absent in this checkout, so I used the available directive copies under `redteam/tmp_prompts`. F1 required adding an integral-side check for positive normalized mouth sources, with the selected repair using the beta-like family `sigma_a(z) = (a+1)(z/L)^a/L`. F2 required making the Mathematica branch calculation independent by deriving `g_-` and `g_+` with `Solve`, then checking those roots against the SymPy closed forms while removing the old balance-relation tautologies.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:74 | `int sigma_a dz = 1` | yes | positive-normalized source family | PASS |
| 2 | moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:79 | uniform moment `g_a(0)=2/pi` | yes | concrete cosine-moment representation | PASS |
| 3 | moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:85 | peaked proxy `abs(g_a(100)) < 0.05` | yes, but too weak | only proximity to zero, not `0 <= g <= 1` | FINDING |
| 4 | moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:87 | uniform lower bound | yes | lower bound for one positive source | PASS |
| 5 | moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:88 | uniform upper bound | yes | upper bound for one positive source | PASS |
| 6 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:57 | derive branches with `Solve` | yes | independent branch construction | PASS |
| 7 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:63 | `g_-` closed-form match | yes | lower branch value | PASS |
| 8 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:65 | `g_+` closed-form match | yes | upper branch value | PASS |
| 9 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:86 | `int sigma_a dz = 1` | yes | positive-normalized source family | PASS |
| 10 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:92 | uniform moment `g_a(0)=2/Pi` | yes | concrete cosine-moment representation | PASS |
| 11 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:93 | peaked limit `a -> Infinity` is `0` | yes | endpoint moment | PASS |
| 12 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:95 | uniform lower bound | yes | lower bound for one positive source | PASS |
| 13 | moving_throat_pde_stage125_positive_source_theorem_mathematica_audit.wl:96 | uniform upper bound | yes | upper bound for one positive source | PASS |

## Findings (only if verdict = FINDINGS)
### R1 — insufficient_verification
- **Where:** `moving_throat_pde_stage125_positive_source_theorem_sympy_audit.py:85`
- **What's wrong:** `expect_true("g[peaked@L proxy a=100] < 0.05", bool(abs(g_a_large) < sp.Rational(1, 20)))` does not assert the paper-side bound `0 <= g[sigma] <= 1`. It would pass for a positive-source moment of `-0.01`, which violates the lower-bound half of the theorem. The saved SymPy output shows the current value is positive (`0.0153964...`) and the check passes, but the assertion is shaped as a smallness test rather than a range test.
- **What it should be:** At minimum, assert the bounded proxy directly, e.g. `0 <= g_a_large < 1/20`, and preferably add an intermediate finite-`a` range check or a direct pointwise-to-integral lemma tying the kernel bound to `0 <= g[sigma] <= 1`.

## Bottom line
The Mathematica F2 repair is non-tautological: `g_±` are now obtained by `Solve` and checked against closed forms, with output showing both closed-form residuals are `0`. The integral-family F1 repair is mostly real, and the saved outputs corroborate the new normalization and uniform-moment checks. The blocker is the SymPy peaked-source assertion: it accepts small negative moments, so it can pass while the paper’s lower-bound claim is false for that positive normalized source.
