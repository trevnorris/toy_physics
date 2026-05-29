---
stage: 169
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py, moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.txt, moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl, moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.txt, stage_169.md, stage_169.tex]
---
# Codex review — stage 169

## What the edit was supposed to do
The directive identified the old `eps_perp - Xi_perp*Igrp` assertion as a harmless distributive-law tautology and required new per-coefficient numeric checks of `Xi_perp` against the paper's Family-1 weights. It also accepted the Mathematica file as a documented policy mirror rather than forcing an independent re-derivation. Finally, it required both script banners and saved transcripts to say Stage 169 instead of the stale Stage 152 label.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:133 | extract XiT/Xiv/XiL coefficients from `Xi_perp` after `g_num,r_num` substitution | yes; any wrong coefficient or wrong imported numeric branch exceeding `1e-12` trips the loop | yes; checks the stated quadratic grouped-transport weights carried into `eps_perp` | PASS |
| 2 | moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:136 | compare computed coefficients to fixed paper literals | yes; `got` is computed from `Xi_perp`, `want` is an external literal, not `got` reused | yes; anchors the printed Family-1 combination to the paper-side numbers | PASS |
| 3 | moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:105 | `Coefficient[xiPerp /. {g -> gNum, r -> rNum}, ...]` checks | yes; copied wrong weights would fail against the fixed paper literals | yes, though the engine is still a policy mirror rather than an independent derivation | PASS |
| 4 | moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:108 | tolerance loop with `fail/pass` for the three coefficients | yes; `Abs[got - want] > 10^-12` calls `fail` | yes; saved output has PASS lines for xiT, xiv, xiL | PASS |
| 5 | moving_throat_pde_stage169_no_linear_p2_scalar_slippage_sympy_audit.py:30 / moving_throat_pde_stage169_no_linear_p2_scalar_slippage_mathematica_audit.wl:31 | Stage banner corrected to 169 | yes as an output-label check; stale text would remain visible in transcripts | housekeeping, not mathematical | PASS |

## Bottom line
PASS. I tried to break the new checks by treating the old `eps_perp - Xi_perp*Igrp` line as the load-bearing proof; it is still tautological, but the directive explicitly left it in place and the new coefficient comparisons no longer rely on it. I checked that the coefficient values are extracted from the symbolic `Xi_perp` and compared to fixed paper literals, so changing `g`, `r`, or any of the three weighting formulas would exceed the tolerance. The SymPy transcript prints the three coefficient residuals, including nonzero truncation residuals for Xiv and XiL below `1e-12`; the Mathematica transcript prints corresponding PASS lines. The `.wl` remains a mirror, but that was the directive's explicit F2 policy disposition, not an unreviewed mathematical assertion.
