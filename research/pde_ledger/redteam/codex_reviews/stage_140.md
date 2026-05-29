---
stage: 140
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [stage_140.md, stage_140.tex, moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py, moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.txt, moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl, moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.txt]
---
# Codex review — stage 140

## What the edit was supposed to do
F1 required correcting the Mathematica transcript banner from Stage 123 to Stage 140, without changing the underlying math. F2 required turning the previously print-only Section-3 numerical values into enforced checks: `That_nat`, `That_comp`, and the fractional traction enhancement must match the notes' boxed numerics to `1e-12`. The intent was to make regressions in the `M_s` literals or the `sqrt(9 M_s/20)` conversion fail instead of silently printing changed values.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:26 | banner reports Stage 140 | Yes; stale Stage 123 would show in saved transcript | Yes; corroborates correct stage transcript identity | PASS |
| 2 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py:25 | `That_nat` matches `0.866512630228382` | Yes; changing `Ms_nat` or the `9/20` formula changes the computed value | Yes; checks the recorded numerical fixed point implied by `M_s=(20/9) That^2` | PASS |
| 3 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py:26 | `That_comp` matches `0.901484054174206` | Yes; changing `Ms_comp` or the `9/20` formula changes the computed value | Yes; checks the compensated numerical fixed point | PASS |
| 4 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_sympy_audit.py:27 | fractional enhancement matches `0.0403588161624` | Yes; depends on both computed hats and the ratio law | Yes; checks the printed gain enhancement value | PASS |
| 5 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:55-60 | `That_nat` numeric PASS gate | Yes; nonzero diff beyond tolerance calls `fail` and exits | Yes; same numerical fixed-point deliverable as SymPy | PASS |
| 6 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:55-61 | `That_comp` numeric PASS gate | Yes; nonzero diff beyond tolerance calls `fail` and exits | Yes; same compensated fixed-point deliverable as SymPy | PASS |
| 7 | moving_throat_pde_stage140_selfmatched_mouth_susceptibility_mathematica_audit.wl:55-62 | fractional enhancement numeric PASS gate | Yes; nonzero diff beyond tolerance calls `fail` and exits | Yes; same enhancement deliverable as SymPy | PASS |

## Bottom line
PASS. I checked for `X - X` construction by tracing each edited numeric check: the computed quantities come from the `M_s` literals and `sqrt(9 M_s/20)`, while the targets are separately typed note decimals, so the assertions are not tautological. Perturbing either `M_s` literal, the `9/20` factor, or the target decimal would break the corresponding check beyond the stated tolerance. The saved outputs corroborate the edits: SymPy reaches `numeric fixed-point checks PASS`, Mathematica prints the Stage 140 banner and all three new numeric `PASS:` lines, while the existing susceptibility residual remains `0`.
