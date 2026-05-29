---
stage: 121
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [audit_directive_121.md, stage_121.tex, moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py, moving_throat_pde_stage121_geometric_r_selection_sympy_audit.txt, moving_throat_pde_stage121_geometric_r_selection.md]
---
# Codex review — stage 121

## What the edit was supposed to do
The directive required no script change for F1 after the notes-side `168π²` typo was resolved to `100π²`. For F2, it required four new `expect_zero` anchors: the explicit squared branch law, the `L/a=37/20` symbolic value, the squared `r_c(F1)` value, and the existence-threshold root. For F3, it required an explicit `Omega_W` substitution check under `L_W=L`.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:25 | `r_geom^2 = 12 L^2/(pi^2 a^2)-1` | yes: wrong tube coefficient or inversion leaves nonzero residual | yes: boxed `r_geom(L/a)` law | PASS |
| 2 | moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:41 | `r_F1 = sqrt(4107-100*pi^2)/(10*pi)` | yes: wrong aspect ratio, sign, or coefficient fails | yes: preferred `L/a=37/20` value in stage notes/directive | PASS |
| 3 | moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:51 | `r_c(F1)=4107/(100*pi^2)-1` | yes: wrong square or target coefficient fails | yes: squared preferred-branch value | PASS |
| 4 | moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:61 | `Omega_W(L_W=L)=pi*c_s/(2L)` | yes: wrong source pole or target denominator fails | yes: boxed mixed-tube pole consequence | PASS |
| 5 | moving_throat_pde_stage121_geometric_r_selection_sympy_audit.py:69 | `r_geom` vanishes at `pi/(2 sqrt(3))` | yes: wrong threshold or branch law fails | yes: existence-boundary claim | PASS |

## Bottom line
PASS. I tried to break the checks by looking for `X-X` construction, copied targets from the same intermediate, sign loss from the squared `r_geom` assertion, unjustified positivity assumptions, and output omissions. The squared branch-law check would miss a global negative branch by itself, but the `r_F1` assertion against the positive target would fail under that branch; `L,a,L_W,r,c_s` positivity is consistent with the geometric setup. The saved output corroborates each edited assertion with `= 0` on lines 6, 10, 13, 15, and 17, and there is no Mathematica engine to compare because the stage card explicitly says none.
