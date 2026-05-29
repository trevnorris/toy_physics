---
stage: 135
reviewer: codex
verdict: PASS
findings_count: 0
files_reviewed: [moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py, moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.txt, moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.wl, moving_throat_pde_stage135_outlet_consistent_mouth_closure_mathematica_audit.txt, stage_135.md, stage_135.tex]
---
# Codex review — stage 135

## What the edit was supposed to do
The directive identified the old SymPy closure residual assertion as tautological because `Sigma_star` was solved from the same closure equation being asserted. The fix was to keep that residual only as a printed sanity probe, delete its load-bearing `raise`, and add substantive SymPy assertions matching the Mathematica coverage. The new checks had to cover the outlet substitution, the range `0 < S_q(Pi_*) < 1`, numerical anchors for `Sigma_m^*`, `M_s^*`, `M_q^*`, and the mixed-lane correction.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:40 | `M_s + M_q*S_q -> Sigma_m*(4 - S_q)` after independent gain-symbol substitution | Yes; wrong sign or factor in either gain leaves a nonzero residual | Yes; matches gain-pair outlet consistency | PASS |
| 2 | moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:58 | Assert `0 < S_q(Pi_star) < 1` | Yes; wrong `S_q` or `Pi_star` can move it outside the interval | Yes; checks finite susceptibility branch range | PASS |
| 3 | moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:69 | Assert `Sigma_m^*` numeric anchor | Yes; wrong `Pi_star`, `S_q`, or closure relation changes `Sigma_star` | Yes; matches quoted `Sigma_m^*` stage output | PASS |
| 4 | moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:71 | Assert `M_s^*` numeric anchor | Yes; wrong outlet factor or `Sigma_star` changes it | Yes; checks `M_s = 4 Sigma_m` carry-forward value | PASS |
| 5 | moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:73 | Assert `M_q^*` numeric anchor | Yes; wrong sign or `Sigma_star` changes it | Yes; checks `M_q = -Sigma_m` carry-forward value | PASS |
| 6 | moving_throat_pde_stage135_outlet_consistent_mouth_closure_sympy_audit.py:81 | Assert mixed-lane correction `M_q^* S_q(Pi_*)` | Yes; wrong `S_q`, `M_q`, or numeric branch changes it | Yes; checks the finite mouth-profile correction value | PASS |

## Bottom line
PASS. I tried to break the new checks by treating the outlet gain sign/factor, the susceptibility value, the hard-coded branch point, and the numeric targets as independently wrong; each edited assertion would fail under those perturbations rather than reducing to `X - X`. The saved SymPy output corroborates the substitution residual as `0`, the interval check as `True`, the anchored values, and the mixed correction, with the old closure residual left only as a print. The Mathematica file is not an edited target here, but its saved transcript independently shows the same claim coverage with explicit `PASS` lines for the corresponding checks.
