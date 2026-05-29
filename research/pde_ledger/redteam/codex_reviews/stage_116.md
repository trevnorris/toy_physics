---
stage: 116
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [audit_directive_116.md, verify_directive_116.md, stage_116.tex, moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py, moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.txt, moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.wl, moving_throat_pde_stage116_dn_mixed_tube_realization_mathematica_audit.txt]
---
# Codex review — stage 116

## What the edit was supposed to do
The directive required replacing hardcoded `kappa0_from_tube = 4 L_W^2/(pi^2 a^2)` with a D/N eigenvalue derivation from `q''+k^2 q=0`, `q(0)=0`, `q'(L_W)=0`, giving `k_W=pi/(2L_W)`. It also required re-anchoring `kappa_c` through the geometric tube expression evaluated at `L_W_required`, extracting `gamma0` from `D_bare`, and adding the missing SymPy assertion for the paper's tube-length law. The listed `redteam/directives/stage_116.md` path was absent; I used the matching stage-116 directive under `redteam/tmp_prompts`.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | sympy:31-47, wl:40-52 | D/N trial, BCs, and `kappa0_derived` | yes | yes, D/N half-wave normalization | PASS |
| 2 | sympy:59-62, wl:60 | tube-length law | yes | yes, card quote `L_W=pi a sqrt((1+r_c)/3)/2` | PASS |
| 3 | sympy:65-69, wl:62-64 | geometric kappa round-trip at `L_W_required` | no, solved from same target | partly | FINDING |
| 4 | sympy:73-74, wl:69-70 | final `gamma_c - 1/9` | no, `gamma0_bare` is literal | yes, but tautologically | FINDING |
| 5 | sympy:80-83,89-92, wl:73-80 | extract `gamma0`/`kappa0` from `D_bare` | yes, sign/coefficient errors fail | yes, bare outgoing scale | PASS |
| 6 | sympy:84-88, wl:76-77 | scaled canonical branch renormalizes to canonical | no, constructed as `(1+r_c)*canonical` | adjacent only | FINDING |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:65`
- **What's wrong:** `kappa0_bare_geom = sp.simplify(4 * L_W_required**2 / (sp.pi**2 * a**2))` is checked against `(1+r_c)/3` after `L_W_required` was defined by solving `kappa0_from_tube == (1+r_c)/3`. The Mathematica lines 62-64 mirror this. This is a round-trip through the same equation, so a wrong geometric prefactor would still satisfy this check once `L_W_required` is solved from that wrong prefactor.
- **What it should be:** Count the eigenvalue-derived `kappa0_derived` check and the independent tube-length-law check as the substantive assertions; do not treat this solve/substitute round-trip as independent evidence.

### R2 — tautological_check
- **Where:** `moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:70`
- **What's wrong:** `gamma0_bare = sp.simplify((1 + r_c) / 9)` followed by `gamma_c = gamma0_bare / (1+r_c)` and `final gamma_c - 1/9` proves only `((1+r_c)/9)/(1+r_c)=1/9`. The later coefficient extraction is substantive, but this final `gamma_c` assertion does not use it.
- **What it should be:** Define `gamma_c` from the extracted coefficient, e.g. `gamma_c_from_D = gamma0_from_D/(1+r_c)`, then assert `gamma_c_from_D - 1/9 == 0`.

### R3 — tautological_check
- **Where:** `moving_throat_pde_stage116_dn_mixed_tube_realization_sympy_audit.py:84`
- **What's wrong:** `D_bare = (1+r_c)*(1 - z**2/3 - I*z**5/9)` and `D_final = D_bare/(1+r_c)` are built from the exact expected canonical branch, so `D_final - canonical` is `X - X` by construction. Mathematica lines 72 and 76-77 repeat the same construction.
- **What it should be:** Either derive `D_bare` from an upstream expression and factor it, or rely on coefficient-level checks (`z^2`, `z^5`) and stop presenting this quotient as an independent assertion.

## Bottom line
The D/N eigenvalue block and the explicit tube-length-law assertion are real checks, and the saved outputs corroborate them with zero residuals/PASS lines. But several edited assertions still have no physics-level failure mode: the kappa round-trip solves from the same target it rechecks, `gamma_c` is still proved from a literal `gamma0_bare`, and the canonical renormalization check divides a deliberately scaled canonical polynomial by its own scale. The most important problem is that these tautological checks remain mixed into the proof chain, so a PASS would overstate what the stage now verifies.
