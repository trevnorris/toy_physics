---
stage: 137
reviewer: codex
verdict: FINDINGS
findings_count: 3
files_reviewed: [moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py, moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.txt, moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl, moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.txt, stage_137.md, stage_137.tex]
---
# Codex review — stage 137

## What the edit was supposed to do
The directive required new non-vacuous anchors for the explicit gain pair \(M_s,M_q\), a Schur/static-limit check for the imported core susceptibility form, and an outlet-consistency reduction. It also required the Mathematica Schur check to use `Series` while SymPy used `sp.limit`. F3 separately required a matrix-Schur reconstruction using `M_core` / `Inverse` to address the hardcoded \(\rho_c,\sigma_c\); F4 was a paper-misalignment hold on the stage-number banner.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:26 | `M_s` equals paper closed form | Yes, for sign/factor errors in `rho_c` or `Ms` propagation | Yes, matches quoted gain | PASS |
| 2 | moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:27 | `M_q` equals paper closed form | Yes, for sign/factor errors in `sigma_c` or `Mq` propagation | Yes, matches quoted gain | PASS |
| 3 | moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:48-49 | `mS,mQ` equal paper closed forms | Yes, same algebraic anchor in second CAS | Yes, matches quoted gain | PASS |
| 4 | moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:37 | static limit equals `rho_c - sigma_c` | Not for wrong `rho_c`/`sigma_c`; both sides move together | Only adjacent to the closure claim | FINDING |
| 5 | moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:55 | Series static term equals `rho_c - sigma_c` | Not for wrong `rhoC`/`sigmaC`; both sides move together | Only adjacent to the closure claim | FINDING |
| 6 | moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:45 | outlet residual at `S_q = 0`, `Pi = M_s` | No for `M_q` errors; the mixed channel is removed | No, does not test gain-pair outlet consistency | FINDING |
| 7 | moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:58-59 | outlet residual at `sqVar = 0`, `piVar = mS` | No for `mQ` errors; the mixed channel is removed | No, does not test gain-pair outlet consistency | FINDING |
| 8 | both scripts | F3 matrix-Schur reconstruction | No check present | No, hardcoded-result issue remains | FINDING |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:37`; `moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:55`
- **What's wrong:** `assert sp.simplify(static_limit - (rho_c - sigma_c)) == 0` and `expectZero["Schur static limit equals rho_c - sigma_c", staticLimit - (rhoC - sigmaC)];` compare a limit/series built from the same `rho_c,sigma_c` symbols against those same symbols. If `sigma_c` is hardcoded with a wrong factor, both sides change together and the check still prints/passses.
- **What it should be:** Build `rho_c,sigma_c` from an independent core Schur object or compare the full susceptibility expression against an independently constructed source formula, not against its own static symbols.

### R2 — tautological_check
- **Where:** `moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:45`; `moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:58-59`
- **What's wrong:** `assert sp.simplify(outlet_residual.subs(Sq_var, 0).subs(Pi_var, Ms)) == 0` and the Mathematica equivalent set `S_q=0`, which deletes the entire `M_q` term. A flipped sign, missing denominator, or even `M_q = 0` would still pass, contrary to the directive's claim that this catches `M_q` sign convention errors.
- **What it should be:** Check outlet consistency with a nonzero independently defined `S_q(Pi)` branch or a recorded fixed point; do not use a substitution that eliminates the mixed-channel gain.

### R3 — insufficient_verification
- **Where:** `moving_throat_pde_stage137_core_to_mouth_gain_map_sympy_audit.py:9-10`; `moving_throat_pde_stage137_core_to_mouth_gain_map_mathematica_audit.wl:34-35`
- **What's wrong:** `rho_c` and `sigma_c` remain directly assigned, and there is no `M_core`, `sp.Matrix`, `Inverse[...]`, `rho_c_schur`, or `sigma_c_schur` block in either script. The saved outputs also lack the required F3 reproduction line/PASS, so the hardcoded-result finding was not applied.
- **What it should be:** Add the requested independent matrix-Schur reconstruction, or explicitly block F3; as written, the direct gain checks only prove consistency with hardcoded inputs.

## Bottom line
The saved outputs corroborate that the added checks execute, but the important new checks are not adversarial. By inspection, changing the hardcoded `sigma_c` factor would not break the Schur static-limit checks, and flipping or deleting `M_q` would not break the outlet checks because `S_q=0` removes that term. The F3 matrix-Schur derivation is absent entirely, so the audit still does not independently verify the core residues behind the paper's gain map.
