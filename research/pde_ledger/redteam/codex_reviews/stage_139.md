---
stage: 139
reviewer: codex
verdict: FINDINGS
findings_count: 4
files_reviewed: [moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py, moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.txt, moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl, moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.txt, stage_139.md, stage_139.tex]
---
# Codex review — stage 139

## What the edit was supposed to do
The directive says the old audits only asserted `R_q^comp = 1/4`, an algebraic identity that did not exercise Family-1 content or the paper card’s checklist. F1 required adding numerical assertions for `r_F1`, `R_q^nat`, the natural and compensated mouth gains, outlet consistency, and compensated `R_q`. F2 required Mathematica to derive the compensated `g_minus` branch with `Solve` rather than transliterating the SymPy formula. F3/F4 required provenance comments for imported literals and the correct Stage 139 Mathematica banner.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | py:38 / wl:71 | `r_F1` vs boxed literal | Yes | Yes, Family-1 branch value | PASS |
| 2 | py:41 / wl:72 | `R_q^nat` vs boxed literal | Yes | Yes, natural branch ratio | PASS |
| 3 | py:44 / wl:73 | `M_s^nat,*` vs boxed literal | Yes | Yes, natural gain value | PASS |
| 4 | py:45 / wl:74 | `M_q^nat,*` vs boxed literal | Yes | Yes, natural gain value | PASS |
| 5 | py:48 / wl:75 | `M_s^comp,*` vs boxed literal | Yes | Yes, compensated gain value | PASS |
| 6 | py:49 / wl:76 | `M_q^comp,*` vs boxed literal | Yes | Yes, compensated gain value | PASS |
| 7 | py:54 / wl:77 | natural outlet consistency | No | Only tautological outlet form | FINDING R1 |
| 8 | py:55 / wl:78 | compensated outlet consistency | No | Only tautological outlet form | FINDING R1 |
| 9 | py:58 / wl:79 | `R_q^comp = 1/4` | No, built into `g_minus`/Solve condition | Adjacent, not independent | FINDING R2 |
| 10 | wl:49 | `g_minus` Solve result vs closed form | Yes | Partially, branch selection only | PASS |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:54-55`; `moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:77-78`
- **What's wrong:** `assert abs(Pi_star - (Ms_nat + Mq_nat * Sq_star)) < tol_algebraic` and the analogous checks use `Ms = Pi/(1 - Rq*Sq)` and `Mq = -Rq*Ms` from the same script. The residual reduces to `Pi_star - Pi_star` by construction for both branches.
- **What it should be:** Derive or import the gain pair independently from the outlet equation, then check `Pi_* = M_s + M_q S_q(Pi_*)`.

### R2 — tautological_check
- **Where:** `moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:15-16,58`; `moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:42,50,79`
- **What's wrong:** SymPy defines `g_minus = rF - sqrt(1 + rF**2)/2`, making `R_q^comp = 1/4` algebraically for any `rF`. Mathematica derives `gMinus` by solving the same condition `(gc - rF)^2 == (1 + rF^2)/4`, then checks `R_q^comp - 1/4`.
- **What it should be:** Derive `g_c` from an independent compensated-branch condition, then compute and test `R_q`.

### R3 — paper_misalignment
- **Where:** `stage_139.tex:Checks`; `moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:8-18`; `moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:31-52`
- **What's wrong:** The paper card requires checking the self-matched susceptibility closure before using the one-scalar branch law. The scripts only hardcode `S_q(Pi_*)` and never evaluate any susceptibility-closure residual.
- **What it should be:** Add an independent check of the susceptibility closure at `Pi_star`, or explicitly verify the imported upstream audit result rather than only copying the literal.

### R4 — other
- **Where:** `moving_throat_pde_stage139_family1_actual_mouth_gains_sympy_audit.py:5-8`; `moving_throat_pde_stage139_family1_actual_mouth_gains_mathematica_audit.wl:28-32`
- **What's wrong:** The directive required provenance comments saying `r_F1` came from Stage 223 and `Pi_*`, `S_q(Pi_*)` came from Stage 236. The edited scripts instead say Stage 121 and Stage 134, with no applied/rework note explaining the deviation.
- **What it should be:** Use the correct provenance, or document the orchestrator rework that justifies different stage numbers.

## Bottom line
The saved outputs show the new checks ran and passed, including `all assertions passed` and the Mathematica `PASS:` lines, but the most important added paper-facing checks are not trustworthy. I tried to break them from the source by treating `Rq`/`Sq` as wrong inputs: the outlet residuals still collapse to zero because the gains are defined from that same equation, and the compensated `R_q` check is built into the selected `g_minus` condition.
