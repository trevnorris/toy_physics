---
stage: 142
reviewer: codex
verdict: FINDINGS
findings_count: 2
files_reviewed: [moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py, moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.txt, moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl, moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.txt, stage_142.md, stage_142.tex]
---
# Codex review — stage 142

## What the edit was supposed to do
The directive required F1-F4 to be applied: add a non-tautological canonical-point check beyond the algebraic \(R_q(g_-)=1/4\), anchor five recorded numerical targets, add Mathematica-side checks to reduce transliteration risk, and add provenance comments for \(r_{F1}\) and \(S_q\). F5 was a held paper-misalignment item and was not supposed to be changed by that directive. The paper card requires this stage to audit the coupled mouth fixed point/gain-selection ledger, especially \(\Pi=\Sigma_0[1-R_q(\Pi)\mathcal S_q(\Pi)]\) and the recorded numerical fixed point.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:75-78; moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:109 | \(R_q(\Pi_*)=1/4\) after solving \(g_\Pi(\Pi_*)=g_-\) | Only on solver/numeric residual; wrong self-consistent \(r\) or \(g_\Pi\) still passes | No, it follows from the same definitions and solve equation | FINDING |
| 2 | moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:80; moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:110 | \(g_-^{F1}\) numeric target | Yes, against external decimal target | Yes | PASS |
| 3 | moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:81; moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:111 | \(\Pi_*\) numeric target | Yes, against external decimal target | Yes | PASS |
| 4 | moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:82-84; moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:112-114 | \(S_q(\Pi_*)\), \(\Sigma_0(\Pi_*)\), \(\widehat T(\Pi_*)\) targets | Yes, against external decimal targets | Yes | PASS |
| 5 | moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:64-73 | \(g_\Pi\) closed form compared to `Series[gPi]` | Not for copied analytic formula errors; it compares expression to its own series | No, it checks local Taylor truncation only | FINDING |
| 6 | moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:77-79 | \(100\pi^2(1+r^2)=4107\) | Yes for many wrong \(r\) magnitudes, though sign-blind | Adjacent upstream anchor | PASS |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage142_selfconsistent_mouth_branch_sympy_audit.py:75-78`; `moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:109`
- **What's wrong:** `Rq_star_residual = abs(float(Rq_star - sp.Rational(1,4)))` and `expectApprox["R_q(Pi_*) numeric = 1/4", rQStar, 1/4, 10^-20];` are determined by the same construction: \(\Pi_*\) is solved from \(g_\Pi(\Pi_*)=g_-\), while \(R_q=(g_\Pi-r)^2/(1+r^2)\) and \(g_-=r-\sqrt{1+r^2}/2\). For any self-consistent but wrong \(r\) or copied wrong \(g_\Pi\), this still gives \(R_q(\Pi_*)=1/4\). The SymPy output also shows a nonzero residual `1.945...E-18` and the code uses `1e-15`, not the directive’s `1e-20`.
- **What it should be:** Treat this only as a redundant solver-consistency check, not as the non-tautological canonical anchor. The real non-tautological checks are the external decimal targets, or an \(R_q\) check evaluated at an independently supplied \(\Pi_*\) target rather than at a point solved from the same equations.

### R2 — transliteration
- **Where:** `moving_throat_pde_stage142_selfconsistent_mouth_branch_mathematica_audit.wl:64-73`
- **What's wrong:** `gPiSeries = Normal[Series[gPi, {piM, 0, 4}]];` followed by comparing `gPi` to `gPiSeries` is not an independent derivation. A typo copied into `gPi` is also copied into its Taylor series, so the residuals in the saved output only show that an expression is close to its own truncated series at `0.1,0.2,0.3`.
- **What it should be:** Compare `Series[gPi]` coefficients against independently entered expected coefficients, or encode \(g_\Pi\) a second way from the upstream mouth-source law and compare those two independent expressions.

## Bottom line
The numerical target checks added for F2 are useful and corroborated by the saved outputs, but the two checks that were supposed to remove tautology/transliteration risk do not do that. I tried to break them conceptually by perturbing \(r_{F1}\) or replacing \(g_\Pi\) with a self-consistent copied typo: the \(R_q(\Pi_*)\) check still follows from the solve equation, and the Mathematica series check still compares the typo to its own series. The most important problem is that the added “independent” evidence is still built from the same intermediates it claims to verify.
