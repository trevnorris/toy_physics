---
stage: 170
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [stage_170.md, stage_170.tex, moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py, moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.txt, moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl, moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.txt]
---
# Codex review — stage 170

## What the edit was supposed to do
The directive required adding the missing weak-axisymmetric signature check to both engines after the user chose direction (a). The intended check was to feed lane-scaled grouped defects \(\delta D_{A,n}=\epsilon\lambda_A D_n^{(1)}\), \(\delta N_{A,0}=\epsilon\lambda_A N_0^{(1)}\), with \(\lambda=(1,1/2,-1)\), through the verified outlet maps and verify the \(\delta\kappa_W,\delta\gamma_W\) signature plus the scalar amplitudes \(\kappa_1,\gamma_1\). The directive also required the Mathematica engine to stop line-by-line mirroring SymPy for the earlier linearization/inversion path, and required both banners to say Stage 170.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:160 | \(\delta\kappa_W^{(A)}=\epsilon\lambda_A\kappa_1\) | No, not against a wrong Sec. 5 map/amplitude copied into both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 2 | moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:161 | \(\delta\gamma_W^{(A)}=\epsilon\lambda_A\gamma_1\) | No, not against a wrong Sec. 5 map/amplitude copied into both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 3 | moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:163 | \(\kappa\) signature \(21=\frac12 20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 4 | moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:164 | \(\kappa\) signature \(22=-20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 5 | moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:165 | \(\gamma\) signature \(21=\frac12 20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 6 | moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:166 | \(\gamma\) signature \(22=-20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 7 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:150 | \(\delta\kappa_W^{20}=\epsilon\kappa_1\) | No, same formula is used to define both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 8 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:151 | \(\delta\kappa_W^{21}=\epsilon\kappa_1/2\) | No, same formula is used to define both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 9 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:152 | \(\delta\kappa_W^{22}=-\epsilon\kappa_1\) | No, same formula is used to define both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 10 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:153 | \(\delta\gamma_W^{20}=\epsilon\gamma_1\) | No, same formula is used to define both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 11 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:154 | \(\delta\gamma_W^{21}=\epsilon\gamma_1/2\) | No, same formula is used to define both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 12 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:155 | \(\delta\gamma_W^{22}=-\epsilon\gamma_1\) | No, same formula is used to define both sides | Only tests homogeneity of a re-hardcoded target map | FINDING |
| 13 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:156 | \(\kappa\) signature \(21=\frac12 20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 14 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:157 | \(\kappa\) signature \(22=-20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 15 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:158 | \(\gamma\) signature \(21=\frac12 20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |
| 16 | moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:159 | \(\gamma\) signature \(22=-20\) | Only if the lane constants are mistyped | Partially, but only after hardcoding a linear target map | WEAK |

## Findings (only if verdict = FINDINGS)
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage170_linear_grouped_outlet_map_sympy_audit.py:144-161`; `moving_throat_pde_stage170_linear_grouped_outlet_map_mathematica_audit.wl:138-155`
- **What's wrong:** The new Sec. 5 block defines `kappa_map`/`gamma_map` directly as the claimed closed-form outlet maps, then defines `kappa1`/`gamma1` using the same closed forms, then checks lane-scaled outputs against those amplitudes. For example, `kappa_map(dD2_, dD0_)` is `3*(1 - sigma)*(dD2_ + dD0_/9)/(sigma*D0)`, and `kappa1` is the same expression with `D2_1,D0_1`; therefore `kappa_map(eps_l*lam*D2_1, eps_l*lam*D0_1) - eps_l*lam*kappa1` vanishes by construction. I tried to break the check by imagining the Sec. 5 target map coefficient or sign being wrong but copied into both `kappa_map` and `kappa1`; the added checks would still pass, while the paper-side amplitude would be wrong. The saved outputs do show all these residuals as `0`/`PASS`, but that corroborates the tautology, not the paper claim.
- **What it should be:** The Sec. 5 lane outputs should be obtained by substituting the lane defects into the already-derived map objects, e.g. SymPy `dkappa_from_du2.subs({dD2: eps_l*lam*D2_1, dD0: eps_l*lam*D0_1})` and `dgamma_from_dP0.subs({dN0: eps_l*lam*N0_1, dD0: eps_l*lam*D0_1})`, then compared to independently written target amplitudes `kappa1` and `gamma1`. The Mathematica block should do the analogous substitution into `dkappaFromdu2` and `dgammaFromdP0`.

## Bottom line
The Mathematica F2 rework is meaningfully different from the SymPy path for the earlier linearization/inversion, and both transcripts show the Stage 170 banner and all residuals passing. The added weak-axisymmetric Sec. 5 checks, however, do not actually feed the lane defects through the previously derived outlet maps; they re-hardcode the target map and the target amplitude from the same formula. The most important problem is that a wrong Sec. 5 amplitude copied into both `kappa_map`/`gamma_map` and `kappa1`/`gamma1` would still print `0`, so the new assertion is not a trustworthy non-tautological verification of the paper-card deliverable.
