---
stage: 166
reviewer: codex
verdict: FINDINGS
findings_count: 1
files_reviewed: [moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.py, moving_throat_pde_stage166_bundle_inversion_four_drifts_sympy_audit.txt, moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl, moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.txt, stage_166.md, stage_166.tex]
---
# Codex review — stage 166

## What the edit was supposed to do
The directive required explicit assertions for the general inversion formulas \( \delta\ln\rho_w=\frac12\delta\ln\Theta_w \), \( \delta\ln a=\frac12\delta\ln K_s-\frac14\delta\ln\Theta_w \), \( \delta\ln c_s=\frac12\delta\ln K_s-\frac14\delta\ln\Theta_w+\frac15\delta\ln P_0 \), and \( \delta\ln Z_q=\delta\ln K_q-\frac25\delta\ln P_0 \), because the old checks only printed these formulas or checked frozen/forward substitutions. It also required fixing the wrong Stage 149 banner to Stage 166 and adding a Mathematica matrix-inverse route so the second engine was not just a Solve transliteration. The paper card’s output is the inversion of \((\Theta_w,K_s,K_q,P_0)\) into the four remaining microscopic drifts.

## Edited-assertion verdict table
| # | file:line | edited check (short) | can it fail? | exercises paper claim? | verdict |
|---|-----------|----------------------|--------------|------------------------|---------|
| 1 | sympy audit:53 | `drho general` vs `dTheta/2` | Yes; changing the `dTheta=2 drho` law or target slope breaks it. | Yes, M1. | PASS |
| 2 | sympy audit:54 | `da general` vs `dKs/2 - dTheta/4` | Yes; wrong `Ks` or `Theta` slope breaks it. | Yes, M2. | PASS |
| 3 | sympy audit:55 | `dcs general` vs `dKs/2 - dTheta/4 + dP/5` | Yes; wrong `P0` or shared drift slope breaks it. | Yes, M3. | PASS |
| 4 | sympy audit:56 | `dZ general` vs `dKq - 2*dP/5` | Yes; wrong `Kq`/`P0` coupling breaks it. | Yes, M4. | PASS |
| 5 | mathematica audit:48 | `drho general` vs `dTheta/2` | Yes; it compares Solve output to an external target. | Yes, M1. | PASS |
| 6 | mathematica audit:49 | `da general` vs `dKs/2 - dTheta/4` | Yes; it is not built as `daSol-daSol`. | Yes, M2. | PASS |
| 7 | mathematica audit:50 | `dcs general` vs `dKs/2 - dTheta/4 + dP/5` | Yes; coefficient/sign errors survive to the residual. | Yes, M3. | PASS |
| 8 | mathematica audit:51 | `dZ general` vs `dKq - 2*dP/5` | Yes; a wrong `dP` coefficient fails. | Yes, M4. | PASS |
| 9 | mathematica audit:69 | matrix `drho` target | Yes; a wrong first matrix row or target slope fails. | Yes, independent route to M1. | PASS |
| 10 | mathematica audit:70 | matrix `da` target | Yes; a wrong `Ks`/`Theta` row or target slope fails. | Yes, independent route to M2. | PASS |
| 11 | mathematica audit:71 | matrix `dcs` target | Yes; a wrong `P0` coupling fails. | Yes, independent route to M3. | PASS |
| 12 | mathematica audit:72 | matrix `dZ` target | Yes; a wrong `Kq`/`P0` coupling fails. | Yes, independent route to M4. | PASS |
| 13 | mathematica audit:76 | matrix round-trip | No; `solVec` is defined as `Inverse[Mmat] . observables`, so this is `Mmat.Inverse[Mmat].v - v`. | No; it only checks internal inverse consistency. | FINDING |

## Findings
### R1 — tautological_check
- **Where:** `moving_throat_pde_stage166_bundle_inversion_four_drifts_mathematica_audit.wl:76`
- **What's wrong:** `expectZero["matrix round-trip", Total[(Mmat . solVec - {dTheta, dKs, dKq, dP})^2]];` is tautological because `solVec` was defined at line 68 as `Inverse[Mmat] . {dTheta, dKs, dKq, dP}`. For any invertible but incorrectly transcribed `Mmat`, this round-trip remains zero; changing a forward-map coefficient in `Mmat` does not break this assertion as long as the matrix stays invertible. It therefore does not do what the comment claims: confirm `Mmat` was transcribed correctly from `eq1`–`eq4`.
- **What it should be:** Replace or supplement it with a real transcription check, e.g. compare `Mmat . {drho, da, dcs, dZ}` against the independently written forward laws `{2*drho, drho + 2*da, dZ + 2*dcs - 2*da, 5*(dcs - da)}`.

## Bottom line
The new general inversion assertions are meaningful and the saved outputs corroborate them: SymPy prints the four `= 0` lines, and Mathematica prints matching `PASS` lines plus the matrix target `PASS` lines. The blocker is the added Mathematica round-trip assertion: I tried to break it by perturbing the forward matrix in thought, and it still passes because it inverts whatever matrix it was given. Since an edited assertion is tautological, this stage cannot receive a trustworthy PASS.
