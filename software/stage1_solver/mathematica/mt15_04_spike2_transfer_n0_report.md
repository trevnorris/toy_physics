# MT15-04 Spike-2 Transfer N0 Report

Status: Spike-2 only. No V2-22B packet, no `R_norm`, no V2-22B extension, and no `S_eta^(A)` return source.

**MACHINERY LABEL:** the coefficients below are method-demonstration numbers on the non-self-consistent M1b smoke background. They are not physical branch values; M1c is required for that.

## Frozen Open DtN

`boundary_class=open_impedance`. The frozen outgoing branch is `Frozen target-blind D2 branch: Y_out,l=2(omega)=i Gamma_port omega^5 on physical tangential lanes at w=L and r=R_exit; Gamma_port=sigma_Q^can=4 a^5/(27 c_s^5).`.

Implementation: `K_ret=K_cons+i Gamma_port omega^5 B_tan`, where `B_tan` is assembled from the physical tangential electric-field traces `E_r/E_E/E_B` at `w=L` and `E_E/E_B/E_w` at `r=R_exit`. The mouth has no injected Maxwell drive; only a scalar gauge anchor is added. This is an open impedance/DtN loss term, not a hard cap.

## Gamma Port

The compact retarded quadrupole denominator gives `sigma_Q^can=9/(8 Omega_Q^5)` with `Omega_Q=3 c_s/(2a)`, hence `Gamma_port=4 a^5/(27 c_s^5)`. The older V2-09 `a^5/(27 c_s^5)` is the residue-reduced wall coefficient after the minimal isotropic pole weight `c1=1/4`; it is not the physical outgoing-port denominator normalization used for Spike-2 extraction.

- `Gamma_port`: `0.14814814814814817`.
- legacy V2-09 reduced coefficient: `0.037037037037037035`.
- factor gap recorded: `4.000000`.
- chi_Q reference check: `1.0`, computed after the derivation on `Sigma_ref=-i sigma_Q^can omega^5`; it was not used as a fit knob.

The Gamma-port value identities are construction checks, not runtime physics gates. `Gamma_port` being derived rather than fitted is established by code structure (`Gamma_port` depends only on `c_s` and `a`) together with the transliteration-fidelity audit; the runtime identity table only records the deterministic formula/value ledger.

## Forward Source

The chain is `eta -> delta V_conf -> delta psi -> delta J_psi`: the wall mode is the M1b weak-form wall eigen-shape, `delta V_conf=-4 V_radial r^4/R0^5 eta`, the M1b smoke BdG amplitude/phase response is solved on the transfer grid, and the Step-8c phasor Frechet current forms `j_eta`. The external source uses `delta A_i=0` and `A_i0=0`; the solved Maxwell field is the response, not an extra source. `S_eta^(A)` is not used.

Current Frechet check against a central finite difference of the nonlinear gauge-covariant current: max relative error `3.802185155263707e-10`.

## Coefficients

| grid | h | Z0 | Z2 | Z4 | N0 | N2 | N4 | max residual |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 7x7 | 0.250000 | 2.713075432486854e-6 | -1.4337939203737365e-6 | 4.7185230702415016e-7 | 5.396455366219856e-6 | -2.758126168323624e-6 | 9.388951783904308e-7 | 8.178243854926708e-16 |
| 9x9 | 0.200000 | 2.6902620249726532e-6 | -1.4003128613686532e-6 | 4.580104505776022e-7 | 5.117894519839075e-6 | -2.632742996398106e-6 | 9.009338162931466e-7 | 1.1911631422396414e-15 |
| 11x11 | 0.166667 | 2.6733784164555218e-6 | -1.3784738025851937e-6 | 4.4921775805827964e-7 | 4.947131039276808e-6 | -2.5567729980147607e-6 | 8.774217407344629e-7 | 1.579750979448191e-15 |

Primary reported values are the finest mesh row (`11x11`, omega window `0.040000`):

- `Z0,Z2,Z4 = 2.6733784164555218e-6, -1.3784738025851937e-6, 4.4921775805827964e-7`.
- `N0,N2,N4 = 4.947131039276808e-6, -2.5567729980147607e-6, 8.774217407344629e-7`.

Mesh convergence orders from successive-difference stabilization: `Z0 -> 1.348930628692854`, `N0 -> 2.1930146181374166`.

Z0 order diagnostic: A diagnostic-only 13x13 grid leaves the reported 7x7/9x9/11x11 rows unchanged and gives a 9x9/11x11/13x13 Z0 order of 1.92845075131403. The Z0 successive differences are 2.2813407514200722e-8, 1.6883608517131443e-8, 1.1878678826418488e-8, while the largest Green-solve residual is 2.224704811797484e-15 and the conditioned residual-scale reference is 2.1832207148604144e-14. These gaps rule out round-off or omega-fit floor control. The added point sharpens Z0 from the 7x7/9x9/11x11 order but does not make the original sub-2 estimate disappear into noise, so this is a resolved coarse/pre-asymptotic conservative spatial/boundary-discretization effect on the tiny M1b smoke Z0 rather than a persistent proof of sub-2 asymptotic order. Since Z0 is fitted from Re Sigma_cons before the imaginary open-DtN loss is added, the open-port DtN is not the cause; plausible contributors are the near-r_min 1/r^2 coefficients and boundary closures. Carry this forward for M1c because Z0 precision feeds D0/R_norm.

Omega-window refinement on the `9x9` grid:

| omega window | Z0 | N0 | max residual | positive flux |
| --- | --- | --- | --- | --- |
| 0.070000 | 2.6902620197586555e-6 | 5.117894509240859e-6 | 1.2030986596204143e-15 | 1.2772914794668074e-12 |
| 0.052500 | 2.690262024198248e-6 | 5.117894518265619e-6 | 1.0266810858605531e-15 | 3.035492578602038e-13 |
| 0.040000 | 2.6902620249726532e-6 | 5.117894519839075e-6 | 1.1911631422396414e-15 | 7.799631267440177e-14 |

Omega-window orders from successive-difference stabilization: `Z0 -> 6.069972791922893`, `N0 -> 6.071624145042677`.

## Validation

Genuine pass checks:

| check | status |
| --- | --- |
| current_frechet_matches_step8c | PASS |
| outgoing_flux_positive | PASS |
| open_not_hard_cap | PASS |
| pure_gauge_zero_physical_transfer | PASS |
| basis_invariance | PASS |
| v2_09_regression | PASS |
| green_residuals_small | PASS |
| bdg_residuals_small | PASS |
| N0_positive | PASS |

Construction identity checks:

| identity check | value | note |
| --- | --- | --- |
| gamma_port_derived_not_fitted | TRUE | Construction identity: Gamma_port is computed solely from c_s and a through sigma_Q^can; this records code structure plus the transliteration-fidelity audit, not a runtime physics failure mode. |
| gamma_port_uses_compact_sigma_Q_can | TRUE | Construction identity: deterministic equality Gamma_port=4 a^5/(27 c_s^5). |
| legacy_factor_four_recorded | TRUE | Construction identity: deterministic factor-four ledger against the V2-09 residue-reduced coefficient. |

- Basis invariance: max relative `Z` delta `3.149126899656948e-16`, max relative `N` delta `7.057516292941943e-16`.
- V2-09 regression: relative `Z` error `2.179611917552392e-16`, relative `N0` error `2.8119093166491344e-16`.
- Outgoing flux: `4.089965809183173e-14`; hard-cap/open boundary norm ratio `3.080820709977718e-7`.
- Pure gauge: physical field norm `2.2024049865573226e-15`, outgoing flux `-2.9619777649206203e-28`.

## Term Map

| term | implementation | source |
| --- | --- | --- |
| open exit/radial DtN | K_ret=K_cons+i Gamma_port omega^5 B_tan, with B_tan assembled from E_r/E_E/E_B traces at w=L and E_E/E_B/E_w traces at r=R_exit | Decision 05 D2; scratch log D2 lines 9617-9648 |
| outgoing flux | Phi_out=(Gamma_port omega^5/2mu0) A^dagger B_tan A >= 0, equivalent to D2 Poynting form under frozen DtN branch | scratch log D2 physical flux lines 9635-9644 |
| wall drive | delta V_conf=-4 V_radial r^4/R0^5 chi_eta(w) | M1b mt15_02_bdg_wall_derivation.wls driveCoeff/delta V ledger |
| BdG response | M1b smoke L_plus/L_phase response on the same VSH transfer grid; no self-consistent M1c claim | mt15_02_bdg_wall_derivation.wls finite-dimensional BdG sector |
| phasor current | Step-8c Frechet formula with delta A_i=0 for the external j_eta source and A_i0=0 on the M1b smoke background | step8c report lines 31-32; scratch log D3 lines 9652-9687 |
| S_eta^(A) | not used; forward source only | Decision 04 and scratch log D3 lines 9684-9687 |
| self-energy extraction | Sigma_A=<j_eta,A>, Z_n from Re Sigma_cons, N_n from -Im Sigma_ret/(Gamma_port omega^5) | Decision 05 D4; scratch log D4 lines 9690-9722 |
| Gamma_port | Gamma_port=sigma_Q^can=9/(8 Omega_Q^5)=4a^5/(27c_s^5), not fitted to chi_Q | compact lines 2559-2568; notes lines 41400-41415 |

Diagnostics JSON: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_04_spike2_transfer_n0/mt15_04_spike2_transfer_n0_diagnostics.json`.