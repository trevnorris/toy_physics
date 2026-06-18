# MT15-06 M1c-Prep Cross-Engine Background Report

**LOUD LABEL:** real self-consistent torch WP1 background, but PLACEHOLDER unfrozen params. This is MACHINERY, NOT physical; it proves the cross-engine pipeline only. No Gate A, no frozen committed packet, no physical claim.

## Torch Background

- Background JSON: `/var/projects/toy_physics/software/stage1_solver/runs/m1c_prep_background_export/background_bundles/m1c_prep_wp1_background_latest.json`
- Content hash: `c2a5e04bf95a9d06`
- Grid: `10x8`
- Coupled stationary residual linf: `9.840448234044175e-9`
- Gauge field linf `(A_00,A_r0,A_w0)`: `0.011264993140235586, 0., 0.`

## Import Fidelity

- Smooth-field interpolation order for the primary extraction: `1`.
- Conservative rho projection: finite-volume overlap in the imported `r^2 dr dw` measure.
- Source-point self-match only: the interpolant reproduces its own source grid nodes, so this is trivially exact by construction and is NOT an interpolation-convergence study.
- Source-point self-match `psi_R0` max abs delta: `0.`
- Source-point self-match `R0(w)` max abs delta: `0.`
- BdG projected rho mass relative error: `1.7439342462105672e-16`

## Interpolation-Method Sensitivity

Diagnostic-only re-extraction on the same torch background using order-2 imported smooth-field interpolation; the primary order-1 extraction above remains unchanged.

| coefficient | order1 primary | order2 diagnostic | relative delta | relative abs delta |
| --- | --- | --- | --- | --- |
| K | 3.8382460113833874 | 3.841222154038521 | 0.0007753913236168853 | 0.0007753913236168853 |
| M | 1.077913573060997 | 1.0787493778931572 | 0.000775391323616697 | 0.000775391323616697 |
| B0 | 0.041439858685538886 | 0.04156471976322465 | 0.003013067168815824 | 0.003013067168815824 |
| B2 | 0.002433488571316712 | 0.0024382584630569157 | 0.001960104434607176 | 0.001960104434607176 |
| B4 | 0.00015666195504379751 | 0.0001567693864894058 | 0.0006857532550148559 | 0.0006857532550148559 |
| Z0 | 0.000018411351589721476 | 0.000018379842539751127 | -0.0017113925513181443 | 0.0017113925513181443 |
| Z2 | -0.000010347056376519815 | -0.000010326957720268063 | 0.0019424516036619943 | 0.0019424516036619943 |
| Z4 | 3.1487586245414714e-6 | 3.128499414911011e-6 | -0.006434030691511185 | 0.006434030691511185 |
| N0 | 0.000015772025302770683 | 0.00001552611992443404 | -0.015591236611409905 | 0.015591236611409905 |
| N2 | -6.613588503095876e-6 | -6.40454220067724e-6 | 0.03160860436369441 | 0.03160860436369441 |
| N4 | 1.8408150408808162e-6 | 1.7376619696071449e-6 | -0.056036629961646406 | 0.056036629961646406 |

- Max relative abs delta: `0.056036629961646406`
- Carry-forward note: Order2-vs-order1 interpolation-method sensitivity is a diagnostic on the same coarse 10x8 torch background. It is distinct from, and smaller than, the unresolved background-resolution error of that coarse background. The production background-resolution term and full interpolation error belong in the M1c-run section J error budget; this M1c-prep machinery result is not resolution-converged.

## Gauge Activation

The imported background activates scalar `A_00` in the BdG `q A0` diagonal. The isotropic WP1 solve returns spatial gauges with linf values shown above; they are numerically zero at this resolution, so the spatial-gauge sign flip is a near-null sensitivity check rather than a large effect.

- BdG active-minus-zero-A0 deltas `(B0,B2,B4)`: `-0.000025832958400237582, -3.1928592257959333e-6, -3.094033149404047e-7`
- Spatial gauge sign-flip relative deltas `(Z0,N0)`: `0., 0.`

## Coefficients

| coefficient | value |
| --- | --- |
| K | 3.8382460113833874 |
| M | 1.077913573060997 |
| B0 | 0.041439858685538886 |
| B2 | 0.002433488571316712 |
| B4 | 0.00015666195504379751 |
| Z0 | 0.000018411351589721476 |
| Z2 | -0.000010347056376519815 |
| Z4 | 3.1487586245414714e-6 |
| N0 | 0.000015772025302770683 |
| N2 | -6.613588503095876e-6 |
| N4 | 1.8408150408808162e-6 |

## A-not-zero Revalidation

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

- BdG max relative eigen residual: `1.430724378737486e-13`
- Stieltjes `B0*B4-B2^2`: `5.701826516861627e-7`
- Current Frechet max relative error: `8.479926814494943e-11`
- Green residual max: `1.1003929175399927e-15`
- Pure gauge physical field norm / flux: `2.731377278937409e-15 / -2.4904021667792283e-27`
- Basis invariance max relative Z/N: `7.366015586981875e-16 / 2.0112880819057522e-15`
- V2-09 relative Z/N errors: `2.179611917552392e-16 / 2.8119093166491344e-16`
- N0: `0.000015772025302770683`

## Direct Formula Preview

These preview values are diagnostics only and are not exported inside the V2-22B packet.

- `D0 = 3.796787741346259`
- `R_norm = -0.0455630514263164` (MACHINERY, NOT physical)
- `R_pole = -3.5007754832235976`
- `P2 = 6.220903751447901e-7`
- `P4 = 5.028801634082825e-7`

## Artifacts

- Direct-derived packet: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_06_m1c_prep_crossengine/mt15_06_v2_22b_m1c_prep_direct_derived_packet.json`
- Diagnostics JSON: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/mt15_06_m1c_prep_crossengine/mt15_06_m1c_prep_crossengine_diagnostics.json`
- Packet content hash: `a4151954a3b5ce74ddca7ddbb87a178853694a61108c18e44db0021bc106e63f`
- Diagnostics content hash: `fa0729c99b12ce864c00f3107d3cd4d8c1d78515d9770a05a2bf3a10d11f7173`