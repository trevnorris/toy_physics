# M1c Physical Frozen Derived Chain Report

**LOUD LABEL:** PHYSICAL frozen decision-07 effective-closure branch. The freeze hash was written before the WP1 solve; the result is conditional on the frozen target-blind posits and modest CPU resolution.

## Torch Background

- Background JSON: `/var/projects/toy_physics/software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/wp1_background_10x8.json`
- Content hash: `6fab0c757e344c565f1903f4dc809ef2cab2a48443c1b56328dc31610f60e176`
- Grid: `10x8`
- Coupled stationary residual linf: `4.5892669664482355e-9`
- Gauge field linf `(A_00,A_r0,A_w0)`: `0.016195563892853235, 0., 0.`

## Import Fidelity

- Smooth-field interpolation order for the primary extraction: `1`.
- Conservative rho projection: finite-volume overlap in the imported `r^2 dr dw` measure.
- Source-point self-match only: the interpolant reproduces its own source grid nodes, so this is trivially exact by construction and is NOT an interpolation-convergence study.
- Source-point self-match `psi_R0` max abs delta: `5.551115123125783e-17`
- Source-point self-match `R0(w)` max abs delta: `0.`
- BdG projected rho mass relative error: `1.7439342475037903e-16`

## Interpolation-Method Sensitivity

Diagnostic-only re-extraction on the same torch background using order-2 imported smooth-field interpolation; the primary order-1 extraction above remains unchanged.

| coefficient | order1 primary | order2 diagnostic | relative delta | relative abs delta |
| --- | --- | --- | --- | --- |
| K | 4.060384354904323 | 4.060384354904323 | 0. | 0. |
| M | 0.5258858064705113 | 0.5258858064705113 | 0. | 0. |
| B0 | 0.004659423689433822 | 0.004659434120662928 | 2.238738050338406e-6 | 2.238738050338406e-6 |
| B2 | 0.000515232316768914 | 0.0005152342236586393 | 3.701028959695537e-6 | 3.701028959695537e-6 |
| B4 | 0.00005749912947089773 | 0.00005749942343049828 | 5.1124182793032915e-6 | 5.1124182793032915e-6 |
| Z0 | 2.5002239096227484e-6 | 2.51711874947392e-6 | 0.006757330727919082 | 0.006757330727919082 |
| Z2 | -1.7477512766478496e-6 | -1.7624919115242904e-6 | -0.008434057564937359 | 0.008434057564937359 |
| Z4 | 6.714367676645645e-7 | 6.778284195356489e-7 | 0.009519365305710407 | 0.009519365305710407 |
| N0 | 2.674541095707681e-6 | 2.6753511114777304e-6 | 0.0003028615904797882 | 0.0003028615904797882 |
| N2 | -1.5324096460422127e-6 | -1.528642022193304e-6 | 0.0024586270770609354 | 0.0024586270770609354 |
| N4 | 5.228908789315485e-7 | 5.196447384806385e-7 | -0.006208064783120649 | 0.006208064783120649 |

- Max relative abs delta: `0.009519365305710407`
- Carry-forward note: Order2-vs-order1 interpolation-method sensitivity is a diagnostic on the same frozen physical torch background. It is distinct from the WP1 background-resolution term and is carried explicitly into the M1c section J numerical-only budget.

## Gauge Activation

The imported background activates scalar `A_00` in the BdG `q A0` diagonal. The isotropic WP1 solve returns spatial gauges with linf values shown above; they are numerically zero at this resolution, so the spatial-gauge sign flip is a near-null sensitivity check rather than a large effect.

- BdG active-minus-zero-A0 deltas `(B0,B2,B4)`: `-5.8083120077889044e-6, -1.3377597358199544e-6, -2.2692803154130407e-7`
- Spatial gauge sign-flip relative deltas `(Z0,N0)`: `0., 0.`

## Coefficients

| coefficient | value |
| --- | --- |
| K | 4.060384354904323 |
| M | 0.5258858064705113 |
| B0 | 0.004659423689433822 |
| B2 | 0.000515232316768914 |
| B4 | 0.00005749912947089773 |
| Z0 | 2.5002239096227484e-6 |
| Z2 | -1.7477512766478496e-6 |
| Z4 | 6.714367676645645e-7 |
| N0 | 2.674541095707681e-6 |
| N2 | -1.5324096460422127e-6 |
| N4 | 5.228908789315485e-7 |

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

- BdG max relative eigen residual: `1.9046513612373904e-13`
- Stieltjes `B0*B4-B2^2`: `2.4484657354608302e-9`
- Current Frechet max relative error: `2.050267008458272e-10`
- Green residual max: `1.3508640790632167e-15`
- Pure gauge physical field norm / flux: `2.4468926817079183e-15 / -1.4500007483594576e-27`
- Basis invariance max relative Z/N: `6.845950329670515e-16 / 1.5726479248121067e-15`
- V2-09 relative Z/N errors: `2.179611917552392e-16 / 2.8119093166491344e-16`
- N0: `2.674541095707681e-6`

## Direct Formula Preview

These values are physical-run diagnostics computed after the freeze; they are not exported inside the V2-22B handoff packet because the validator forbids target residual fields there.

- `D0 = 4.055722430990979`
- `R_norm = -10.799999340551249` (PHYSICAL frozen diagnostic; V2 remains authoritative)
- `R_pole = -0.8310527171393044`
- `P2 = -2.0665687783479044e-7`
- `P4 = 6.419188997684018e-8`

## Artifacts

- Direct-derived packet: `/var/projects/toy_physics/software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_v2_22b_physical_frozen_packet.json`
- Diagnostics JSON: `/var/projects/toy_physics/software/stage1_solver/frozen/m1c/834835c999efd572e7aba46e04831a201ef51a3a7c88820459964b44b303d5e8/m1c_physical_derived_diagnostics.json`
- Packet content hash: `e25b699964ea27bbd0c294d354f9ea03307886d405dbaaedf72291fd1b58766d`
- Diagnostics content hash: `837da5268028070c5e25e6de3d235389a64d836bc52d1f12cd244c1b3adc4bd8`