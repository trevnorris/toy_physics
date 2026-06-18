# Path-A B2a BdG Derivation

Overall B2a gate: PASS
Validated bundle: `software/stage1_solver/runs/patha_b2a_bdg_derivation/bundles/patha_b2a_validated_bdg_bundle_tau_1.json`
Bundle content hash: `58ae88e8e6479a96435192f7f60a1c37616b888e62bbe0d9ebb6070fe5a098dd`

## Load-Bearing Sources

- Decision-12 requires Reading A, deriving spectra on the Path-A background, not inheriting the M1b packet: `software/stage1_solver/decisions/12_pathA_b2_derive_backreaction_bundle.md:18` and `:69`.
- Decision-11 freezes the Hooke family and geometry `a=1`, `L=37/20`: `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md:45` and `:74`.
- Decision-11 clarifies numeric `varpi`/mode data are derived outputs, not frozen numbers: `software/stage1_solver/decisions/11_pathA_gate_a_freeze_sheet.md:166`.
- Canonical EOS is `P=K rho^5`, `U=(K/4)rho^5`, `h=(5K/4)rho^4`: `notes/moving_throat_pde_program_compact.md:576`.
- Canonical BdG and wall source are the compact BdG skeleton and `delta V_conf`: `notes/moving_throat_pde_program_compact.md:1406` and `:1424`.

## Background

`tau=1` closed residual Linf `2.749412086889e-10` versus old smoke residual `2.433925092213e+02`.
`tau=2` closed residual Linf `2.747768956812e-10`.
Frozen geometry used by the exported background: `a=1.0`, `L=1.85`.

## Derived Tau-1 Bundle

| mode | varpi | coupling | overlap | lambda_B | eigen_residual |
| --- | --- | --- | --- | --- | --- |
| patha_b2a_bdg_mode_1 | 4.528666496834e+00 | 2.645837444635e-02 | 8.811078907386e-01 | 3.002852967776e-02 | 9.630620076447e-14 |
| patha_b2a_bdg_mode_2 | 5.959660024989e+00 | 1.136118164667e-02 | 3.848036469310e-01 | 2.952462051043e-02 | 6.880002450397e-14 |
| patha_b2a_bdg_mode_3 | 1.010910694416e+01 | 3.280848284762e-03 | 1.131408301979e-01 | 2.899791595149e-02 | 3.271470116656e-14 |
| patha_b2a_bdg_mode_4 | 1.375819948807e+01 | 1.330317549368e-02 | 8.811155400730e-01 | 1.509810562707e-02 | 2.759972027811e-14 |
| patha_b2a_bdg_mode_5 | 1.518919329243e+01 | 5.712318547954e-03 | 3.848060343616e-01 | 1.484466988006e-02 | 2.884357447664e-14 |
| patha_b2a_bdg_mode_6 | 1.657300997187e+01 | 8.404195176168e-04 | 3.015605634539e-02 | 2.786901271145e-02 | 1.805961692080e-14 |
| patha_b2a_bdg_mode_7 | 1.933863972443e+01 | 1.649597993859e-03 | 1.131421356235e-01 | 1.457987322555e-02 | 2.233718507644e-14 |
| patha_b2a_bdg_mode_8 | 2.471812192080e+01 | 6.897649967821e-04 | 2.618664141890e-02 | 2.634033840950e-02 | 1.228896741163e-14 |
| patha_b2a_bdg_mode_9 | 2.565341680133e+01 | 4.952086633955e-03 | 8.811023450162e-01 | 5.620330784460e-03 | 1.405918103858e-14 |
| patha_b2a_bdg_mode_10 | 2.580254285789e+01 | 4.230973678342e-04 | 3.019492874994e-02 | 1.401219957623e-02 | 1.563550480503e-14 |
| patha_b2a_bdg_mode_11 | 2.708441005691e+01 | 2.126375502062e-03 | 3.847950530769e-01 | 5.525994903155e-03 | 1.360756526877e-14 |
| patha_b2a_bdg_mode_12 | 3.123385682499e+01 | 6.140799340971e-04 | 1.131442400229e-01 | 5.427407828917e-03 | 1.337768756797e-14 |
| patha_b2a_bdg_mode_13 | 3.374707442850e+01 | 2.490953560781e-04 | 1.015868840492e-02 | 2.452042489632e-02 | 9.004425052534e-15 |
| patha_b2a_bdg_mode_14 | 3.394765489397e+01 | 3.467180812073e-04 | 2.618094680489e-02 | 1.324314524571e-02 | 9.998446522710e-15 |
| patha_b2a_bdg_mode_15 | 3.769776007820e+01 | 1.572915223741e-04 | 3.015492838911e-02 | 5.216113278215e-03 | 9.879568720845e-15 |
| patha_b2a_bdg_mode_16 | 3.914830204221e+01 | 2.497481586368e-03 | 8.811083173683e-01 | 2.834477370305e-03 | 8.660831774628e-15 |
| patha_b2a_bdg_mode_17 | 4.057929438014e+01 | 1.072394014492e-03 | 3.847968510752e-01 | 2.786909538100e-03 | 1.009291594827e-14 |
| patha_b2a_bdg_mode_18 | 4.277602696637e+01 | 2.446104743180e-04 | 1.084858420714e-02 | 2.254768637525e-02 | 6.489289882992e-15 |
| patha_b2a_bdg_mode_19 | 4.297660746677e+01 | 1.251022148570e-04 | 1.014777637403e-02 | 1.232804214894e-02 | 8.024450375197e-15 |
| patha_b2a_bdg_mode_20 | 4.472874182218e+01 | 3.096976369511e-04 | 1.131444520591e-01 | 2.737188004494e-03 | 8.168513938693e-15 |
| patha_b2a_bdg_mode_21 | 4.584287214802e+01 | 1.290933820085e-04 | 2.618530856900e-02 | 4.929992773174e-03 | 8.763735459458e-15 |
| patha_b2a_bdg_mode_22 | 5.092113898366e+01 | 9.647161500029e-05 | 4.678186679125e-03 | 2.062158302292e-02 | 6.548693862692e-15 |
| patha_b2a_bdg_mode_23 | 5.119264522220e+01 | 7.932886202519e-05 | 3.015584690468e-02 | 2.630629551740e-03 | 9.683053728584e-15 |
| patha_b2a_bdg_mode_24 | 5.200556002853e+01 | 1.229612437781e-04 | 1.084639431450e-02 | 1.133660092127e-02 | 6.351850997364e-15 |
| patha_b2a_bdg_mode_25 | 5.291973547543e+01 | 1.472116905058e-03 | 8.811045940952e-01 | 1.670762943383e-03 | 7.961383661238e-15 |
| patha_b2a_bdg_mode_26 | 5.435072629380e+01 | 6.321191294286e-04 | 3.847991672903e-01 | 1.642724785191e-03 | 8.160892183960e-15 |
| patha_b2a_bdg_mode_27 | 5.487182475907e+01 | 4.659639195446e-05 | 1.015321555173e-02 | 4.589323620389e-03 | 9.351441822726e-15 |
| patha_b2a_bdg_mode_28 | 5.738504223996e+01 | 8.103865247906e-05 | 4.266798604668e-03 | 1.899284685957e-02 | 7.762064982640e-15 |
| patha_b2a_bdg_mode_29 | 5.850017489780e+01 | 1.825489709561e-04 | 1.131443708612e-01 | 1.613416289000e-03 | 7.105829647468e-15 |
| patha_b2a_bdg_mode_30 | 5.933775734714e+01 | 6.510403299168e-05 | 2.618476709420e-02 | 2.486332330453e-03 | 7.271166579295e-15 |

`B0=3.901043467984e-05`, `B2=1.773353799229e-06`, `B4=8.407280302307e-08`.
Structural value `B0*B4-B2^2=1.349328934410e-13`.

## Modal Truncation

Exported mode count `K=30` was checked against `100` positive modes with tolerance `1.0e-04`.
Modal gate criterion: max_n |B_n(K)-B_n(all-positive)|/|B_n(all-positive)| <= 1.0e-04 at the export grid and finer confirmation grid.
| grid | M | B0 | B2 | B4 | rel_prev_max | rel_all_B0 | rel_all_B2 | rel_all_B4 | rel_all_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10x10 | 3 | 3.787336713451e-05 | 1.767701946899e-06 | 8.404386288118e-08 | - | 3.003713088641e-02 | 3.197350629269e-03 | 3.443459648998e-04 | 3.003713088641e-02 |
| 10x10 | 5 | 3.894974960870e-05 | 1.773254269586e-06 | 8.407261405648e-08 | 2.763515773542e-02 | 1.571892301029e-03 | 5.619061371694e-05 | 2.247934831068e-06 | 1.571892301029e-03 |
| 10x10 | 8 | 3.896037602437e-05 | 1.773284362432e-06 | 8.407270225271e-08 | 2.727493097523e-04 | 1.298714258736e-03 | 3.921954165535e-05 | 1.198885299198e-06 | 1.298714258736e-03 |
| 10x10 | 15 | 3.900463505030e-05 | 1.773350338861e-06 | 8.407280088976e-08 | 1.134711961322e-03 | 1.625286308102e-04 | 2.013685649016e-06 | 2.565017967833e-08 | 1.625286308102e-04 |
| 10x10 | 30 | 3.901043467984e-05 | 1.773353799229e-06 | 8.407280302307e-08 | 1.486686724530e-04 | 1.383579544145e-05 | 6.236892798227e-08 | 2.756329326113e-10 | 1.383579544145e-05 |
| 10x10 | all_positive | 3.901097442023e-05 | 1.773353909832e-06 | 8.407280304624e-08 | 1.383560401487e-05 | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 |

Finer-grid modal confirmation `12x12`:
| grid | M | B0 | B2 | B4 | rel_prev_max | rel_all_B0 | rel_all_B2 | rel_all_B4 | rel_all_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12x12 | 3 | 3.773888506555e-05 | 1.756221234423e-06 | 8.326146992917e-08 | - | 2.950389930233e-02 | 3.089392355864e-03 | 3.271019537599e-04 | 2.950389930233e-02 |
| 12x12 | 5 | 3.879385775832e-05 | 1.761554364695e-06 | 8.328853334693e-08 | 2.719432285750e-02 | 1.507237881645e-03 | 5.252530802431e-05 | 2.059968264583e-06 | 1.507237881645e-03 |
| 12x12 | 8 | 3.880450694239e-05 | 1.761583792522e-06 | 8.328861748030e-08 | 2.744316294369e-04 | 1.232392618461e-03 | 3.581910724594e-05 | 1.049823631564e-06 | 1.232392618461e-03 |
| 12x12 | 15 | 3.884695341353e-05 | 1.761644073668e-06 | 8.328870334441e-08 | 1.092658945248e-03 | 1.383870883945e-04 | 1.599194704369e-06 | 1.890108868203e-08 | 1.383870883945e-04 |
| 12x12 | 30 | 3.885201613232e-05 | 1.761646841785e-06 | 8.328870491046e-08 | 1.303077495538e-04 | 8.061305930666e-06 | 2.786920538100e-08 | 9.845566797074e-11 | 8.061305930666e-06 |
| 12x12 | all_positive | 3.885232933031e-05 | 1.761646890880e-06 | 8.328870491866e-08 | 8.061240946536e-06 | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 | 0.000000000000e+00 |

## Dual Engine

Scope: the Python and Mathematica paths independently assemble the BdG matrix, solve the eigensystem, and perform the overlap quadrature. The closed background and B1 wall mode `chi` are shared Python inputs by directive pathA_09; those common-mode inputs were validated earlier by the Path-A closed background checks in chunks 1b/1c and the B1 analytic `chi` oracle.

| tau | grid | varpi_abs | varpi_rel | coupling_abs | coupling_rel | moments_abs | moments_rel |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1.000000000000e+00 | 6x6 | 3.623767952377e-13 | 7.620011948101e-15 | 1.239676764020e-15 | 6.395173076911e-11 | 2.778268066994e-19 | 2.084337840890e-14 |
| 1.000000000000e+00 | 8x8 | 2.984279490192e-13 | 8.722258220544e-15 | 8.380943779685e-15 | 2.059027116712e-10 | 1.084202172486e-18 | 7.236922156299e-14 |
| 1.000000000000e+00 | 10x10 | 6.110667527537e-13 | 3.844022746716e-14 | 1.771044248755e-15 | 4.185902308541e-12 | 2.419126097358e-18 | 2.139361909275e-13 |
| 2.000000000000e+00 | 10x10 | 4.973799150321e-13 | 2.295125062329e-14 | 1.007006977805e-15 | 2.394823674225e-12 | 1.551764359370e-18 | 9.252067007830e-14 |
| 1.000000000000e+00 | 12x12 | 4.263256414561e-13 | 3.545022593392e-14 | 4.422460661568e-15 | 1.020868045325e-11 | 1.077425908907e-18 | 1.364983445372e-13 |

Max dual-engine diffs: varpi abs/rel `6.111e-13`/`3.844e-14`, c abs/rel `8.381e-15`/`2.059e-10`, B abs/rel `2.419e-18`/`2.139e-13`.
Dual-engine criterion: AND gate: varpi abs/rel <= 1e-9/1e-9, coupling abs/rel <= 1e-12/1e-8, B moments abs/rel <= 1e-16/1e-9.
Max eigen residual across exported packets: `1.776797789129e-13`.

## Grid Convergence

Criterion: finest-vs-previous max rel change: primary-three varpi <= 0.08 and B moments <= 0.20; full exported varpi max is recorded as an error bar.
| grid | varpi | B0 | B2 | B4 | rel_change_varpi_max | rel_change_varpi_primary3_max | rel_change_B0 | rel_change_B2 | rel_change_B4 | rel_change_B_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 6x6 | [4.491991328594e+00, 5.902131029806e+00, 9.751479882869e+00, 1.293123885097e+01, 1.434137698479e+01, 1.501061847995e+01, 1.819072691754e+01, 2.026975725586e+01, 2.195084063066e+01, 2.336097226513e+01, 2.344986601235e+01, 2.411910574245e+01, 2.721032740708e+01, 2.792878427402e+01, 2.870900526468e+01, 2.933888747143e+01, 3.246946776244e+01, 3.255835483647e+01, 3.318826585835e+01, 3.544715250320e+01, 3.685732229621e+01, 3.772860831474e+01, 3.844741146160e+01, 4.070664607955e+01, 4.157796310995e+01, 4.370655722475e+01, 4.596577965549e+01, 4.755593530637e+01, 5.122491342382e+01, 5.507423679392e+01] | 4.030919759559e-05 | 1.858051573127e-06 | 8.953056578981e-08 | - | - | - | - | - | - |
| 8x8 | [4.517219980893e+00, 5.941589197903e+00, 9.994467727500e+00, 1.349919068701e+01, 1.492355971957e+01, 1.606097089230e+01, 1.897643814210e+01, 2.321700121359e+01, 2.449001960775e+01, 2.504294148902e+01, 2.591438735095e+01, 2.996726670490e+01, 3.037303156944e+01, 3.219897192177e+01, 3.586106300015e+01, 3.603377029251e+01, 3.643953495697e+01, 3.728542808231e+01, 3.935500238516e+01, 4.049241313072e+01, 4.133830949065e+01, 4.318980084206e+01, 4.542150595408e+01, 4.563016940650e+01, 4.705452289346e+01, 4.740481349655e+01, 4.947438401746e+01, 5.034583142481e+01, 5.039698563768e+01, 5.110741325115e+01] | 3.941943525890e-05 | 1.799994574235e-06 | 8.577076252223e-08 | 1.994379735882e-01 | 2.431223465379e-02 | 2.257166625674e-02 | 3.225398549718e-02 | 4.383548842310e-02 | 4.383548842310e-02 |
| 10x10 | [4.528666496834e+00, 5.959660024989e+00, 1.010910694416e+01, 1.375819948807e+01, 1.518919329243e+01, 1.657300997187e+01, 1.933863972443e+01, 2.471812192080e+01, 2.565341680133e+01, 2.580254285789e+01, 2.708441005691e+01, 3.123385682499e+01, 3.374707442850e+01, 3.394765489397e+01, 3.769776007820e+01, 3.914830204221e+01, 4.057929438014e+01, 4.277602696637e+01, 4.297660746677e+01, 4.472874182218e+01, 4.584287214802e+01, 5.092113898366e+01, 5.119264522220e+01, 5.200556002853e+01, 5.291973547543e+01, 5.435072629380e+01, 5.487182475907e+01, 5.738504223996e+01, 5.850017489780e+01, 5.933775734714e+01] | 3.901043467984e-05 | 1.773353799229e-06 | 8.407280302307e-08 | 1.518296388475e-01 | 1.134019229375e-02 | 1.048438917479e-02 | 1.502282004687e-02 | 2.019629937517e-02 | 2.019629937517e-02 |
No Richardson order is claimed on this non-doubling 6->8->10 ladder.

## Tau Sensitivity

`tau=1` vs `tau=2` max relative movement `1.020122727435e-02` (max relative movement in {varpi,B0,B2,B4} >= 1e-5).
Varpi max abs/rel `1.201701373716e-04`/`2.543471465665e-05`; B max abs/rel `3.939356828084e-07`/`1.020122727435e-02`.
Background movement: rho0 max abs/rel `1.306729603811e-06`/`2.400452138295e-04`, R0(w) max abs/rel `1.166200769463e-03`/`1.167548610810e-03`, A00 max abs/rel `2.229677904537e-08`/`6.509365079410e-05`.
B2c design input: doubling tau moves the matter background at sub-percent scale in this frozen Hooke family; tau leverage on `R_norm` is expected to be dominated by exact `K=tau*kappahat` and the wall/Maxwell sectors, not by the BdG `B_n` drift.

## B1 Consumer Check

`patha_extraction.bdg_moments` applied to MMA-engine `bdg_modes[]` versus Python-engine converged `B_n` max abs/rel `2.433e-18`/`2.144e-13`.

## Error Budget

Recorded B-moment spatial ladder rel errors: B0 `1.048e-02`, B2 `1.502e-02`, B4 `2.020e-02`; max `2.020e-02`.
Recorded B-moment finer-grid confirmation rel errors: B0 `4.061e-03`, B2 `6.602e-03`, B4 `9.326e-03`; max `9.326e-03`.
Recorded B-moment modal truncation rel errors at K: B0 `1.384e-05`, B2 `6.237e-08`, B4 `2.756e-10`; max `1.384e-05`.
Recorded varpi spatial ladder max rel: exported `1.518e-01`, primary-three `1.134e-02`; finer-grid confirmation exported `1.268e-01`, primary-three `6.217e-03`.

## Gates And Wrong Answers Caught

| gate | status | catches |
| --- | --- | --- |
| self_consistency | PASS | old smoke profile, stale background bundle, or M1a L=2/R_exit geometry |
| dual_engine_agreement | PASS | engine-specific matrix assembly, interpolation, normalization, overlap bug, or tiny-B abs-only false pass |
| eigen_residual | PASS | wrong eigenvector/eigenvalue pair exported from either engine |
| grid_convergence | PASS | under-resolved BdG spectrum or moments being shipped as converged |
| modal_truncation | PASS | under-truncated B_n moment sum with a fat high-mode tail |
| structural_sanity | PASS | sign/index/exponent bug in B0/B2/B4 moment assembly |
| tau_sensitivity | PASS | stale background reuse instead of re-solving at the requested tau |
| b1_consumer_moments | PASS | MMA packet whose B1-consumed modes no longer reproduce the independently assembled Python converged B_n |

## Adaptations

- The old inline smoke `rho0` and M1a `L=2.0`, `R_exit=1.65` were replaced by the exported closed Path-A background and frozen `a=1`, `L=1.85`.
- The BdG operator form was kept to the canonical matter-sector terms: l=2 kinetic, confinement, quintic enthalpy plus `rho dh/drho`, `-mu`, `q A0`, and the wall-drive kernel.
- The v2_22a overlap algebra is used with B1's flat-normalized wall mode. The Path-A bundle records both `lambda_B`, `overlap_I_eta_phi`, and `coupling` so B1 can either adapt profiles or call `bdg_moments` directly.
- The radial volume convention follows the Path-A finite-volume grid (`4*pi*r^2 dr`) and the closed wall-return source, rather than the old smoke script's ad hoc radial volume normalization.
- No Maxwell transfer, `R_norm`, `R_pole`, `P2`, `P4`, or root-find is assembled in this chunk.

## Files Created Or Modified

- `software/stage1_solver/src/stage1_solver/patha_b2a_bdg.py`
- `software/stage1_solver/mathematica/mt15_02_bdg_wall_derivation.wls`
- `software/stage1_solver/mathematica/mt15_02_patha_b2a_bdg_wall_derivation.wls`
- `software/stage1_solver/tests/test_patha_b2a_bdg.py`
- `software/stage1_solver/reports/patha_b2a_bdg_derivation_report.md`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/backgrounds/patha_b2a_closed_background_tau_1_61184d4ac66cebaf163fecd8bcbb4b691145f4cf97c5f85c82f0e1b2209e11c3.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/backgrounds/patha_b2a_closed_background_tau_1_latest.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/backgrounds/patha_b2a_closed_background_tau_2_356975cfc6afbcd8d563f1063140b0783d08798ba5fac5860330641d75d46468.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/backgrounds/patha_b2a_closed_background_tau_2_latest.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/bundles/patha_b2a_validated_bdg_bundle_tau_1.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/mathematica/patha_b2a_mma_tau_1_nr_10_nw_10.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/mathematica/patha_b2a_mma_tau_1_nr_12_nw_12.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/mathematica/patha_b2a_mma_tau_1_nr_6_nw_6.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/mathematica/patha_b2a_mma_tau_1_nr_8_nw_8.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/mathematica/patha_b2a_mma_tau_2_nr_10_nw_10.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/python/patha_b2a_python_tau_1_nr_10_nw_10.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/python/patha_b2a_python_tau_1_nr_12_nw_12.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/python/patha_b2a_python_tau_1_nr_6_nw_6.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/python/patha_b2a_python_tau_1_nr_8_nw_8.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/python/patha_b2a_python_tau_2_nr_10_nw_10.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/wall_inputs/patha_b2a_wall_input_tau_1.json`
- `software/stage1_solver/runs/patha_b2a_bdg_derivation/wall_inputs/patha_b2a_wall_input_tau_2.json`
