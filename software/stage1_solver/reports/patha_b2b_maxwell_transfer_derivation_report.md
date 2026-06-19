# Path-A B2b Maxwell Transfer Derivation

Overall B2b gate: PASS
Validated bundle: `software/stage1_solver/runs/patha_b2b_maxwell_transfer/bundles/patha_b2b_validated_maxwell_transfer_tau_1.json`
Bundle content hash: `138b17080cd1a973f7a9c1fcb6a6dc209808f948140b5b781b5d5e255ad4157b`

## Scope

This run derives only `{Z0,Z2,Z4,N0,N2,N4}` on the Path-A closed background and emits the B1 `direct_coefficients` shape. No `D0`, `P0`, `R_norm`, `R_pole`, `P2`, `P4`, root-find, or `mt15_05` preview was run.

## Canonical Sources

- Five-lane VSH Maxwell operator and open DtN design: `software/stage1_solver/decisions/05_mixed_maxwell_spike_design.md` D1-D4.
- Localized Maxwell PDE and mixed invariants: `notes/moving_throat_pde_program_compact.md:592`, `:674`, `:769`.
- Compact outgoing normalization: `Gamma_port=4a^5/(27 c_s^5)` from `notes/moving_throat_pde_program_compact.md:2559`.
- The old U/W formulas were used only as a regression reference, not as the live extraction path: `research/pde_audit/notes/stage_v2_09_maxwell_mixed_kernel_derivation.md:117` and `:299`.

## Background

`tau=1` closed residual Linf `2.749412086889e-10` versus old smoke residual `2.433925092213e+02`.
`tau=2` closed residual Linf `2.747768956812e-10`.
Frozen geometry: `a=1.0`, `L=1.85`. Background hash: `61184d4ac66cebaf163fecd8bcbb4b691145f4cf97c5f85c82f0e1b2209e11c3`.
Rest-background gauge-field linf: `A_00=7.724847250734e-04`, `A_r0=0.000000000000e+00`, `A_w0=0.000000000000e+00`. The B2a rest defect has `A_r0 == A_w0 == 0`, as expected for no rest current; this calibration transfer therefore runs through scalar `A_00`, kinematic `j_E`, and `q*delta_rho` channels. Vector-current channels are held-out non-rest/excited-defect coverage.

## Derived Tau-1 Transfer

`Z0=2.395089299987e-08`, `Z2=-1.527859528310e-08`, `Z4=5.016962667932e-09`.
`N0=2.157666812539e-08`, `N2=-5.934003963233e-09`, `N4=3.697738812793e-09`.
`Gamma_port=1.481481481481e-01` using `4a^5/(27c_s^5)` at `a=1`, `c_s=1`.
Relative movement from the old unconverged 7x7/window-0.036/radial-scale-1.0 packet:

| Z0 | Z2 | Z4 | N0 | N2 | N4 |
| --- | --- | --- | --- | --- | --- |
| 1.003678606736e-02 | 5.777991280520e-03 | 4.243387178788e-01 | 4.443773657806e-01 | 1.549588843476e+00 | 5.304668907569e-01 |

Direct coefficients emitted for B1 shape:

| K | M | B0 | B2 | B4 | Z0 | Z2 | Z4 | N0 | N2 | N4 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 7.720944358201e+00 | 1.000000000000e+00 | 3.901043467984e-05 | 1.773353799229e-06 | 8.407280302307e-08 | 2.395089299987e-08 | -1.527859528310e-08 | 5.016962667932e-09 | 2.157666812539e-08 | -5.934003963233e-09 | 3.697738812793e-09 |

## Dual Engine

Correctness dual-engine comparison is primary Python second-order interior-grid FD versus a genuinely different Python staggered cell-centered high-order discretization. The second engine uses different point locations, cell quadrature, and fourth-order interior derivative rows with a second-order boundary closure; it does not read or reuse primary packet outputs. The old Python-Wolfram files remain useful only as a transcription check of the primary discretization.

| tau | grid | window | radial_scale | abs | rel |
| --- | --- | --- | --- | --- | --- |
| 1.000000000000e+00 | 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 9.603190625643e-10 | 4.368649426740e-02 |
| 2.000000000000e+00 | 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 9.493375510608e-10 | 4.367890574413e-02 |
Max dual-engine coefficient abs/rel diff: `9.603190625643e-10`/`4.368649426740e-02`. Criterion: AND gate: max coefficient abs <= 2.0e-09 and rel <= 5.0e-02.
Per-coefficient dual diffs:

| coefficient | max_abs | max_rel | passed |
| --- | --- | --- | --- |
| Z0 | 9.603190625643e-10 | 3.854967233861e-02 | 1 |
| Z2 | 6.979596766624e-10 | 4.368649426740e-02 | 1 |
| Z4 | 1.867000190383e-10 | 3.589145147707e-02 | 1 |
| N0 | 1.396391535945e-10 | 6.605575385946e-03 | 1 |
| N2 | 4.227451183920e-11 | 7.073718538306e-03 | 1 |
| N4 | 9.751019204886e-11 | 2.719969560053e-02 | 1 |

Independent staggered-high-order engine grid convergence:

| grid | window | radial_scale | Z0 | Z2 | Z4 | N0 | N2 | N4 | rel_change_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 31x13 | 2.800000000000e-02 | 5.000000000000e+00 | 2.539017365771e-08 | -1.631432725964e-08 | 5.268791178793e-09 | 2.168313087372e-08 | -6.060208411119e-09 | 3.634281583499e-09 | - |
| 43x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.494942011640e-08 | -1.597829459966e-08 | 5.186080214658e-09 | 2.147026909954e-08 | -5.977632889889e-09 | 3.594038697061e-09 | 2.103057105925e-02 |
| 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.491121206243e-08 | -1.597655495976e-08 | 5.203662686970e-09 | 2.143704181935e-08 | -5.976278475073e-09 | 3.600228620744e-09 | 3.378864728642e-03 |
Independent omega-window convergence:

| grid | window | radial_scale | Z0 | Z2 | Z4 | N0 | N2 | N4 | rel_change_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 47x17 | 5.800000000000e-02 | 5.000000000000e+00 | 2.491121185972e-08 | -1.597627214460e-08 | 5.108634022695e-09 | 2.143704172268e-08 | -5.976143407572e-09 | 3.554748347174e-09 | - |
| 47x17 | 4.600000000000e-02 | 5.000000000000e+00 | 2.491121201313e-08 | -1.597645137716e-08 | 5.154024829043e-09 | 2.143704179566e-08 | -5.976228663037e-09 | 3.576341203075e-09 | 8.806866061791e-03 |
| 47x17 | 3.600000000000e-02 | 5.000000000000e+00 | 2.491121205304e-08 | -1.597652613898e-08 | 5.184433837260e-09 | 2.143704181482e-08 | -5.976264556789e-09 | 3.590941989845e-09 | 5.865444361307e-03 |
| 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.491121206243e-08 | -1.597655495976e-08 | 5.203662686970e-09 | 2.143704181935e-08 | -5.976278475073e-09 | 3.600228620744e-09 | 3.695252914495e-03 |

Independent outward-radial convergence:

| grid | window | radial_scale | Z0 | Z2 | Z4 | N0 | N2 | N4 | rel_change_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 27x17 | 2.800000000000e-02 | 3.000000000000e+00 | 2.493353072004e-08 | -1.598304853814e-08 | 4.898450815296e-09 | 2.146266038813e-08 | -6.081467070698e-09 | 3.820714200099e-09 | - |
| 37x17 | 2.800000000000e-02 | 4.000000000000e+00 | 2.491566847819e-08 | -1.597802153592e-08 | 5.094815077625e-09 | 2.144174794873e-08 | -5.998970129142e-09 | 3.660686564853e-09 | 4.371519724789e-02 |
| 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.491121206243e-08 | -1.597655495976e-08 | 5.203662686970e-09 | 2.143704181935e-08 | -5.976278475073e-09 | 3.600228620744e-09 | 2.091749905650e-02 |
Independent full-axis max relative error bar `2.091749905650e-02` under criterion: >=3 levels, every coefficient's final relative increment <= 5.00%, and every coefficient's increments shrink by factor <= 0.98; no max-over-coefficients aggregation is used for pass/fail.

## Convergence, Window, Truncation

Chosen resolution: `{'grid': '47x17', 'window': 0.028, 'radial_scale': 5.0}`.
Criterion: >=3 levels, every coefficient's final relative increment <= 5.00%, and every coefficient's increments shrink by factor <= 0.98; no max-over-coefficients aggregation is used for pass/fail.

Mesh sweep:

| grid | window | radial_scale | Z0 | Z2 | Z4 | N0 | N2 | N4 | rel_change_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 31x13 | 2.800000000000e-02 | 5.000000000000e+00 | 2.456743476261e-08 | -1.576066831335e-08 | 5.073894610628e-09 | 2.232332002863e-08 | -6.150464585554e-09 | 3.777148345395e-09 | - |
| 43x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.404966503987e-08 | -1.536185435626e-08 | 5.012294045536e-09 | 2.166483917256e-08 | -5.967027463344e-09 | 3.697425367548e-09 | 3.074179285029e-02 |
| 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.395089299987e-08 | -1.527859528310e-08 | 5.016962667932e-09 | 2.157666812539e-08 | -5.934003963233e-09 | 3.697738812793e-09 | 5.565129432811e-03 |

Omega-window sweep:

| grid | window | radial_scale | Z0 | Z2 | Z4 | N0 | N2 | N4 | rel_change_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 47x17 | 5.800000000000e-02 | 5.000000000000e+00 | 2.395089280231e-08 | -1.527831965343e-08 | 4.924349414195e-09 | 2.157666801878e-08 | -5.933855024430e-09 | 3.647594654946e-09 | - |
| 47x17 | 4.600000000000e-02 | 5.000000000000e+00 | 2.395089295182e-08 | -1.527849433589e-08 | 4.968587888899e-09 | 2.157666809927e-08 | -5.933949060227e-09 | 3.671411194409e-09 | 8.903631312154e-03 |
| 47x17 | 3.600000000000e-02 | 5.000000000000e+00 | 2.395089299072e-08 | -1.527856719622e-08 | 4.998223469033e-09 | 2.157666812039e-08 | -5.933988626546e-09 | 3.687505833501e-09 | 5.929222716223e-03 |
| 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.395089299987e-08 | -1.527859528310e-08 | 5.016962667932e-09 | 2.157666812539e-08 | -5.934003963233e-09 | 3.697738812793e-09 | 3.735168096653e-03 |

Radial truncation sweep:

| grid | window | radial_scale | Z0 | Z2 | Z4 | N0 | N2 | N4 | rel_change_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 27x17 | 2.800000000000e-02 | 3.000000000000e+00 | 2.397397797193e-08 | -1.530312581864e-08 | 4.707900010139e-09 | 2.160576377333e-08 | -6.066328512081e-09 | 3.918409285616e-09 | - |
| 37x17 | 2.800000000000e-02 | 4.000000000000e+00 | 2.395892735685e-08 | -1.528609944665e-08 | 4.904547150661e-09 | 2.158505435383e-08 | -5.964028329213e-09 | 3.760336887751e-09 | 4.203676494502e-02 |
| 47x17 | 2.800000000000e-02 | 5.000000000000e+00 | 2.395089299987e-08 | -1.527859528310e-08 | 5.016962667932e-09 | 2.157666812539e-08 | -5.934003963233e-09 | 3.697738812793e-09 | 2.240708666023e-02 |

Primary per-coefficient outward-radial increments:

| from | to | Z0 | Z2 | Z4 | N0 | N2 | N4 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000000000e+00 | 4.000000000000e+00 | 6.281840110083e-04 | 1.113846736709e-03 | 4.009486186612e-02 | 9.594332802764e-04 | 1.715286669033e-02 | 4.203676494502e-02 |
| 4.000000000000e+00 | 5.000000000000e+00 | 3.354512494268e-04 | 4.911553330155e-04 | 2.240708666023e-02 | 3.886711511755e-04 | 5.059714514138e-03 | 1.692874432925e-02 |

Independent per-coefficient outward-radial increments:

| from | to | Z0 | Z2 | Z4 | N0 | N2 | N4 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 3.000000000000e+00 | 4.000000000000e+00 | 7.169079918670e-04 | 3.146198172309e-04 | 3.854198029522e-02 | 9.753141138276e-04 | 1.375185069775e-02 | 4.371519724789e-02 |
| 4.000000000000e+00 | 5.000000000000e+00 | 1.788919681678e-04 | 9.179551949844e-05 | 2.091749905650e-02 | 2.195325931056e-04 | 3.796953934525e-03 | 1.679280692367e-02 |

Recorded per-coefficient relative error bars:

| Z0 | Z2 | Z4 | N0 | N2 | N4 |
| --- | --- | --- | --- | --- | --- |
| 4.123939763204e-03 | 5.449393194564e-03 | 2.240708666023e-02 | 4.086406977178e-03 | 5.565129432811e-03 | 1.692874432925e-02 |

Primary max relative error bar `2.240708666023e-02`; independent full-axis max relative error bar `2.091749905650e-02`. Cross-engine converged-value max relative difference is `4.368649426740e-02`; interpret that against these convergence bars, not as precision beyond the demonstrated discretization spread.

## Tau Sensitivity

`tau=1` vs `tau=2` primary transfer max abs/rel movement `2.974392410314e-10`/`1.397791445582e-02` (tau=1 vs tau=2 max relative movement in Z/N coefficients >= 1e-6).
`tau=1` vs `tau=2` independent transfer max abs/rel movement `2.974520885938e-10`/`1.407085437240e-02`.
B2a input hashes moved from `58ae88e8e6479a96435192f7f60a1c37616b888e62bbe0d9ebb6070fe5a098dd` to `257197d56fb27785381b11ef8dc2c70e175a1d59b68dd050ce31328a429782da`; background hashes moved from `61184d4ac66cebaf163fecd8bcbb4b691145f4cf97c5f85c82f0e1b2209e11c3` to `356975cfc6afbcd8d563f1063140b0783d08798ba5fac5860330641d75d46468`.

## Basis Invariance And Physicality

Basis rotation check: max relative `Z` delta `1.370409911516e-16`, max relative `N` delta `3.167017596832e-16`.
Port-leak negative control: fixed-port relative movement `9.879334500494e-02`; failed as expected: `True`. Toggling to the fixed-port branch changes the tested base value by relative `9.223054445192e-01` and gives tested extraction movement `9.879334500494e-02`.
Outgoing physicality: `N0_positive=True`, `Sigma_cons_real=True`, `minus_Im_Sigma_ret_positive=True`, min `-Im Sigma_ret=2.768572431447e-18`, min flux `1.384286215723e-18`.

## N-Channel Robustness

Minimum `-Im Sigma_ret=2.768572431447e-18`; solve/imaginary floor `1.905662928519e-23`; signal/floor `1.452813291382e+05`.
Scaled N fit condition `3.303145292285e+01` with fit RMS `1.085277843974e-19` and `N0/fit_RMS=1.988123893359e+11`.
Individual N-term fit-floor checks over the fitted omega window:

| coefficient | window_contribution | to_fit_rms | above_floor |
| --- | --- | --- | --- |
| N0 | 2.157666812539e-08 | 1.988123893359e+11 | 1 |
| N2 | 6.152612669239e-12 | 5.669159011583e+07 | 1 |
| N4 | 3.975206726389e-15 | 3.662847028952e+04 | 1 |

## B1 Consumer Compatibility

Required direct keys present: `True`; finite: `True`; `N0>0`: `True`. Cross-engine Z/N abs/rel through the coefficient dict: `9.603190625643e-10`/`4.368649426740e-02`.
Validated the direct_coefficients input contract only. The actual lane_extract call is not made here because that function immediately assembles D0/P0/R_norm-scope derived quantities.

## Gates And Wrong Answers Caught

| gate | status | catches |
| --- | --- | --- |
| self_consistency | PASS | running the old smoke R0/Z(w) background, stale tau bundle, or wrong M1a geometry |
| dual_engine_agreement | PASS | a discretization-specific VSH assembly, Green solve, source normalization, or abs-only tiny-number false pass |
| basis_invariance | PASS | basis-dependent U/W-port leakage or rotating the operator without the source/boundary quadratic |
| outgoing_physicality | PASS | wrong DtN sign, bound non-radiating solution, anti-radiating solution, or N0<=0 |
| N_channel_robustness | PASS | shipping an N-channel below the solve/fit noise floor or from an ill-conditioned omega fit |
| convergence_truncation | PASS | under-resolved mesh, unstable omega-window fit, or radial truncation dependence |
| tau_sensitivity | PASS | reusing a frozen tau=1 background/BdG packet instead of re-solving at tau=2 |
| b1_consumer_compatibility | PASS | missing direct-coefficient keys, non-finite scalar, N0<=0, or cross-engine mismatch through B1 input shape |
| no_cant_fail_gates | PASS | report-only validations without a concrete negative-control failure mode |

Negative controls proving gates can fail:

| gate | failed_as_expected | wrong_answer |
| --- | --- | --- |
| self_consistency | 1 | a stale smoke-level background or wrong-geometry fallback |
| dual_engine_agreement | 1 | a second discretization that lands on a materially different coefficient |
| basis_invariance | 1 | a posited U/W-style fixed-port extraction that is not rebased with the physical lanes |
| outgoing_physicality | 1 | an anti-radiating or non-positive-N0 transfer |
| convergence | 1 | a grow-one/shrink-another sequence hidden by a max-over-coefficients convergence check |
| tau_sensitivity | 1 | a frozen tau=1 transfer relabelled as tau=2 |
| b1_consumer_compatibility | 1 | a B1 direct-coefficient lane with a cross-engine mismatch |

## Dual-Engine Feasibility Note

A genuine independent l=2 transfer discretization was built: the primary engine uses the original interior-node second-order FD grid; the independent engine uses staggered cell centers, cell-volume quadrature, and fourth-order interior derivative stencils with a separate boundary closure. Both read the same Path-A background and B2a BdG JSON inputs and use the same canonical self-energy definitions. The old Mathematica worker is now labelled transcription-only because it uses the same stencil/grid/quadrature as the primary path.

## Adaptations

- The old hardcoded smoke `R0(w)` and `Z(w)` were replaced by the exported Path-A `A0={A_00,A_r0,A_w0}`, `Z_w`, and frozen geometry.
- The VSH operator and DtN/self-energy forms are retained from Spike-1/Spike-2; only the background interpolation and forward source changed.
- The live transfer uses `Sigma=<j,Gj>` and never constructs U/W mixed ports. The U/W formulas remain only as the cited regression identity.
- The forward source is the canonical decision-05 D3 Frechet source over the B2a BdG-mode response; the B2a overlaps `c_j` carry the chi coupling. It is not a separate live B1-chi adapter-overlap path.

## Files Created Or Modified

- `software/stage1_solver/src/stage1_solver/patha_b2b_maxwell.py`
- `software/stage1_solver/mathematica/mt15_03_spike1_vsh_maxwell_operator.wls`
- `software/stage1_solver/mathematica/mt15_04_spike2_transfer_n0.wls`
- `software/stage1_solver/mathematica/mt15_03_patha_b2b_vsh_maxwell_operator.wls`
- `software/stage1_solver/mathematica/mt15_04_patha_b2b_maxwell_transfer.wls`
- `software/stage1_solver/directives/pathA_11_chunk_b2b_maxwell_transfer.md`
- `software/stage1_solver/tests/test_patha_b2b_maxwell.py`
- `software/stage1_solver/reports/patha_b2b_maxwell_transfer_derivation_report.md`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/bundles/patha_b2b_validated_maxwell_transfer_tau_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_11_nw_11_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_13_nw_11_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_13_nw_13_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_16_nw_13_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_19_nw_15_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_22_nw_17_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_27_nw_17_w_0.028_rs_3.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_31_nw_13_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_37_nw_17_w_0.028_rs_4.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_39_nw_15_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_43_nw_16_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_43_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_47_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_47_nw_17_w_0.036_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_47_nw_17_w_0.046_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_47_nw_17_w_0.058_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_7_nw_7_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_1_nr_9_nw_9_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_2_nr_13_nw_13_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_2_nr_22_nw_17_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/independent/patha_b2b_independent_tau_2_nr_47_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_5_nw_5_w_0.036_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_6_nw_6_w_0.036_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_7_nw_7_w_0.036_rs_0.8.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_7_nw_7_w_0.036_rs_0.9.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_7_nw_7_w_0.036_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_7_nw_7_w_0.046_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_1_nr_7_nw_7_w_0.058_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/mathematica/patha_b2b_mma_tau_2_nr_7_nw_7_w_0.036_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/operator/patha_b2b_spike1_operator_tau_1_nr_7_nw_7.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_11_nw_11_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_11_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.028_rs_1.25.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.028_rs_1.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.028_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.036_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.046_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_13_nw_13_w_0.058_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_16_nw_13_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_19_nw_15_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_2.25.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_2.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_3.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_4.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.036_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.046_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_22_nw_17_w_0.058_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_27_nw_11_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_27_nw_17_w_0.028_rs_3.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_31_nw_13_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_37_nw_13_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_37_nw_17_w_0.028_rs_4.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_39_nw_15_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_43_nw_16_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_43_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_47_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_47_nw_17_w_0.036_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_47_nw_17_w_0.046_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_47_nw_17_w_0.058_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_57_nw_17_w_0.028_rs_6.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_7_nw_7_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_7_nw_7_w_0.036_rs_1.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_1_nr_9_nw_9_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_2_nr_13_nw_13_w_0.028_rs_1.75.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_2_nr_22_nw_17_w_0.028_rs_2.5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_2_nr_47_nw_17_w_0.028_rs_5.json`
- `software/stage1_solver/runs/patha_b2b_maxwell_transfer/python/patha_b2b_python_tau_2_nr_7_nw_7_w_0.036_rs_1.json`
