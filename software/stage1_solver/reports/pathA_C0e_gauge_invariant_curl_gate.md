# Path-A C0e Gauge-Invariant Curl Gate

Primary C0e-4 verdict: **GAUGE_FRAMING_REFUTED**

## Scope And Artifact Policy

Diagnosis only. C0e uses saved C0b states and assembled Jacobians, recomputes fresh dense SVD modes, runs one read-only `J*dx=-F` solve per labeled artifact, and evaluates trial residuals only on temporary tensors. It implements no deflation, gauge fix, stencil change, xi reassembly, or recrawl.

```yaml
primary_artifact: {'tau': 0.02899, 'tau_source_label': 'deepest_saved_stall_attempt', 'state_artifact': '/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/states/attempt_001_tau_0p02899.npz', 'matrix_path': '/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/matrices/attempt_tau_0p02899_bt_0.npz'}
artifact_policy: {'same_state_and_matrix_used_for_c0e_0_through_c0e_3': True, 'preferred_primary': 'deepest saved stall/source closest to tau≈0.029 with assembled J', 'both_tau_sources_labeled_when_available': True, 'available_artifact_count': 2}
scope_guard: {'diagnosis_only': True, 'single_read_only_linear_solve_per_artifact': True, 'trial_residual_evaluations_on_temporary_tensors_only': True, 'state_advance_or_newton_write': False, 'gauge_fix_or_deflation_implemented': False, 'stencil_change_implemented': False, 'recrawl_implemented': False, 'changed_xi_reassembly': False, 'faithful_operators_touched_by_c0e': False, 'frozen_physics_touched_by_c0e': False, 'physical_export_guard_touched_by_c0e': False}
faithful_operator_boundary: {'physical_residual_entry_point': 'stage1_solver.coupled_branch.patha_closed_branch_residual', 'operators_modified_by_c0': False, 'frozen_physics_modified_by_c0': False, 'physical_export_guard_modified_by_c0': False}
```

## Operator Sources

```yaml
curl_fraction: {'field_strength': 'stage1_solver.operators.tensor_center_gradient_r(aw) - tensor_center_gradient_w(ar)', 'localization_weight': 'stage1_solver.coupled_branch.localization_weight_torch(grid.w_centers, branch)', 'ratio': '||Z*F_rw[v]|| / ||A[v]||', 'norms': ['unweighted_code_space', 'cell_volume_weighted']}
coupled_generator: {'real_lane_formula': 'delta psi_real=-(q/hbar)*lambda*psi_imag; delta psi_imag=+(q/hbar)*lambda*psi_real; delta a0=0; delta ar=d_r(lambda); delta aw=d_w(lambda)', 'q_over_hbar_source': 'branch.gauge_charge / branch.hbar from stage1_solver.coupled_branch.coupled_pde_residual alpha and q*A0 matter coupling', 'q_over_hbar_value': 0.35}
maxwell_operator: {'source': 'stage1_solver.operators.localized_maxwell_operator', 'curl_rows': 'ar=-d_w(Z*F_rw), aw=radial_divergence_from_center_flux(Z*F_rw)', 'gauge_penalty_rows': '(1/xi)*grad(Z*axisymmetric_vector_divergence(A))', 'xi': 1.0, 'sponge_gauge_strength': 0.0}
```

## Artifact Sources

| tau | source | converged | state | matrix | verdict |
| --- | --- | --- | --- | --- | --- |
| 2.899000e-02 | deepest_saved_stall_attempt | false | /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/states/attempt_001_tau_0p02899.npz | /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/matrices/attempt_tau_0p02899_bt_0.npz | GAUGE_FRAMING_REFUTED |
| 3.000000e-02 | converged_tau_source | true | /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/states/attempt_000_tau_0p03.npz | /var/projects/toy_physics/software/stage1_solver/runs/pathA_C0b_wall_diagnosis/matrices/attempt_tau_0p03_bt_0.npz | GAUGE_FRAMING_REFUTED |

## Artifact tau=2.899000000000e-02 (deepest_saved_stall_attempt)

### C0e-0 Read-Only Newton Framing

```yaml
solve_method: dense_np_linalg_solve
initial_residual_l2: 0.00023760662696478132
initial_residual_linf: 5.322040363520651e-05
full_step_residual_l2_ratio: 5.079172580992675
full_step_residual_linf_ratio: 8.018483255583421
gauge_removed_step_residual_l2_ratio: 5.079172580992675
gauge_removed_step_residual_linf_ratio: 8.018483255583421
step_norm: 0.10258689022246326
near_null_component_norm: 4.356036079167163e-12
complement_component_norm: 0.10258689022246326
near_null_component_fraction: 1.803014446021408e-21
linear_solve_relative_residual: 7.906368184950328e-13
trial_residual_eval_scope: temporary_tensors_only_no_state_advance_no_line_search_no_newton_write
```

### Controls

```yaml
status: PASS
reason: None
bands: {'unweighted': {'status': 'SEPARATED', 'positive_min': 8.914242170483634e-16, 'positive_max': 1.292970967505837e-15, 'negative_min': 8.278057249635262, 'negative_max': 15.66630495789788, 'geometric_separator': 1.034566947622512e-07, 'separation_ratio': 6402353539000010.0, 'separation_log10': 15.806339652281952}, 'cell_volume_weighted': {'status': 'SEPARATED', 'positive_min': 1.074390273260633e-15, 'positive_max': 1.391812884137857e-15, 'negative_min': 8.669247827457557, 'negative_max': 17.917239778515754, 'geometric_separator': 1.0984521301376566e-07, 'separation_ratio': 6228745204372516.0, 'separation_log10': 15.794400565807795}}
negative_capture_max_observed: 0.45603900781055506
```
| name | family | curl_unw | curl_w | P_G | P_cpl | non_A | construction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| smooth_sin_poly_gradient | positive_gradient | 8.914242e-16 | 1.074390e-15 | 1.000000e+00 | 1.000000e+00 | 2.220446e-16 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| smooth_mixed_gradient | positive_gradient | 1.247053e-15 | 1.391813e-15 | 1.000000e+00 | 1.000000e+00 | 0.000000e+00 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| high_k_sine_gradient | positive_gradient | 1.240151e-15 | 1.284269e-15 | 1.000000e+00 | 1.000000e+00 | 0.000000e+00 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| checkerboard_gradient | positive_gradient | 1.292971e-15 | 1.381820e-15 | 1.000000e+00 | 1.000000e+00 | 0.000000e+00 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| smooth_stream_sin | negative_stream_function_transverse | 8.278057e+00 | 8.669248e+00 | 3.212945e-01 | 3.212684e-01 | 0.000000e+00 | independent_stream_function_A_not_random_minus_own_gradient_projection |
| smooth_stream_mixed | negative_stream_function_transverse | 8.503628e+00 | 1.010160e+01 | 4.560390e-01 | 4.559988e-01 | 1.110223e-16 | independent_stream_function_A_not_random_minus_own_gradient_projection |
| high_k_stream_sine | negative_stream_function_transverse | 1.477800e+01 | 1.615710e+01 | 3.186413e-01 | 3.186200e-01 | 0.000000e+00 | independent_stream_function_A_not_random_minus_own_gradient_projection |
| checkerboard_stream | negative_stream_function_transverse | 1.566630e+01 | 1.791724e+01 | 1.006014e-01 | 1.006007e-01 | 2.220446e-16 | independent_stream_function_A_not_random_minus_own_gradient_projection |

### Phase Mode Separate Gate

```yaml
mode_index: 0
sigma: 1.43587972095615e-12
classification: PHASE_GAUGE_CONFIRMED_BY_PHASE_AND_COUPLED_CAPTURE
phase_capture_fraction: 0.9941818585990793
coupled_capture_fraction: 0.9999865144306872
spatial_a_norm_unweighted: 0.07627646175430372
z_f_rw_norm_unweighted: 0.012603332949842101
spatial_a_norm_cell_volume_weighted: 0.027627166091042203
z_f_rw_norm_cell_volume_weighted: 0.005948691935595193
note: The curl ratio is denominator-noise for the phase mode and is not used as its gate.
```

### C0e-1 Curl Gate For Maxwell-Lane Modes

| mode | sigma | curl_unw | curl_w | outcome | margin_log10 | P_G | P_cpl | A_energy | ar | aw |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.720320e-11 | 9.325885e-02 | 1.257907e-01 | TRANSVERSE | 5.954931e+00 | 9.938824e-01 | 9.991596e-01 | 9.945887e-01 | 8.738668e-03 | 9.858501e-01 |
| 2 | 2.569643e-10 | 1.441708e+00 | 1.442006e+00 | TRANSVERSE | 7.118186e+00 | 8.299733e-01 | 8.301451e-01 | 9.998231e-01 | 9.551637e-01 | 4.465939e-02 |
| 3 | 5.917734e-10 | 5.293367e-02 | 5.365212e-02 | TRANSVERSE | 5.688806e+00 | 9.995554e-01 | 9.997692e-01 | 9.997856e-01 | 9.864249e-01 | 1.336071e-02 |
| 4 | 2.808634e-08 | 1.085534e-02 | 1.009633e-02 | TRANSVERSE | 4.963382e+00 | 9.989212e-01 | 9.989316e-01 | 9.999850e-01 | 2.225656e-01 | 7.774194e-01 |

Boundary/interior curl split and norm-specific margins are in the machine JSON for each mode.

### C0e-3 Mechanism Evidence

```yaml
primary_mechanism: SMOOTH_K2
secondary_flags: ['NONZERO_UNEXPLAINED_JV_REMAINDER', 'LARGE_COMPONENT_CANCELLATION_RELATIVE_TO_NEAR_NULL_JV']
lambda_frequency_counts: {'SMOOTH': 3, 'MIXED': 1}
row_scaling_ratio: 0.5514997971090752
max_unexplained_fraction_of_component_sum_norm: 0.0007138950947225955
max_gauge_penalty_to_assembled_jv_ratio: 18228180612.447304
```

Assembled Maxwell-row scaling:

```yaml
{'row_norm_kind': 'sparse_row_abs_max', 'groups': {'matter_psi_rows': {'count': 512, 'min': 135.5785391369993, 'median': 173.17238607306473, 'max': 226.73060892055955, 'rms': 176.83330216592032}, 'gauss_a0_rows': {'count': 256, 'min': 150.71097934758558, 'median': 351.26474931868756, 'max': 493.14997163667454, 'rms': 347.66753998397354}, 'maxwell_ar_rows': {'count': 256, 'min': 80.40302539598561, 'median': 97.55902941100211, 'max': 334.7181512267305, 'rms': 134.55349200642257}, 'maxwell_aw_rows': {'count': 256, 'min': 78.32186920836185, 'median': 92.91656294114222, 'max': 369.99931293464823, 'rms': 135.20537184674174}, 'wall_rows': {'count': 16, 'min': 2.1684265887509127, 'median': 4.316379053459011, 'max': 10.307048748295095, 'rms': 4.823322913428753}, 'mass_row': {'count': 1, 'min': 0.014800261631798196, 'median': 0.014800261631798196, 'max': 0.014800261631798196, 'rms': 0.014800261631798196}}, 'maxwell_araw_rms': 134.87982574747443, 'non_maxwell_rms': 244.56913031429096, 'maxwell_to_non_maxwell_rms_ratio': 0.5514997971090752}
```

| mode | Jv | curl | penalty | current | matter | gauss | wall_mass | unexplained | lambda_class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.721116e-11 | 3.021108e-01 | 3.026399e-01 | 6.771445e-04 | 1.089083e-13 | 3.574279e-14 | 2.513849e-14 | 3.370443e-14 | SMOOTH |
| 2 | 2.569656e-10 | 4.684035e+00 | 4.684015e+00 | 1.624648e-04 | 1.164439e-13 | 2.241502e-14 | 7.739548e-14 | 8.338920e-14 | MIXED |
| 3 | 5.917823e-10 | 1.719825e-01 | 1.719409e-01 | 6.856761e-05 | 9.306852e-14 | 1.675892e-14 | 1.996418e-14 | 3.848071e-14 | SMOOTH |
| 4 | 2.808633e-08 | 3.545073e-02 | 3.543322e-02 | 3.079045e-04 | 2.686157e-13 | 3.178224e-14 | 3.080415e-14 | 4.122096e-14 | SMOOTH |

### C0e-4 Verdict And C0f Recommendation

```yaml
verdict: GAUGE_FRAMING_REFUTED
verdict_support: {'controls_status': 'PASS', 'control_bands': {'unweighted': {'status': 'SEPARATED', 'positive_min': 8.914242170483634e-16, 'positive_max': 1.292970967505837e-15, 'negative_min': 8.278057249635262, 'negative_max': 15.66630495789788, 'geometric_separator': 1.034566947622512e-07, 'separation_ratio': 6402353539000010.0, 'separation_log10': 15.806339652281952}, 'cell_volume_weighted': {'status': 'SEPARATED', 'positive_min': 1.074390273260633e-15, 'positive_max': 1.391812884137857e-15, 'negative_min': 8.669247827457557, 'negative_max': 17.917239778515754, 'geometric_separator': 1.0984521301376566e-07, 'separation_ratio': 6228745204372516.0, 'separation_log10': 15.794400565807795}}, 'mode_outcomes': {1: 'TRANSVERSE', 2: 'TRANSVERSE', 3: 'TRANSVERSE', 4: 'TRANSVERSE'}, 'mode_curl_support': [{'mode_index': 1, 'sigma': 4.720320480768486e-11, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 0.09325885396184329, 'weighted_curl_fraction': 0.12579074706184312, 'band_margin_log10_min': 5.954931474758163, 'a_only_p_g_fraction': 0.9938823913577117, 'coupled_capture_fraction': 0.9991595729239308}, {'mode_index': 2, 'sigma': 2.569643078156986e-10, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 1.4417077614825824, 'weighted_curl_fraction': 1.4420062852022166, 'band_margin_log10_min': 7.118186017940804, 'a_only_p_g_fraction': 0.8299733285268662, 'coupled_capture_fraction': 0.8301450585484521}, {'mode_index': 3, 'sigma': 5.917734218459453e-10, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 0.052933665716501946, 'weighted_curl_fraction': 0.05365212130303726, 'band_margin_log10_min': 5.6888057624360675, 'a_only_p_g_fraction': 0.9995554478155377, 'coupled_capture_fraction': 0.9997691748835861}, {'mode_index': 4, 'sigma': 2.8086338422231964e-08, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 0.010855343855719549, 'weighted_curl_fraction': 0.010096327914675012, 'band_margin_log10_min': 4.963382312024457, 'a_only_p_g_fraction': 0.9989211771219381, 'coupled_capture_fraction': 0.9989315894944502}], 'phase_mode': {'mode_index': 0, 'phase_capture_fraction': 0.9941818585990793, 'coupled_capture_fraction': 0.9999865144306872}, 'mechanism': {'primary_mechanism': 'SMOOTH_K2', 'secondary_flags': ['NONZERO_UNEXPLAINED_JV_REMAINDER', 'LARGE_COMPONENT_CANCELLATION_RELATIVE_TO_NEAR_NULL_JV'], 'lambda_frequency_counts': {'SMOOTH': 3, 'MIXED': 1}, 'row_scaling_ratio': 0.5514997971090752, 'max_unexplained_fraction_of_component_sum_norm': 0.0007138950947225955, 'max_gauge_penalty_to_assembled_jv_ratio': 18228180612.447304, 'evidence_note': 'ODD_EVEN_DECOUPLING is emitted only when the lambda-preimage frequency evidence is checkerboard/high-k and the row/block residual budget is present.'}, 'required_maxwell_mode_count': 4, 'observed_maxwell_mode_count': 4, 'mode_2_outcome': 'TRANSVERSE', 'reason': 'one_of_modes_1_3_4_has_real_curl'}
```

Stop before C0f: the gauge reading is refuted by curl-carrying modes among 1/3/4.

## Artifact tau=3.000000000000e-02 (converged_tau_source)

### C0e-0 Read-Only Newton Framing

```yaml
solve_method: dense_np_linalg_solve
initial_residual_l2: 1.4697080797396755e-09
initial_residual_linf: 3.3879743455145217e-10
full_step_residual_l2_ratio: 0.24701723069563955
full_step_residual_linf_ratio: 0.36687658891223535
gauge_removed_step_residual_l2_ratio: 0.24701723069563955
gauge_removed_step_residual_linf_ratio: 0.36687658891223535
step_norm: 3.843088140138664e-08
near_null_component_norm: 3.714852760156458e-20
complement_component_norm: 3.843088140138664e-08
near_null_component_fraction: 9.343778183571955e-25
linear_solve_relative_residual: 7.032364967260572e-14
trial_residual_eval_scope: temporary_tensors_only_no_state_advance_no_line_search_no_newton_write
```

### Controls

```yaml
status: PASS
reason: None
bands: {'unweighted': {'status': 'SEPARATED', 'positive_min': 8.914242170483634e-16, 'positive_max': 1.292970967505837e-15, 'negative_min': 8.278057249635262, 'negative_max': 15.66630495789788, 'geometric_separator': 1.034566947622512e-07, 'separation_ratio': 6402353539000010.0, 'separation_log10': 15.806339652281952}, 'cell_volume_weighted': {'status': 'SEPARATED', 'positive_min': 1.074390273260633e-15, 'positive_max': 1.391812884137857e-15, 'negative_min': 8.669247827457557, 'negative_max': 17.917239778515754, 'geometric_separator': 1.0984521301376566e-07, 'separation_ratio': 6228745204372516.0, 'separation_log10': 15.794400565807795}}
negative_capture_max_observed: 0.45603900781055506
```
| name | family | curl_unw | curl_w | P_G | P_cpl | non_A | construction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| smooth_sin_poly_gradient | positive_gradient | 8.914242e-16 | 1.074390e-15 | 1.000000e+00 | 1.000000e+00 | 2.220446e-16 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| smooth_mixed_gradient | positive_gradient | 1.247053e-15 | 1.391813e-15 | 1.000000e+00 | 1.000000e+00 | 0.000000e+00 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| high_k_sine_gradient | positive_gradient | 1.240151e-15 | 1.284269e-15 | 1.000000e+00 | 1.000000e+00 | 0.000000e+00 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| checkerboard_gradient | positive_gradient | 1.292971e-15 | 1.381820e-15 | 1.000000e+00 | 1.000000e+00 | 0.000000e+00 | held_out_discrete_gradient_not_used_to_define_mode_classification |
| smooth_stream_sin | negative_stream_function_transverse | 8.278057e+00 | 8.669248e+00 | 3.212945e-01 | 3.212682e-01 | 0.000000e+00 | independent_stream_function_A_not_random_minus_own_gradient_projection |
| smooth_stream_mixed | negative_stream_function_transverse | 8.503628e+00 | 1.010160e+01 | 4.560390e-01 | 4.559984e-01 | 1.110223e-16 | independent_stream_function_A_not_random_minus_own_gradient_projection |
| high_k_stream_sine | negative_stream_function_transverse | 1.477800e+01 | 1.615710e+01 | 3.186413e-01 | 3.186199e-01 | 0.000000e+00 | independent_stream_function_A_not_random_minus_own_gradient_projection |
| checkerboard_stream | negative_stream_function_transverse | 1.566630e+01 | 1.791724e+01 | 1.006014e-01 | 1.006007e-01 | 2.220446e-16 | independent_stream_function_A_not_random_minus_own_gradient_projection |

### Phase Mode Separate Gate

```yaml
mode_index: 0
sigma: 7.434693840833268e-14
classification: PHASE_GAUGE_CONFIRMED_BY_PHASE_AND_COUPLED_CAPTURE
phase_capture_fraction: 0.9999999938582824
coupled_capture_fraction: 0.9999999999992311
spatial_a_norm_unweighted: 7.836902765234002e-05
z_f_rw_norm_unweighted: 2.481051695832616e-07
spatial_a_norm_cell_volume_weighted: 2.8046196560794556e-05
z_f_rw_norm_cell_volume_weighted: 1.1730622921637655e-07
note: The curl ratio is denominator-noise for the phase mode and is not used as its gate.
```

### C0e-1 Curl Gate For Maxwell-Lane Modes

| mode | sigma | curl_unw | curl_w | outcome | margin_log10 | P_G | P_cpl | A_energy | ar | aw |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.700699e-11 | 7.703935e-02 | 1.037776e-01 | TRANSVERSE | 5.871954e+00 | 9.995140e-01 | 9.993824e-01 | 9.999985e-01 | 7.729896e-03 | 9.922686e-01 |
| 2 | 2.622571e-10 | 1.443001e+00 | 1.444095e+00 | TRANSVERSE | 7.118815e+00 | 8.298041e-01 | 8.298119e-01 | 9.999865e-01 | 9.562146e-01 | 4.377191e-02 |
| 3 | 5.926281e-10 | 3.755573e-02 | 3.806719e-02 | TRANSVERSE | 5.539770e+00 | 9.998625e-01 | 9.998828e-01 | 9.999790e-01 | 9.864547e-01 | 1.352430e-02 |
| 4 | 2.811268e-08 | 1.231383e-02 | 1.145779e-02 | TRANSVERSE | 5.018320e+00 | 9.989308e-01 | 9.989298e-01 | 9.999964e-01 | 2.223673e-01 | 7.776291e-01 |

Boundary/interior curl split and norm-specific margins are in the machine JSON for each mode.

### C0e-3 Mechanism Evidence

```yaml
primary_mechanism: SMOOTH_K2
secondary_flags: ['NONZERO_UNEXPLAINED_JV_REMAINDER', 'LARGE_COMPONENT_CANCELLATION_RELATIVE_TO_NEAR_NULL_JV']
lambda_frequency_counts: {'SMOOTH': 3, 'MIXED': 1}
row_scaling_ratio: 0.5514874810907664
max_unexplained_fraction_of_component_sum_norm: 0.000596234840114065
max_gauge_penalty_to_assembled_jv_ratio: 17877834910.56677
```

Assembled Maxwell-row scaling:

```yaml
{'row_norm_kind': 'sparse_row_abs_max', 'groups': {'matter_psi_rows': {'count': 512, 'min': 135.56797209758844, 'median': 173.17125663665433, 'max': 226.9090127663893, 'rms': 176.8447380429161}, 'gauss_a0_rows': {'count': 256, 'min': 150.71097934758558, 'median': 351.26474931868756, 'max': 493.14997163667454, 'rms': 347.66753998397354}, 'maxwell_ar_rows': {'count': 256, 'min': 80.40302497876003, 'median': 97.55902759473614, 'max': 334.7181457071445, 'rms': 134.55348991463302}, 'maxwell_aw_rows': {'count': 256, 'min': 78.32186920836185, 'median': 92.91656295430431, 'max': 369.9993074150622, 'rms': 135.20536980047694}, 'wall_rows': {'count': 16, 'min': 2.2439737034331624, 'median': 4.4550322452931335, 'max': 10.66633333459402, 'rms': 4.9822083027324995}, 'mass_row': {'count': 1, 'min': 0.014899116231355477, 'median': 0.014899116231355477, 'max': 0.014899116231355477, 'rms': 0.014899116231355477}}, 'maxwell_araw_rms': 134.8798236785083, 'non_maxwell_rms': 244.57458836914043, 'maxwell_to_non_maxwell_rms_ratio': 0.5514874810907664}
```

| mode | Jv | curl | penalty | current | matter | gauss | wall_mass | unexplained | lambda_class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 4.700842e-11 | 2.502299e-01 | 2.507591e-01 | 6.784216e-04 | 1.088284e-13 | 3.209766e-14 | 1.178562e-14 | 2.802851e-14 | SMOOTH |
| 2 | 2.622579e-10 | 4.688619e+00 | 4.688604e+00 | 1.605666e-04 | 1.332031e-13 | 2.792875e-14 | 1.297356e-14 | 3.705641e-14 | MIXED |
| 3 | 5.926320e-10 | 1.220333e-01 | 1.219909e-01 | 6.924915e-05 | 1.036178e-13 | 2.642352e-14 | 1.699012e-14 | 3.211352e-14 | SMOOTH |
| 4 | 2.811267e-08 | 4.017326e-02 | 4.016159e-02 | 3.093763e-04 | 2.594196e-13 | 3.262459e-14 | 1.219110e-14 | 4.040935e-14 | SMOOTH |

### C0e-4 Verdict And C0f Recommendation

```yaml
verdict: GAUGE_FRAMING_REFUTED
verdict_support: {'controls_status': 'PASS', 'control_bands': {'unweighted': {'status': 'SEPARATED', 'positive_min': 8.914242170483634e-16, 'positive_max': 1.292970967505837e-15, 'negative_min': 8.278057249635262, 'negative_max': 15.66630495789788, 'geometric_separator': 1.034566947622512e-07, 'separation_ratio': 6402353539000010.0, 'separation_log10': 15.806339652281952}, 'cell_volume_weighted': {'status': 'SEPARATED', 'positive_min': 1.074390273260633e-15, 'positive_max': 1.391812884137857e-15, 'negative_min': 8.669247827457557, 'negative_max': 17.917239778515754, 'geometric_separator': 1.0984521301376566e-07, 'separation_ratio': 6228745204372516.0, 'separation_log10': 15.794400565807795}}, 'mode_outcomes': {1: 'TRANSVERSE', 2: 'TRANSVERSE', 3: 'TRANSVERSE', 4: 'TRANSVERSE'}, 'mode_curl_support': [{'mode_index': 1, 'sigma': 4.7006989691328044e-11, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 0.07703935140343174, 'weighted_curl_fraction': 0.10377756471239048, 'band_margin_log10_min': 5.871954018352952, 'a_only_p_g_fraction': 0.9995139675622529, 'coupled_capture_fraction': 0.9993823808912822}, {'mode_index': 2, 'sigma': 2.6225709826789087e-10, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 1.4430011623033967, 'weighted_curl_fraction': 1.444094795377305, 'band_margin_log10_min': 7.118814567376489, 'a_only_p_g_fraction': 0.829804081276686, 'coupled_capture_fraction': 0.8298119432852828}, {'mode_index': 3, 'sigma': 5.926280559311367e-10, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 0.03755573117334073, 'weighted_curl_fraction': 0.03806719219672659, 'band_margin_log10_min': 5.539769709424326, 'a_only_p_g_fraction': 0.9998624973765375, 'coupled_capture_fraction': 0.999882807964325}, {'mode_index': 4, 'sigma': 2.8112679473453323e-08, 'outcome': 'TRANSVERSE', 'unweighted_curl_fraction': 0.012313826924187777, 'weighted_curl_fraction': 0.011457792720669381, 'band_margin_log10_min': 5.018319825914713, 'a_only_p_g_fraction': 0.9989308103957472, 'coupled_capture_fraction': 0.9989297855681731}], 'phase_mode': {'mode_index': 0, 'phase_capture_fraction': 0.9999999938582824, 'coupled_capture_fraction': 0.9999999999992311}, 'mechanism': {'primary_mechanism': 'SMOOTH_K2', 'secondary_flags': ['NONZERO_UNEXPLAINED_JV_REMAINDER', 'LARGE_COMPONENT_CANCELLATION_RELATIVE_TO_NEAR_NULL_JV'], 'lambda_frequency_counts': {'SMOOTH': 3, 'MIXED': 1}, 'row_scaling_ratio': 0.5514874810907664, 'max_unexplained_fraction_of_component_sum_norm': 0.000596234840114065, 'max_gauge_penalty_to_assembled_jv_ratio': 17877834910.56677, 'evidence_note': 'ODD_EVEN_DECOUPLING is emitted only when the lambda-preimage frequency evidence is checkerboard/high-k and the row/block residual budget is present.'}, 'required_maxwell_mode_count': 4, 'observed_maxwell_mode_count': 4, 'mode_2_outcome': 'TRANSVERSE', 'reason': 'one_of_modes_1_3_4_has_real_curl'}
```

Stop before C0f: the gauge reading is refuted by curl-carrying modes among 1/3/4.

## Machine Artifact

`/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0e_curl_gate/pathA_C0e_gauge_invariant_curl_gate.json`

