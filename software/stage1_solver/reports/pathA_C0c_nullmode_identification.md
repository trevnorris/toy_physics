# Path-A C0c Nullmode Identification

Verdict: **MIXED**

## Scope

Diagnosis only: generators are evaluated on existing C0b states and saved Jacobians. `coupled_branch.py`, `operators.py`, frozen physics, and `physical_export_permitted` are untouched.

## Generator Inventory

| name | class | status | norm | description |
| --- | --- | --- | --- | --- |
| phase | GAUGE_PHASE | EXACT_SYMMETRY | 1.106963e+00 | global U(1) phase: (-psi_imag, psi_real, 0, 0, 0, 0, 0) |
| translation_r | TRANSLATION | PROBE | 1.109158e+00 | broken radial-translation probe: finite-difference d_r(state) |
| translation_w | TRANSLATION | PROBE | 3.325876e+00 | broken axial-translation probe: finite-difference d_w(state) |
| dilation_r | DILATION | PROBE | 1.150248e+00 | broken radial dilation probe: r*d_r(field lanes) |
| maxwell_residual_gauge | MAXWELL_RESIDUAL_GAUGE | PROBE | 1.182943e+01 | Maxwell-sector pure-gradient probe: delta A=(0, d_r lambda, d_w lambda) |

## Annihilation And Equivariance

### tau=3.000000000000e-02 (converged=True, residual=3.387974e-10)
| generator | status | test | assembled | jvp | crosscheck | equiv_jvp | gate | reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| phase | MEASURED | root_annihilation | 1.485839e-12 | 1.485829e-12 | 2.807855e-17 | - | true | - |
| translation_r | MEASURED | root_annihilation | 3.033432e-02 | 3.033182e-02 | 2.154307e-04 | - | false | - |
| translation_w | MEASURED | root_annihilation | 3.783344e-03 | 3.783343e-03 | 1.599927e-05 | - | false | - |
| dilation_r | MEASURED | root_annihilation | 4.431382e-02 | 4.431274e-02 | 1.350332e-04 | - | false | - |
| maxwell_residual_gauge | MEASURED | root_annihilation | 4.508338e-03 | 4.508338e-03 | 3.643549e-17 | - | false | - |

### tau=2.899000000000e-02 (converged=False, residual=5.322040e-05)
| generator | status | test | assembled | jvp | crosscheck | equiv_jvp | gate | reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| phase | MEASURED | nonroot_phase_equivariance_identity | 2.098478e-07 | 2.098478e-07 | 2.736154e-17 | 0.000000e+00 | - | - |
| translation_r | NOT_MEASURED | not_measured_nonroot_probe | - | - | - | - | - | nonconverged_state_and_probe_not_an_exact_symmetry |
| translation_w | NOT_MEASURED | not_measured_nonroot_probe | - | - | - | - | - | nonconverged_state_and_probe_not_an_exact_symmetry |
| dilation_r | NOT_MEASURED | not_measured_nonroot_probe | - | - | - | - | - | nonconverged_state_and_probe_not_an_exact_symmetry |
| maxwell_residual_gauge | NOT_MEASURED | not_measured_nonroot_probe | - | - | - | - | - | nonconverged_state_and_probe_not_an_exact_symmetry |

## Near-Null SVD Overlap

### tau=3.000000000000e-02

| mode | sigma | class | phase | tr | tw | dil | maxwell | residual | psi_re | psi_im | a0 | ar | aw | r0 | mu |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 7.434694e-14 | GAUGE_PHASE | 1.000000e+00 | 3.492718e-14 | 3.935185e-13 | 2.639141e-14 | 1.511049e-05 | 5.913390e-09 | 2.836855e-27 | 1.000000e+00 | 1.567147e-28 | 2.949980e-10 | 5.846706e-09 | 7.273224e-25 | 4.583897e-27 |
| 1 | 4.700699e-11 | UNEXPLAINED_STIFFNESS | 7.737668e-05 | 1.902133e-14 | 1.947792e-13 | 1.449331e-14 | 1.111623e-01 | 9.876429e-01 | 8.029988e-28 | 1.499120e-06 | 1.841054e-30 | 7.729896e-03 | 9.922686e-01 | 1.809274e-25 | 1.556852e-27 |
| 2 | 2.622571e-10 | UNEXPLAINED_STIFFNESS | 4.264932e-06 | 7.581350e-15 | 8.031111e-14 | 5.751321e-15 | 9.718889e-03 | 9.999055e-01 | 1.384099e-28 | 1.351818e-05 | 9.206175e-30 | 9.562146e-01 | 4.377191e-02 | 2.945885e-26 | 1.380658e-28 |
| 3 | 5.926281e-10 | UNEXPLAINED_STIFFNESS | 1.167517e-05 | 6.590076e-15 | 6.528069e-14 | 5.187538e-15 | 5.710172e-01 | 6.739393e-01 | 8.124176e-29 | 2.101846e-05 | 1.015027e-30 | 9.864547e-01 | 1.352430e-02 | 1.430914e-26 | 3.586010e-28 |
| 4 | 2.811268e-08 | UNEXPLAINED_STIFFNESS | 2.610379e-07 | 4.299256e-15 | 4.153583e-14 | 3.391606e-15 | 7.709786e-01 | 4.055920e-01 | 3.805748e-29 | 3.625923e-06 | 1.049565e-29 | 2.223673e-01 | 7.776291e-01 | 5.754529e-27 | 6.104839e-29 |

### tau=2.899000000000e-02

| mode | sigma | class | phase | tr | tw | dil | maxwell | residual | psi_re | psi_im | a0 | ar | aw | r0 | mu |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1.435880e-12 | - | 9.970867e-01 | 1.963854e-13 | 1.664590e-11 | 5.843922e-15 | 3.763009e-03 | 5.803981e-03 | 1.865632e-24 | 9.941819e-01 | 2.036160e-28 | 1.430328e-04 | 5.675066e-03 | 1.548389e-21 | 7.582747e-24 |
| 1 | 4.720320e-11 | - | 7.355168e-02 | 4.888012e-15 | 1.758202e-12 | 9.671879e-15 | 1.047512e-01 | 9.836173e-01 | 1.778127e-26 | 5.411280e-03 | 4.687111e-30 | 8.738668e-03 | 9.858501e-01 | 1.716524e-23 | 8.056863e-26 |
| 2 | 2.569643e-10 | - | 1.278037e-02 | 6.035266e-14 | 5.971658e-12 | 4.365838e-15 | 1.767383e-02 | 9.995243e-01 | 2.321633e-25 | 1.768633e-04 | 4.748169e-30 | 9.551637e-01 | 4.465939e-02 | 1.994851e-22 | 9.725691e-25 |
| 3 | 5.917734e-10 | - | 1.390570e-02 | 1.695008e-14 | 1.528408e-12 | 6.973613e-17 | 5.714848e-01 | 6.732117e-01 | 1.543747e-26 | 2.143707e-04 | 4.797291e-29 | 9.864249e-01 | 1.336071e-02 | 1.306253e-23 | 6.242999e-26 |
| 4 | 2.808634e-08 | - | 3.367799e-03 | 2.460271e-14 | 1.576909e-12 | 4.532633e-15 | 7.711262e-01 | 4.053530e-01 | 1.799315e-26 | 1.496167e-05 | 5.718282e-30 | 2.225656e-01 | 7.774194e-01 | 1.403918e-23 | 7.368619e-26 |

## Dense Sigma Validation

```yaml
status: MATCH
tau: 0.03
method: dense_full_jvp_jacobian_svd
assembled_matrix_method: dense_svd_recomputed_from_saved_c0b_matrix
dense_full_jvp_sigma_min: 3.7817308463154297e-16
dense_full_jvp_sigma_max: 825.9943942481186
assembled_sigma_min: 2.064616665800032e-15
assembled_sigma_max: 825.9943942481183
abs_difference: 1.6864435811684892e-15
tolerance_abs: 1.8340759894111405e-10
trusted_sigma_source: both
elapsed_seconds: 8.428366545005701
chunk_size: 32
```

## Verdict Support

```yaml
thresholds: {'annihilation_rel_max': 1e-08, 'overlap_min': 0.9, 'span_residual_fraction_reference': 0.2}
generator_null_gates: {'phase': True, 'translation_r': False, 'translation_w': False, 'dilation_r': False, 'maxwell_residual_gauge': False}
explained_mode_count: 1
cluster_mode_count: 5
mode_classifications:
  - mode_index: 0, sigma: 7.434693840833268e-14, classification: GAUGE_PHASE, residual: 5.913390443978983e-09, lanes: {'psi_real': 2.83685464103925e-27, 'psi_imag': 0.9999999938582957, 'a0': 1.5671470087636524e-28, 'ar': 2.949979981408009e-10, 'aw': 5.8467064970324326e-09, 'r0': 7.273224323143217e-25, 'mu': 4.583896994797658e-27}, support: {'generator': 'phase', 'overlap': 0.9999999969291412, 'annihilation_gate_pass': True}
  - mode_index: 1, sigma: 4.7006989691328044e-11, classification: UNEXPLAINED_STIFFNESS, residual: 0.9876429467643332, lanes: {'psi_real': 8.029988181532186e-28, 'psi_imag': 1.4991198073132757e-06, 'a0': 1.8410544936760114e-30, 'ar': 0.007729895594373362, 'aw': 0.9922686052858193, 'r0': 1.80927381137277e-25, 'mu': 1.5568523920127209e-27}, support: {'best_generator': 'maxwell_residual_gauge', 'best_overlap': 0.11116225640259612, 'annihilation_gate_pass': False}
  - mode_index: 2, sigma: 2.6225709826789087e-10, classification: UNEXPLAINED_STIFFNESS, residual: 0.9999055431770426, lanes: {'psi_real': 1.3840986340073814e-28, 'psi_imag': 1.351818008692548e-05, 'a0': 9.206175320841355e-30, 'ar': 0.9562145732697735, 'aw': 0.04377190855013947, 'r0': 2.945885101621408e-26, 'mu': 1.380658128081834e-28}, support: {'best_generator': 'maxwell_residual_gauge', 'best_overlap': 0.009718889070657266, 'annihilation_gate_pass': False}
  - mode_index: 3, sigma: 5.926280559311367e-10, classification: UNEXPLAINED_STIFFNESS, residual: 0.6739393245410806, lanes: {'psi_real': 8.124175910251807e-29, 'psi_imag': 2.101845681068869e-05, 'a0': 1.0150268255194537e-30, 'ar': 0.9864546775079923, 'aw': 0.013524304035197277, 'r0': 1.430913637716034e-26, 'mu': 3.58600955696789e-28}, support: {'best_generator': 'maxwell_residual_gauge', 'best_overlap': 0.5710172285689894, 'annihilation_gate_pass': False}
  - mode_index: 4, sigma: 2.8112679473453323e-08, classification: UNEXPLAINED_STIFFNESS, residual: 0.40559198440068767, lanes: {'psi_real': 3.805748005127968e-29, 'psi_imag': 3.6259234598036724e-06, 'a0': 1.0495650841741824e-29, 'ar': 0.22236725410600497, 'aw': 0.7776291199705352, 'r0': 5.7545285884514175e-27, 'mu': 6.104839144856811e-29}, support: {'best_generator': 'maxwell_residual_gauge', 'best_overlap': 0.770978609041291, 'annihilation_gate_pass': False}
```

## Recommended Next Step

Next step: pin or deflate the explained symmetry mode(s), especially the global phase if present, then investigate the remaining field-lane residual subspace before a re-crawl.

## Guard Confirmation

```yaml
faithful_operator_boundary: {'physical_residual_entry_point': 'stage1_solver.coupled_branch.patha_closed_branch_residual', 'operators_modified_by_c0': False, 'frozen_physics_modified_by_c0': False, 'physical_export_guard_modified_by_c0': False}
scope_guard: {'diagnosis_only': True, 'gauge_fix_or_deflation_implemented': False, 'recrawl_implemented': False, 'faithful_operators_touched_by_c0c': False, 'frozen_physics_touched_by_c0c': False, 'physical_export_guard_touched_by_c0c': False}
```

Machine artifact: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0c_nullmode_identification/pathA_C0c_nullmode_identification.json`.
