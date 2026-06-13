# Step 5 Boundary-Control Characterization

Overall engineering gate: PASS
Config hash: `fd9bc1b134624e96`

**Discipline:** engineering smoke, placeholder parameters, stationary branch, target-blind. This report uses raw-field interior diagnostics only.

## Method

Truncation sweep radial extents: `[2.8, 3.2, 3.6]`.
Impedance coefficient scales: `[0.5, 1.0, 2.0]` on fixed `r_max=3.2`.
Fixed resolution: uniform cell spacing 0.05; truncation changes add or remove outer radial cells only.
The acceptance sweep uses radial placements at or beyond 2.8. Nearer stress placements down to 2.4 showed residual boundary-zone contamination and are not treated as adequate placements.
Interior window: `{'r_min': 0.0, 'r_max': 1.2, 'w_min': 0.2, 'w_max': 1.0}`.
Difference norm: sqrt(sum_fields integral_interior (u - u_ref)^2 dV) with dV=4*pi*r^2 dr dw, divided by the reference interior raw-field L2.
Reference scale: interior raw-field L2 magnitude; threshold is a fixed multiple of the step-4 raw-field relative self-difference.
Threshold: 1.391000e-03 relative to the interior signal (4.0 x the step-4 reference floor 3.477501e-04).
Pass criterion: the interior raw-field L2 relative change is the governing boundary-characterization pass criterion; per-surrogate-observable relative changes are diagnostic/informational only and are not thresholded.

Preconditioner:
```yaml
type: colored_sparse_jacobian_lu
side: left
rebuild_policy: once_per_continuation
stencil_radius: 3
color_separation: 7
factorization: splu
diagonal_shift: 0.0
drop_tolerance: 0.0
fill_factor: 10.0
permutation: COLAMD
```

## Solve Rows

| sweep | setting | variable_value | nr | nw | r_max | w_max | dof | wall_clock_seconds | final_residual_linf | converged | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| robin_only_truncation_placement | robin_only_rmax_2.8_nr_56_nw_32 | 2.800000e+00 | 56 | 32 | 2.800000e+00 | 1.600000e+00 | 8961 | 4.628495e+01 | 6.654378e-09 | true | continuation completed |
| robin_only_truncation_placement | robin_only_rmax_3.2_nr_64_nw_32 | 3.200000e+00 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 4.745092e+01 | 6.653254e-09 | true | continuation completed |
| robin_only_truncation_placement | robin_only_rmax_3.6_nr_72_nw_32 | 3.600000e+00 | 72 | 32 | 3.600000e+00 | 1.600000e+00 | 11521 | 4.860831e+01 | 6.653174e-09 | true | continuation completed |
| robin_only_impedance_coefficients | robin_only_alpha_scale_0.5_nr_64_nw_32 | 5.000000e-01 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 4.806175e+01 | 4.926383e-09 | true | continuation completed |
| robin_only_impedance_coefficients | robin_only_alpha_scale_1.0_nr_64_nw_32 | 1.000000e+00 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 4.718961e+01 | 6.653254e-09 | true | continuation completed |
| robin_only_impedance_coefficients | robin_only_alpha_scale_2.0_nr_64_nw_32 | 2.000000e+00 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 4.774126e+01 | 1.105997e-08 | true | continuation completed |
| sponge_truncation_placement | sponge_rmax_2.8_nr_56_nw_32 | 2.800000e+00 | 56 | 32 | 2.800000e+00 | 1.600000e+00 | 8961 | 8.811567e+01 | 1.368292e-11 | true | continuation completed |
| sponge_truncation_placement | sponge_rmax_3.2_nr_64_nw_32 | 3.200000e+00 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 1.038200e+02 | 1.387473e-11 | true | continuation completed |
| sponge_truncation_placement | sponge_rmax_3.6_nr_72_nw_32 | 3.600000e+00 | 72 | 32 | 3.600000e+00 | 1.600000e+00 | 11521 | 1.630251e+02 | 1.388620e-11 | true | continuation completed |
| sponge_impedance_coefficients | sponge_alpha_scale_0.5_nr_64_nw_32 | 5.000000e-01 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 7.950143e+01 | 1.387607e-11 | true | continuation completed |
| sponge_impedance_coefficients | sponge_alpha_scale_1.0_nr_64_nw_32 | 1.000000e+00 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 5.265616e+01 | 1.387473e-11 | true | continuation completed |
| sponge_impedance_coefficients | sponge_alpha_scale_2.0_nr_64_nw_32 | 2.000000e+00 | 64 | 32 | 3.200000e+00 | 1.600000e+00 | 10241 | 5.229882e+01 | 1.387176e-11 | true | continuation completed |

## Robin-Only Diagnostic

Robin-only exceeded the criterion, so the same sweeps were rerun with the smooth sponge layer enabled. The failed diagnostic is retained here.

Robin-only truncation:
| setting | reference_setting | variable_value | interior_relative_l2_change | relative_to_step4_reference_floor | max_interior_surrogate_relative_change | interior_linf_change | converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| robin_only_rmax_2.8_nr_56_nw_32 | robin_only_rmax_3.6_nr_72_nw_32 | 2.800000e+00 | 1.463584e-04 | 4.208724e-01 | 7.531948e-03 | 6.113086e-05 | true |
| robin_only_rmax_3.2_nr_64_nw_32 | robin_only_rmax_3.6_nr_72_nw_32 | 3.200000e+00 | 3.920517e-05 | 1.127395e-01 | 2.017254e-03 | 1.639526e-05 | true |
| robin_only_rmax_3.6_nr_72_nw_32 | robin_only_rmax_3.6_nr_72_nw_32 | 3.600000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |

Robin-only impedance:
| setting | reference_setting | variable_value | interior_relative_l2_change | relative_to_step4_reference_floor | max_interior_surrogate_relative_change | interior_linf_change | converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| robin_only_alpha_scale_0.5_nr_64_nw_32 | robin_only_alpha_scale_1.0_nr_64_nw_32 | 5.000000e-01 | 1.233220e-02 | 3.546283e+01 | 2.595529e-01 | 9.303247e-03 | true |
| robin_only_alpha_scale_1.0_nr_64_nw_32 | robin_only_alpha_scale_1.0_nr_64_nw_32 | 1.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |
| robin_only_alpha_scale_2.0_nr_64_nw_32 | robin_only_alpha_scale_1.0_nr_64_nw_32 | 2.000000e+00 | 2.019158e-02 | 5.806349e+01 | 6.623398e-01 | 1.533832e-02 | true |

Robin-only achieved max interior relative L2 boundary error: 2.019158e-02; threshold ratio 1.451587e+01.
Physical note: the Robin-only impedance sensitivity is driven by the long-range A0 scalar-potential tail: the matter field psi is tightly confined and decays fast, but the A0 Poisson-type potential has a slowly-decaying tail that reaches the exit, so at fixed truncation the interior remains sensitive to how the Robin impedance terminates that tail, whereas moving the truncation outward only lowers the already-small tail amplitude.

## Truncation-Placement Sweep

| setting | reference_setting | variable_value | interior_relative_l2_change | relative_to_step4_reference_floor | max_interior_surrogate_relative_change | interior_linf_change | converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sponge_rmax_2.8_nr_56_nw_32 | sponge_rmax_3.6_nr_72_nw_32 | 2.800000e+00 | 6.850290e-05 | 1.969889e-01 | 1.463913e-02 | 6.172193e-05 | true |
| sponge_rmax_3.2_nr_64_nw_32 | sponge_rmax_3.6_nr_72_nw_32 | 3.200000e+00 | 2.428884e-06 | 6.984566e-03 | 8.254997e-04 | 1.071780e-06 | true |
| sponge_rmax_3.6_nr_72_nw_32 | sponge_rmax_3.6_nr_72_nw_32 | 3.600000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |

Max non-reference truncation relative L2 change: 6.850290e-05.

## Impedance-Coefficient Sweep

| setting | reference_setting | variable_value | interior_relative_l2_change | relative_to_step4_reference_floor | max_interior_surrogate_relative_change | interior_linf_change | converged |
| --- | --- | --- | --- | --- | --- | --- | --- |
| sponge_alpha_scale_0.5_nr_64_nw_32 | sponge_alpha_scale_1.0_nr_64_nw_32 | 5.000000e-01 | 8.066541e-06 | 2.319637e-02 | 1.084896e-04 | 8.558949e-06 | true |
| sponge_alpha_scale_1.0_nr_64_nw_32 | sponge_alpha_scale_1.0_nr_64_nw_32 | 1.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 | true |
| sponge_alpha_scale_2.0_nr_64_nw_32 | sponge_alpha_scale_1.0_nr_64_nw_32 | 2.000000e+00 | 1.561899e-05 | 4.491440e-02 | 2.140471e-04 | 1.659750e-05 | true |

Max non-reference impedance relative L2 change: 1.561899e-05.

## Boundary-Error Metric

Achieved max interior relative L2 boundary error: 6.850290e-05; threshold ratio 4.924721e-02; step-4-floor ratio 1.969889e-01.

## Verdict

Robin-only boundary treatment exceeded the stationary criterion; the smooth sponge layer is required and brought the measured interior boundary error below threshold.
Sponge/absorber used: true; parameters: {'width': 0.4, 'matter_strength': 100.0, 'gauge_strength': 100.0, 'power': 2}.

## Deferral Note

The outgoing-wave reflection coefficient is deferred to the step-8 linearized tangent because the stationary isotropic branch has no propagating waves to reflect.

## Machine-Readable Output

Boundary characterization table: `software/stage1_solver/runs/step5_boundary_characterization/step5_boundary_characterization/boundary_characterization_table.json`

Manifests:
- robin_only_rmax_2.8_nr_56_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/robin_only_truncation_placement/step5_boundary_characterization_engineering_smoke/robin_only_rmax_2.8_nr_56_nw_32/manifest.json`
- robin_only_rmax_3.2_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/robin_only_truncation_placement/step5_boundary_characterization_engineering_smoke/robin_only_rmax_3.2_nr_64_nw_32/manifest.json`
- robin_only_rmax_3.6_nr_72_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/robin_only_truncation_placement/step5_boundary_characterization_engineering_smoke/robin_only_rmax_3.6_nr_72_nw_32/manifest.json`
- robin_only_alpha_scale_0.5_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/robin_only_impedance_coefficients/step5_boundary_characterization_engineering_smoke/robin_only_alpha_scale_0.5_nr_64_nw_32/manifest.json`
- robin_only_alpha_scale_1.0_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/robin_only_impedance_coefficients/step5_boundary_characterization_engineering_smoke/robin_only_alpha_scale_1.0_nr_64_nw_32/manifest.json`
- robin_only_alpha_scale_2.0_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/robin_only_impedance_coefficients/step5_boundary_characterization_engineering_smoke/robin_only_alpha_scale_2.0_nr_64_nw_32/manifest.json`
- sponge_rmax_2.8_nr_56_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/sponge_truncation_placement/step5_boundary_characterization_engineering_smoke/sponge_rmax_2.8_nr_56_nw_32/manifest.json`
- sponge_rmax_3.2_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/sponge_truncation_placement/step5_boundary_characterization_engineering_smoke/sponge_rmax_3.2_nr_64_nw_32/manifest.json`
- sponge_rmax_3.6_nr_72_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/sponge_truncation_placement/step5_boundary_characterization_engineering_smoke/sponge_rmax_3.6_nr_72_nw_32/manifest.json`
- sponge_alpha_scale_0.5_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/sponge_impedance_coefficients/step5_boundary_characterization_engineering_smoke/sponge_alpha_scale_0.5_nr_64_nw_32/manifest.json`
- sponge_alpha_scale_1.0_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/sponge_impedance_coefficients/step5_boundary_characterization_engineering_smoke/sponge_alpha_scale_1.0_nr_64_nw_32/manifest.json`
- sponge_alpha_scale_2.0_nr_64_nw_32: `software/stage1_solver/runs/step5_boundary_characterization/sponge_impedance_coefficients/step5_boundary_characterization_engineering_smoke/sponge_alpha_scale_2.0_nr_64_nw_32/manifest.json`

