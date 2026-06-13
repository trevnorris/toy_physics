# Step 4 Coupled Branch Grid-Convergence Study

Overall engineering gate: PASS
Config hash: `6ca33829da5aec8b`

**Discipline:** engineering smoke, placeholder parameters, not a physical packet, target-blind. Surrogates are raw-field functionals only; no extraction map is used.

## Ladder And Norm

Refinement ratio: `2x` in both grid directions.
Configured levels: `((6, 4), (12, 8), (24, 16), (48, 32), (96, 64))`.
Restriction: Nested finite-volume control volumes; fine fields are volume-averaged onto the next coarser grid before raw-field differences.
Difference norm: sqrt(sum_fields integral (u_h - R u_h/2)^2 dV) with dV=4*pi*r^2 dr dw; Linf is also reported.

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

## Level Performance

| level | grid | dof | spacing | wall_clock_seconds | peak_memory_mb | newton_iterations | final_residual_linf | gmres_max | gmres_mean | converged | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | convergence_l0_nr_6_nw_4 | 121 | 4.000000e-01 | 2.461210e+01 | 6.160117e+02 | 7 | 1.068892e-08 | 9 | 6.857143e+00 | true | continuation completed |
| 1 | convergence_l1_nr_12_nw_8 | 481 | 2.000000e-01 | 3.719478e+01 | 6.160117e+02 | 6 | 9.691578e-09 | 4 | 3.333333e+00 | true | continuation completed |
| 2 | convergence_l2_nr_24_nw_16 | 1921 | 1.000000e-01 | 3.664135e+01 | 6.344414e+02 | 5 | 9.375697e-09 | 3 | 2.600000e+00 | true | continuation completed |
| 3 | convergence_l3_nr_48_nw_32 | 7681 | 5.000000e-02 | 3.718557e+01 | 7.242852e+02 | 5 | 9.305613e-09 | 3 | 2.600000e+00 | true | continuation completed |
| 4 | convergence_l4_nr_96_nw_64 | 30721 | 2.500000e-02 | 4.137755e+01 | 1.093586e+03 | 5 | 9.284885e-09 | 3 | 2.200000e+00 | true | continuation completed |

## Raw-Field Self-Convergence

| coarse_grid | fine_grid | raw_field_l2_change | raw_field_relative_l2_change | raw_field_linf_change | observed_order_from_l2_change | observed_order_from_relative_l2_change |
| --- | --- | --- | --- | --- | --- | --- |
| convergence_l0_nr_6_nw_4 | convergence_l1_nr_12_nw_8 | 2.238856e-02 | 2.238173e-02 | 9.723233e-03 | - | - |
| convergence_l1_nr_12_nw_8 | convergence_l2_nr_24_nw_16 | 5.577530e-03 | 5.575966e-03 | 2.946680e-03 | 2.005063e+00 | 2.005028e+00 |
| convergence_l2_nr_24_nw_16 | convergence_l3_nr_48_nw_32 | 1.392046e-03 | 1.391664e-03 | 7.701851e-04 | 2.002419e+00 | 2.002411e+00 |
| convergence_l3_nr_48_nw_32 | convergence_l4_nr_96_nw_64 | 3.478452e-04 | 3.477501e-04 | 1.970138e-04 | 2.000690e+00 | 2.000688e+00 |

## Surrogate Observable Verdicts

| label | finest_grid | finest_value | last_observed_order | richardson_estimate | finest_error_estimate | verdict | diagnosis |
| --- | --- | --- | --- | --- | --- | --- | --- |
| density mass integral | convergence_l4_nr_96_nw_64 | 1.000000e+00 | 2.005109e+00 | 1.000000e+00 | 2.305045e-12 | solver-floor diagnostic | Newton mass constraint; differences read the solver floor, not discretization. |
| peak density | convergence_l4_nr_96_nw_64 | 2.800545e-01 | 2.005433e+00 | 2.801066e-01 | 5.205938e-05 | expected-order convergence | Pointwise extrema are sensitive to cell-center placement and coarse throat resolution. |
| minimum density | convergence_l4_nr_96_nw_64 | 2.309396e-07 | 2.020453e+00 | 2.452615e-09 | 2.284870e-07 | expected-order convergence | Pointwise extrema are sensitive to cell-center placement and coarse throat resolution. |
| raw field L2 norm | convergence_l4_nr_96_nw_64 | 1.000273e+00 | 2.008130e+00 | 1.000273e+00 | 1.142831e-07 | expected-order convergence | Raw-field integral on the shared coupled branch. |
| A0 L2 norm | convergence_l4_nr_96_nw_64 | 2.337608e-02 | 2.005847e+00 | 2.337118e-02 | 4.898945e-06 | expected-order convergence | Coupled gauge/current response can show reduced order from open Robin boundaries. |
| spatial gauge L2 norm | convergence_l4_nr_96_nw_64 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic | This raw channel is identically zero on the completed placeholder branch; no order is measured. |
| spatial current L2 norm | convergence_l4_nr_96_nw_64 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic | This raw channel is identically zero on the completed placeholder branch; no order is measured. |
| field-energy-like integral | convergence_l4_nr_96_nw_64 | 2.044752e+00 | 2.169728e+00 | 2.045116e+00 | 3.632239e-04 | expected-order convergence | Gradient-weighted integral; boundary and coarse-grid throat terms can reduce order. |
| chemical potential | convergence_l4_nr_96_nw_64 | 2.093577e+00 | 1.998526e+00 | 2.093825e+00 | 2.487751e-04 | expected-order convergence | Raw-field integral on the shared coupled branch. |
| final residual Linf | convergence_l4_nr_96_nw_64 | 9.284885e-09 | 1.757522e+00 | 9.276180e-09 | 8.704916e-12 | solver-floor diagnostic | Newton/GMRES stopping floor; useful as the numerical floor diagnostic. |

## Surrogate Per-Level Table

| label | grid | value | successive_change | observed_order | richardson_estimate | error_estimate | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| density mass integral | convergence_l0_nr_6_nw_4 | 1.000000e+00 | - | - | 1.000000e+00 | 6.901260e-10 | solver-floor diagnostic |
| density mass integral | convergence_l1_nr_12_nw_8 | 1.000000e+00 | 5.363407e-10 | - | 1.000000e+00 | 1.537852e-10 | solver-floor diagnostic |
| density mass integral | convergence_l2_nr_24_nw_16 | 1.000000e+00 | 1.166438e-10 | 2.201040e+00 | 1.000000e+00 | 3.714140e-11 | solver-floor diagnostic |
| density mass integral | convergence_l3_nr_48_nw_32 | 1.000000e+00 | 2.788880e-11 | 2.064352e+00 | 1.000000e+00 | 9.252599e-12 | solver-floor diagnostic |
| density mass integral | convergence_l4_nr_96_nw_64 | 1.000000e+00 | 6.947554e-12 | 2.005109e+00 | 1.000000e+00 | 2.305045e-12 | solver-floor diagnostic |
| peak density | convergence_l0_nr_6_nw_4 | 2.662493e-01 | - | - | 2.801066e-01 | 1.385730e-02 | expected-order convergence |
| peak density | convergence_l1_nr_12_nw_8 | 2.767194e-01 | 1.047014e-02 | - | 2.801066e-01 | 3.387155e-03 | expected-order convergence |
| peak density | convergence_l2_nr_24_nw_16 | 2.792673e-01 | 2.547908e-03 | 2.038896e+00 | 2.801066e-01 | 8.392470e-04 | expected-order convergence |
| peak density | convergence_l3_nr_48_nw_32 | 2.798976e-01 | 6.302238e-04 | 2.015377e+00 | 2.801066e-01 | 2.090232e-04 | expected-order convergence |
| peak density | convergence_l4_nr_96_nw_64 | 2.800545e-01 | 1.569638e-04 | 2.005433e+00 | 2.801066e-01 | 5.205938e-05 | expected-order convergence |
| minimum density | convergence_l0_nr_6_nw_4 | 6.560288e-05 | - | - | 2.452615e-09 | 6.560043e-05 | expected-order convergence |
| minimum density | convergence_l1_nr_12_nw_8 | 1.543808e-05 | 5.016480e-05 | - | 2.452615e-09 | 1.543563e-05 | expected-order convergence |
| minimum density | convergence_l2_nr_24_nw_16 | 3.763385e-06 | 1.167470e-05 | 2.103290e+00 | 2.452615e-09 | 3.760933e-06 | expected-order convergence |
| minimum density | convergence_l3_nr_48_nw_32 | 9.294501e-07 | 2.833935e-06 | 2.042507e+00 | 2.452615e-09 | 9.269974e-07 | expected-order convergence |
| minimum density | convergence_l4_nr_96_nw_64 | 2.309396e-07 | 6.985104e-07 | 2.020453e+00 | 2.452615e-09 | 2.284870e-07 | expected-order convergence |
| raw field L2 norm | convergence_l0_nr_6_nw_4 | 1.000305e+00 | - | - | 1.000273e+00 | 3.198571e-05 | expected-order convergence |
| raw field L2 norm | convergence_l1_nr_12_nw_8 | 1.000281e+00 | 2.445738e-05 | - | 1.000273e+00 | 7.528332e-06 | expected-order convergence |
| raw field L2 norm | convergence_l2_nr_24_nw_16 | 1.000275e+00 | 5.679076e-06 | 2.106542e+00 | 1.000273e+00 | 1.849256e-06 | expected-order convergence |
| raw field L2 norm | convergence_l3_nr_48_nw_32 | 1.000274e+00 | 1.389540e-06 | 2.031049e+00 | 1.000273e+00 | 4.597159e-07 | expected-order convergence |
| raw field L2 norm | convergence_l4_nr_96_nw_64 | 1.000273e+00 | 3.454328e-07 | 2.008130e+00 | 1.000273e+00 | 1.142831e-07 | expected-order convergence |
| A0 L2 norm | convergence_l0_nr_6_nw_4 | 2.470226e-02 | - | - | 2.337118e-02 | 1.331080e-03 | expected-order convergence |
| A0 L2 norm | convergence_l1_nr_12_nw_8 | 2.369120e-02 | 1.011057e-03 | - | 2.337118e-02 | 3.200234e-04 | expected-order convergence |
| A0 L2 norm | convergence_l2_nr_24_nw_16 | 2.345020e-02 | 2.410024e-04 | 2.068745e+00 | 2.337118e-02 | 7.902104e-05 | expected-order convergence |
| A0 L2 norm | convergence_l3_nr_48_nw_32 | 2.339085e-02 | 5.934568e-05 | 2.021832e+00 | 2.337118e-02 | 1.967536e-05 | expected-order convergence |
| A0 L2 norm | convergence_l4_nr_96_nw_64 | 2.337608e-02 | 1.477641e-05 | 2.005847e+00 | 2.337118e-02 | 4.898945e-06 | expected-order convergence |
| spatial gauge L2 norm | convergence_l0_nr_6_nw_4 | 0.000000e+00 | - | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial gauge L2 norm | convergence_l1_nr_12_nw_8 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial gauge L2 norm | convergence_l2_nr_24_nw_16 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial gauge L2 norm | convergence_l3_nr_48_nw_32 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial gauge L2 norm | convergence_l4_nr_96_nw_64 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial current L2 norm | convergence_l0_nr_6_nw_4 | 0.000000e+00 | - | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial current L2 norm | convergence_l1_nr_12_nw_8 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial current L2 norm | convergence_l2_nr_24_nw_16 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial current L2 norm | convergence_l3_nr_48_nw_32 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| spatial current L2 norm | convergence_l4_nr_96_nw_64 | 0.000000e+00 | 0.000000e+00 | - | 0.000000e+00 | 0.000000e+00 | null diagnostic |
| field-energy-like integral | convergence_l0_nr_6_nw_4 | 1.901939e+00 | - | - | 2.045116e+00 | 1.431763e-01 | expected-order convergence |
| field-energy-like integral | convergence_l1_nr_12_nw_8 | 2.011382e+00 | 1.094424e-01 | - | 2.045116e+00 | 3.373396e-02 | expected-order convergence |
| field-energy-like integral | convergence_l2_nr_24_nw_16 | 2.037762e+00 | 2.638069e-02 | 2.052618e+00 | 2.045116e+00 | 7.353273e-03 | expected-order convergence |
| field-energy-like integral | convergence_l3_nr_48_nw_32 | 2.043481e+00 | 5.718989e-03 | 2.205650e+00 | 2.045116e+00 | 1.634284e-03 | expected-order convergence |
| field-energy-like integral | convergence_l4_nr_96_nw_64 | 2.044752e+00 | 1.271060e-03 | 2.169728e+00 | 2.045116e+00 | 3.632239e-04 | expected-order convergence |
| chemical potential | convergence_l0_nr_6_nw_4 | 2.031140e+00 | - | - | 2.093825e+00 | 6.268502e-02 | expected-order convergence |
| chemical potential | convergence_l1_nr_12_nw_8 | 2.077986e+00 | 4.684599e-02 | - | 2.093825e+00 | 1.583903e-02 | expected-order convergence |
| chemical potential | convergence_l2_nr_24_nw_16 | 2.089853e+00 | 1.186676e-02 | 1.981000e+00 | 2.093825e+00 | 3.972273e-03 | expected-order convergence |
| chemical potential | convergence_l3_nr_48_nw_32 | 2.092831e+00 | 2.978190e-03 | 1.994418e+00 | 2.093825e+00 | 9.940838e-04 | expected-order convergence |
| chemical potential | convergence_l4_nr_96_nw_64 | 2.093577e+00 | 7.453087e-04 | 1.998526e+00 | 2.093825e+00 | 2.487751e-04 | expected-order convergence |
| final residual Linf | convergence_l0_nr_6_nw_4 | 1.068892e-08 | - | - | 9.276180e-09 | 1.412741e-09 | solver-floor diagnostic |
| final residual Linf | convergence_l1_nr_12_nw_8 | 9.691578e-09 | 9.973429e-10 | - | 9.276180e-09 | 4.153978e-10 | solver-floor diagnostic |
| final residual Linf | convergence_l2_nr_24_nw_16 | 9.375697e-09 | 3.158807e-10 | 1.658710e+00 | 9.276180e-09 | 9.951716e-11 | solver-floor diagnostic |
| final residual Linf | convergence_l3_nr_48_nw_32 | 9.305613e-09 | 7.008438e-11 | 2.172215e+00 | 9.276180e-09 | 2.943278e-11 | solver-floor diagnostic |
| final residual Linf | convergence_l4_nr_96_nw_64 | 9.284885e-09 | 2.072786e-11 | 1.757522e+00 | 9.276180e-09 | 8.704916e-12 | solver-floor diagnostic |

## Numerical Floor

Solver floor read: residual Linf <= 1.068892e-08, mass-constraint floor <= 2.084683e-09.
Finest raw-field relative self-difference: 3.477501e-04; successive-difference floor reached: false.
Preliminary numerical floor for scalar diagnostics: 1.068892e-08.

## Resolution Sizing

| label | last_observed_order | finest_error_estimate | finest_dof | dof_for_error_1e-03 | dof_for_error_1e-04 | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| density mass integral | 2.005109e+00 | 2.305045e-12 | 30721 | - | - | solver-floor diagnostic |
| peak density | 2.005433e+00 | 5.205938e-05 | 30721 | 30721 | 30721 | expected-order convergence |
| minimum density | 2.020453e+00 | 2.284870e-07 | 30721 | 30721 | 30721 | expected-order convergence |
| raw field L2 norm | 2.008130e+00 | 1.142831e-07 | 30721 | 30721 | 30721 | expected-order convergence |
| A0 L2 norm | 2.005847e+00 | 4.898945e-06 | 30721 | 30721 | 30721 | expected-order convergence |
| spatial gauge L2 norm | - | 0.000000e+00 | 30721 | - | - | null diagnostic |
| spatial current L2 norm | - | 0.000000e+00 | 30721 | - | - | null diagnostic |
| field-energy-like integral | 2.169728e+00 | 3.632239e-04 | 30721 | 30721 | 100877 | expected-order convergence |
| chemical potential | 1.998526e+00 | 2.487751e-04 | 30721 | 30721 | 76478 | expected-order convergence |
| final residual Linf | 1.757522e+00 | 8.704916e-12 | 30721 | - | - | solver-floor diagnostic |

Laptop limit read: max completed `convergence_l4_nr_96_nw_64` / 30721 DOF. Stop reason: completed configured ladder. Next 2x grid would be (192, 128) / 122881 DOF.
Direct sparse LU is the CPU-laptop limiter; the next 2x level is the engineering sizing point for a scalable preconditioner.

## Machine-Readable Output

Convergence table: `software/stage1_solver/runs/step4_convergence_study/step4_coupled_branch_grid_convergence/convergence_table.json`

Manifests:
- convergence_l0_nr_6_nw_4: `software/stage1_solver/runs/step4_convergence_study/step3_coupled_branch_engineering_smoke/convergence_l0_nr_6_nw_4/manifest.json`
- convergence_l1_nr_12_nw_8: `software/stage1_solver/runs/step4_convergence_study/step3_coupled_branch_engineering_smoke/convergence_l1_nr_12_nw_8/manifest.json`
- convergence_l2_nr_24_nw_16: `software/stage1_solver/runs/step4_convergence_study/step3_coupled_branch_engineering_smoke/convergence_l2_nr_24_nw_16/manifest.json`
- convergence_l3_nr_48_nw_32: `software/stage1_solver/runs/step4_convergence_study/step3_coupled_branch_engineering_smoke/convergence_l3_nr_48_nw_32/manifest.json`
- convergence_l4_nr_96_nw_64: `software/stage1_solver/runs/step4_convergence_study/step3_coupled_branch_engineering_smoke/convergence_l4_nr_96_nw_64/manifest.json`

