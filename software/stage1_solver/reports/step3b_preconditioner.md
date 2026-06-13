# Step 3b Coupled Branch JFNK Preconditioner

Overall engineering gate: PASS
Preconditioned config hash: `537b01b3b94ab7c2`
Baseline config hash: `5fc7c70c59e94660`

**Discipline:** engineering smoke, placeholder parameters, not a physical packet, target-blind. No field-to-coefficient export is performed.

## Method

Interface: `solve_newton_jvp` accepts a left-preconditioner factory and passes the resulting SciPy `LinearOperator` as GMRES `M`.  The residual and JVP used by GMRES are unchanged.
Concrete preconditioner: assembled sparse Jacobian from the coupled residual JVP, colored by the radius-3 local stencil, factored with SciPy SuperLU, then used as an inverse preconditioner.  No new dependency was added.

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

## Correctness Preservation

Fixed-grid solution comparison `comparison_nr_10_nw_8` / 401 DOF: linf difference=1.242532e-12, l2 difference=6.796410e-12, tolerance=1.000000e-06, status=PASS.
Coupled residual JVP vs centered finite difference on the preconditioned solve: relative=9.453008e-12, absolute=7.912095e-10, epsilon=1.000000e-05, status=PASS.

| grid | spacing | error | observed_order | reference_norm | spatial_gauge_l2 | spatial_current_l2 |
| --- | --- | --- | --- | --- | --- | --- |
| nr_12_nw_10 | 1.666667e-01 | 7.800679e-01 | - | 3.725578e+01 | 2.123486e-01 | 3.480506e-01 |
| nr_24_nw_20 | 8.333333e-02 | 2.142920e-01 | 1.864022e+00 | 3.770192e+01 | 2.118832e-01 | 3.594497e-01 |
| nr_48_nw_40 | 4.166667e-02 | 5.519830e-02 | 1.956882e+00 | 3.781388e+01 | 2.117665e-01 | 3.624551e-01 |
| nr_96_nw_80 | 2.083333e-02 | 1.384584e-02 | 1.995171e+00 | 3.784190e+01 | 2.117373e-01 | 3.632263e-01 |

MMS checks:
- observed_order: PASS
- final_error: PASS
- cross_sector_gauge_nonzero: PASS
- cross_sector_current_nonzero: PASS

## Before After GMRES Curve

| grid | dof | baseline_gmres_max | preconditioned_gmres_max | baseline_gmres_mean | preconditioned_gmres_mean | baseline_seconds | preconditioned_seconds |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ladder_nr_8_nw_6 | 241 | 52 | 9 | 4.885714e+01 | 7.142857e+00 | 5.986401e+01 | 3.693172e+01 |
| ladder_nr_10_nw_8 | 401 | 71 | 3 | 6.460000e+01 | 2.600000e+00 | 4.442971e+01 | 3.607925e+01 |
| ladder_nr_12_nw_10 | 601 | 87 | 3 | 8.040000e+01 | 2.600000e+00 | 5.486192e+01 | 3.619512e+01 |
| ladder_nr_16_nw_12 | 961 | 128 | 3 | 1.148000e+02 | 2.600000e+00 | 7.976845e+01 | 3.585807e+01 |
| ladder_nr_20_nw_14 | 1401 | 200 | 3 | 1.902000e+02 | 2.600000e+00 | 1.333739e+02 | 3.626942e+01 |

Baseline max-GMRES growth on the shared ladder: 3.846154e+00. Preconditioned growth on the same ladder: 3.333333e-01.

## Extended Preconditioned Ladder

| grid | dof | wall_clock_seconds | peak_memory_mb | newton_iterations | final_residual_linf | gmres_iterations | gmres_max | gmres_mean | converged | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ladder_nr_8_nw_6 | 241 | 3.693172e+01 | 6.165078e+02 | 7 | 1.013057e-08 | [1, 9, 9, 6, 9, 8, 8] | 9 | 7.142857e+00 | True | continuation completed |
| ladder_nr_10_nw_8 | 401 | 3.607925e+01 | 6.165078e+02 | 5 | 9.801077e-09 | [1, 3, 3, 3, 3] | 3 | 2.600000e+00 | True | continuation completed |
| ladder_nr_12_nw_10 | 601 | 3.619512e+01 | 6.165078e+02 | 5 | 9.669862e-09 | [1, 3, 3, 3, 3] | 3 | 2.600000e+00 | True | continuation completed |
| ladder_nr_16_nw_12 | 961 | 3.585807e+01 | 6.213750e+02 | 5 | 9.507307e-09 | [1, 3, 3, 3, 3] | 3 | 2.600000e+00 | True | continuation completed |
| ladder_nr_20_nw_14 | 1401 | 3.626942e+01 | 6.271914e+02 | 5 | 9.430104e-09 | [1, 3, 3, 3, 3] | 3 | 2.600000e+00 | True | continuation completed |
| ladder_nr_24_nw_16 | 1921 | 3.610877e+01 | 6.348359e+02 | 5 | 9.385590e-09 | [1, 2, 3, 3, 3] | 3 | 2.400000e+00 | True | continuation completed |
| ladder_nr_28_nw_20 | 2801 | 3.628103e+01 | 6.481445e+02 | 5 | 9.358243e-09 | [1, 2, 3, 3, 3] | 3 | 2.400000e+00 | True | continuation completed |
| ladder_nr_32_nw_24 | 3841 | 3.609601e+01 | 6.638828e+02 | 5 | 9.340123e-09 | [1, 2, 3, 3, 3] | 3 | 2.400000e+00 | True | continuation completed |
| ladder_nr_40_nw_28 | 5601 | 3.673894e+01 | 6.908789e+02 | 5 | 9.318366e-09 | [1, 2, 2, 3, 3] | 3 | 2.200000e+00 | True | continuation completed |

New maximum converged DOF on this laptop run: 5601. Ladder stop reason: completed configured ladder.

## Main Preconditioned Solve

`solve_nr_10_nw_8` final residual linf=9.843399e-09, wall-clock=1.131815e+02s, peak RSS=6.165078e+02 MB, manifest=`software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/solve_nr_10_nw_8/manifest.json`.

## D2 Trigger Read

The Krylov counts flatten when the linearized coupled operator is preconditioned by an assembled sparse Jacobian inverse.  This supports the D2 premise that the physics residual can remain in the torch-owned operator stack, but direct sparse LU is a CPU smoke method, not the production preconditioner.  The production path should replace this inverse with structured multigrid or a PETSc-class algebraic equivalent.

## Manifests

- coupled MMS nr_12_nw_10: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/mms_coupled_branch_matter_maxwell/nr_12_nw_10/manifest.json`
- coupled MMS nr_24_nw_20: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/mms_coupled_branch_matter_maxwell/nr_24_nw_20/manifest.json`
- coupled MMS nr_48_nw_40: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/mms_coupled_branch_matter_maxwell/nr_48_nw_40/manifest.json`
- coupled MMS nr_96_nw_80: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/mms_coupled_branch_matter_maxwell/nr_96_nw_80/manifest.json`
- comparison unpreconditioned: `software/stage1_solver/runs/step3b_preconditioner/unpreconditioned/step3_coupled_branch_engineering_smoke/comparison_nr_10_nw_8_unpreconditioned/manifest.json`
- comparison preconditioned: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/comparison_nr_10_nw_8_preconditioned/manifest.json`
- main preconditioned solve: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/solve_nr_10_nw_8/manifest.json`
- baseline ladder ladder_nr_8_nw_6: `software/stage1_solver/runs/step3b_preconditioner/unpreconditioned/step3_coupled_branch_engineering_smoke/ladder_nr_8_nw_6/manifest.json`
- baseline ladder ladder_nr_10_nw_8: `software/stage1_solver/runs/step3b_preconditioner/unpreconditioned/step3_coupled_branch_engineering_smoke/ladder_nr_10_nw_8/manifest.json`
- baseline ladder ladder_nr_12_nw_10: `software/stage1_solver/runs/step3b_preconditioner/unpreconditioned/step3_coupled_branch_engineering_smoke/ladder_nr_12_nw_10/manifest.json`
- baseline ladder ladder_nr_16_nw_12: `software/stage1_solver/runs/step3b_preconditioner/unpreconditioned/step3_coupled_branch_engineering_smoke/ladder_nr_16_nw_12/manifest.json`
- baseline ladder ladder_nr_20_nw_14: `software/stage1_solver/runs/step3b_preconditioner/unpreconditioned/step3_coupled_branch_engineering_smoke/ladder_nr_20_nw_14/manifest.json`
- preconditioned ladder ladder_nr_8_nw_6: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_8_nw_6/manifest.json`
- preconditioned ladder ladder_nr_10_nw_8: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_10_nw_8/manifest.json`
- preconditioned ladder ladder_nr_12_nw_10: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_12_nw_10/manifest.json`
- preconditioned ladder ladder_nr_16_nw_12: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_16_nw_12/manifest.json`
- preconditioned ladder ladder_nr_20_nw_14: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_20_nw_14/manifest.json`
- preconditioned ladder ladder_nr_24_nw_16: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_24_nw_16/manifest.json`
- preconditioned ladder ladder_nr_28_nw_20: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_28_nw_20/manifest.json`
- preconditioned ladder ladder_nr_32_nw_24: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_32_nw_24/manifest.json`
- preconditioned ladder ladder_nr_40_nw_28: `software/stage1_solver/runs/step3b_preconditioner/colored_sparse_jacobian_lu/step3_coupled_branch_engineering_smoke/ladder_nr_40_nw_28/manifest.json`

