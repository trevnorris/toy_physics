# Step 3 Coupled Branch Engineering Smoke

Overall engineering gate: PASS
Config hash: `dc3bb753f52df730`

**Discipline:** engineering smoke, placeholder parameters, not a physical packet, target-blind. No field-to-coefficient export is performed.

## Coupled Residual

Matter sector: stationary gauged quintic GNLS on `(r,w)`, `[-hbar^2/(2m) D_iD_i + V_conf + (5K/4)rho^4 + q_star A0 - mu] psi`; `D_i = d_i - i(q_star/hbar) A_i`, expanded with the axisymmetric full radial measure.
Maxwell sector: localized H=Z operator over `(A0, Ar, Aw)` with `d_M(Z F^{MN}) + xi^-1 d^N(Z d.A) = mu0 J_psi^N`; `J_psi^0=q_star rho`, `J_psi^i=q_star j^i`, and `j^i=(hbar/m) Im(psi* D_i psi)`.  No external matter source is added.
Sources: compact lines 556-583, 638-648, 651-659, 674-689; prereg lines 44-56; NONLINEAR_PROTOCOL_V2 lines 35-38.

## Placeholder Parameters

Label: engineering-smoke placeholders only; not a physical packet; target-blind

```yaml
hbar: 1.0
particle_mass: 1.0
gauge_charge: 0.35
mu0: 1.0
xi: 1.0
continuation_K_values: (0.05, 0.15, 0.3, 0.5)
mass: 1.0
localization_width: 0.75
localization_floor: 0.8
localization_amplitude: 0.45
r_mouth: 1.2
r_exit: 0.9
geometry_decay_length: 0.8
radial_wall_strength: 0.65
axial_trap_strength: 0.12
solve_grid: (10, 8)
ladder_levels: ((8, 6), (10, 8), (12, 10), (16, 12), (20, 14))
```

`R_0(w)=R_exit+(R_mouth-R_exit) exp(-(w-w_min)/geometry_decay_length)` is prescribed for this isotropic smoke; it satisfies the mouth value and finite open-exit Robin class.

## Coupled MMS

Full coupled stationary matter plus localized-Maxwell residual on the (r,w) grid: A_i enters D_i and the Maxwell vector blocks are sourced by q_star*j_psi.
Forcing: SymPy evaluated the continuum coupled operator directly: covariant axisymmetric D_iD_i for matter and H=Z localized Maxwell minus the matter charge/current source.
MMS placeholders:
```yaml
hbar: 1.0
particle_mass: 1.0
gauge_charge: 0.35
mu0: 1.0
xi: 1.0
eos_K: 0.45
chemical_potential: 0.9
localization_width: 0.75
localization_floor: 0.8
localization_amplitude: 0.45
r_mouth: 1.2
r_exit: 0.9
geometry_decay_length: 0.8
radial_wall_strength: 0.65
axial_trap_strength: 0.12
```

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

## Newton And Continuation

| eos_K | converged | iterations | initial_residual_norm | final_residual_norm | gmres_iterations | message |
| --- | --- | --- | --- | --- | --- | --- |
| 5.000000e-02 | True | 4 | 2.494727e+00 | 1.189295e-10 | [65, 68, 68, 58] | residual tolerance reached |
| 1.500000e-01 | True | 1 | 3.787584e-04 | 2.185050e-09 | [61] | residual tolerance reached |
| 3.000000e-01 | True | 1 | 5.677387e-04 | 5.175554e-09 | [62] | residual tolerance reached |
| 5.000000e-01 | True | 1 | 7.561868e-04 | 9.840448e-09 | [62] | residual tolerance reached |

Main solve `solve_nr_10_nw_8`: final residual linf=9.840448e-09, wall-clock=6.393921e+01s, peak RSS=6.164414e+02 MB, manifest=`software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/solve_nr_10_nw_8/manifest.json`.

## Jacobian Check

Coupled residual JVP vs centered finite difference: relative=9.113676e-12, absolute=7.628077e-10, epsilon=1.000000e-05, status=PASS.

## Resolution Ladder

| grid | dof | wall_clock_seconds | peak_memory_mb | newton_iterations | final_residual_linf | gmres_iterations | converged | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ladder_nr_8_nw_6 | 241 | 4.787101e+01 | 6.164414e+02 | 7 | 1.013693e-08 | [50, 52, 51, 44, 48, 48, 49] | True | continuation completed |
| ladder_nr_10_nw_8 | 401 | 4.513359e+01 | 6.164414e+02 | 5 | 9.840448e-09 | [71, 67, 61, 62, 62] | True | continuation completed |
| ladder_nr_12_nw_10 | 601 | 5.611007e+01 | 6.164414e+02 | 5 | 9.668407e-09 | [87, 80, 77, 79, 79] | True | continuation completed |
| ladder_nr_16_nw_12 | 961 | 8.112736e+01 | 6.164414e+02 | 5 | 9.517841e-09 | [128, 118, 102, 111, 115] | True | continuation completed |
| ladder_nr_20_nw_14 | 1401 | 1.408827e+02 | 6.164414e+02 | 5 | 9.430273e-09 | [184, 179, 190, 200, 198] | True | continuation completed |

Ladder stop reason: stopped after ladder_nr_20_nw_14 because the level exceeded 120.0s

## STOP And Flag

- isotropic wall-to-branch static balance: STOP_AND_FLAG. compact lines 6783-6785 keep the wall PDE as an effective closure unless S_Sigma[R] is promoted; NONLINEAR_PROTOCOL_V2 lines 14-19 state that the full coupled physical residual equations are not frozen.  For WP1, eta=0, so this smoke prescribes R_0(w) and does not invent a wall force law.

## Manifests

- coupled MMS nr_12_nw_10: `software/stage1_solver/runs/step3_coupled_branch_smoke/mms_coupled_branch_matter_maxwell/nr_12_nw_10/manifest.json`
- coupled MMS nr_24_nw_20: `software/stage1_solver/runs/step3_coupled_branch_smoke/mms_coupled_branch_matter_maxwell/nr_24_nw_20/manifest.json`
- coupled MMS nr_48_nw_40: `software/stage1_solver/runs/step3_coupled_branch_smoke/mms_coupled_branch_matter_maxwell/nr_48_nw_40/manifest.json`
- coupled MMS nr_96_nw_80: `software/stage1_solver/runs/step3_coupled_branch_smoke/mms_coupled_branch_matter_maxwell/nr_96_nw_80/manifest.json`
- main continuation: `software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/solve_nr_10_nw_8/manifest.json`
- ladder ladder_nr_8_nw_6: `software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/ladder_nr_8_nw_6/manifest.json`
- ladder ladder_nr_10_nw_8: `software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/ladder_nr_10_nw_8/manifest.json`
- ladder ladder_nr_12_nw_10: `software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/ladder_nr_12_nw_10/manifest.json`
- ladder ladder_nr_16_nw_12: `software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/ladder_nr_16_nw_12/manifest.json`
- ladder ladder_nr_20_nw_14: `software/stage1_solver/runs/step3_coupled_branch_smoke/step3_coupled_branch_engineering_smoke/ladder_nr_20_nw_14/manifest.json`

