# Step 2 Manufactured-Solution Validation

Overall gate: PASS
Config hash: `4b12c4f897dbe4ec`

## Config Hash Hygiene

`HarnessConfig.config_hash()` excludes `run_root` and `report_path`; manifests still record the full config.

## Step 1 Regression

Linear regression: PASS.
Cubic GPE regression: PASS.
Cubic JVP check: relative=8.542902e-13, absolute=1.143120e-09.

## MMS Solver Controls

```yaml
dtype: float64
device: cpu
seed: 20260612
residual_atol: 1e-11
residual_rtol: 1e-10
step_atol: 1e-13
step_rtol: 1e-12
max_newton_iters: 24
residual_norm: linf
linear_solver: gmres_jvp
gmres_rtol: 1e-11
gmres_atol: 1e-13
gmres_restart: 256
gmres_maxiter: 8
line_search: armijo
line_search_c1: 0.0001
line_search_shrink: 0.5
max_line_search_iters: 16
accept_best_line_search_decrease: True
finite_difference_jvp_epsilon: 0.0001
```

## Operator Forms

Matter: `-(hbar^2/2m) r^-2 d_r(r^2 d_r psi) + V psi + (5K/4)|psi|^8 psi - mu psi`.
Tensor Laplacian: cell FV divergence with radial flux `4*pi*r_face^2*dw*d_r u` and w flux `shell_volume*d_w u`.
Current: `(hbar/m) Im(conj(psi) d_r psi)` with `A_r=0`.
Maxwell: stationary axisymmetric `d_M(Z F^{MN}) + xi^-1 d^N(Z d.A)` with components `A0, Ar, Aw` and `H=Z`.
Wall: flat-`dw` stationary modal `-d_w(T_w d_w eta) + [K_eta + ell(ell+1)T_Omega] eta`.

## mms_quintic_matter_radial

Radial stationary gauged-GNLS operator with A=0 and h(rho)=5K*rho^4/4.
Source: compact lines 556-583 and 638-648.
Manufactured field: `11*(1 - r**2/R**2)**4/50 + 1 + 3*r**2*(1 - r**2/R**2)**4/(100*R**2)`
Forcing: SymPy applied the continuum radial operator to the manufactured field.

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_32 | 6.250000e-02 | 1.864662e-03 | - | 1.902143e+00 |
| nr_64 | 3.125000e-02 | 4.690801e-04 | 1.991008e+00 | 1.900330e+00 |
| nr_128 | 1.562500e-02 | 1.174552e-04 | 1.997724e+00 | 1.899876e+00 |
| nr_256 | 7.812500e-03 | 2.937544e-05 | 1.999427e+00 | 1.899762e+00 |

Finest-grid error: 2.937544e-05 at `nr_256`.

Matter JVP check: relative=1.022082e-11, absolute=1.433641e-07, epsilon=1.000000e-04.

Checks:
- observed_order: PASS
- final_error: PASS
- jacobian_relative: PASS
- jacobian_absolute: PASS

## mms_tensor_laplacian_2d

Genuine 2D tensor-grid scalar Laplacian with nontrivial radial and w dependence and all face boundary operators applied.
Source: Step-1 FV tensor operator carrying compact/prereg full radial measure convention.
Manufactured field: `(1 - r**2/R**2)**4/5 + 192*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(5*(w_max - w_min)**4) + r**2*(1 - r**2/R**2)**4/(25*R**2)`
Forcing: SymPy evaluated r^-2*d_r(r^2*d_r u)+d_w^2 u.

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_24_nw_20 | 8.333333e-02 | 1.760197e-01 | - | 7.779101e+00 |
| nr_48_nw_40 | 4.166667e-02 | 4.537784e-02 | 1.955677e+00 | 7.778874e+00 |
| nr_96_nw_80 | 2.083333e-02 | 1.143358e-02 | 1.988711e+00 | 7.778777e+00 |
| nr_192_nw_160 | 1.041667e-02 | 2.864070e-03 | 1.997139e+00 | 7.778752e+00 |
| nr_384_nw_320 | 5.208333e-03 | 7.163755e-04 | 1.999279e+00 | 7.778746e+00 |
| nr_768_nw_640 | 2.604167e-03 | 1.791163e-04 | 1.999819e+00 | 7.778745e+00 |

Finest-grid error: 1.791163e-04 at `nr_768_nw_640`.

Checks:
- observed_order: PASS
- final_error: PASS

## mms_complex_current_radial

Complex radial manufactured field with nonzero A=0 current.
Source: compact lines 651-659.
Manufactured field: `amplitude=r**2/20 + 3*cos(pi*r**2/R**2)/25 + 1; phase=2*r/5 + 17*sin(17*r/10)/100`
Forcing: SymPy evaluated (hbar/m)*amplitude^2*d_r phase.

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_32 | 6.250000e-02 | 3.500573e-03 | - | 1.290265e+00 |
| nr_64 | 3.125000e-02 | 8.150539e-04 | 2.102624e+00 | 1.289962e+00 |
| nr_128 | 1.562500e-02 | 1.936489e-04 | 2.073452e+00 | 1.289887e+00 |
| nr_256 | 7.812500e-03 | 4.692890e-05 | 2.044895e+00 | 1.289868e+00 |
| nr_512 | 3.906250e-03 | 1.153071e-05 | 2.024995e+00 | 1.289863e+00 |

Finest-grid error: 1.153071e-05 at `nr_512`.

Checks:
- observed_order: PASS
- final_error: PASS
- nonzero_current: PASS

## mms_localized_maxwell_h_equals_z

Stationary axisymmetric localized Maxwell operator on (r,w), components A0/Ar/Aw, H=Z.
Source: compact lines 590-630 and 674-689; prereg section D gauge and mixed-channel rows.
Manufactured field: `A0=(1 - r**2/R**2)**4/5 + 192*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(5*(w_max - w_min)**4); Ar=r*(1 - r**2/R**2)**4*(448*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(25*(w_max - w_min)**4) + 19/100); Aw=(1 - r**2/R**2)**4*(576*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(25*(w_max - w_min)**4) + 11/100) + 384*r**2*(1 - r**2/R**2)**4*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(25*R**2*(w_max - w_min)**4)`
Forcing: SymPy evaluated the displayed H=Z Maxwell equation after stationary axisymmetric reduction.

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nr_24_nw_20 | 8.333333e-02 | 1.836694e-01 | - | 7.520117e+00 |
| nr_48_nw_40 | 4.166667e-02 | 4.621307e-02 | 1.990738e+00 | 7.519158e+00 |
| nr_96_nw_80 | 2.083333e-02 | 1.145813e-02 | 2.011929e+00 | 7.518914e+00 |
| nr_192_nw_160 | 1.041667e-02 | 2.839842e-03 | 2.012489e+00 | 7.518853e+00 |
| nr_384_nw_320 | 5.208333e-03 | 7.058768e-04 | 2.008322e+00 | 7.518837e+00 |

Finest-grid error: 7.058768e-04 at `nr_384_nw_320`.

Checks:
- observed_order: PASS
- final_error: PASS

## mms_wall_s_eta_2

Stationary modal wall operator from S_eta^(2) in densitized flat-dw convention.
Source: research/pde_ledger/notes/stages/moving_throat_pde_stage001_geometry_lift.md lines 198-225.
Manufactured field: `256*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(5*(w_max - w_min)**4) + 1`
Forcing: SymPy varied the stationary modal action to -d_w(T_w d_w eta)+[K_eta+ell(ell+1)T_Omega]eta.
Previous STOP resolved: MMS certifies discretization of the pinned form; these placeholders are not physical constitutive values.
Placeholder caveat: MMS-only structural-certification placeholders; NOT the physical wall packet, which is `free_choice` (prereg §E) and is frozen per-run at solve time.
Placeholder coefficients:
- mu_eta: `1.0`
- T_w(w): `T_w_amp*sin(2*pi*(w - w_min)/(w_max - w_min)) + T_w_base`
- T_Omega(w): `T_Omega_amp*cos(2*pi*(w - w_min)/(w_max - w_min)) + T_Omega_base`
- K_eta(w): `256*K_eta_amp*(w - w_min)**4*(-(w - w_min)/(w_max - w_min) + 1)**4/(w_max - w_min)**4 + K_eta_base`

| grid | spacing | error | observed_order | reference_norm |
| --- | --- | --- | --- | --- |
| nw_32 | 5.000000e-02 | 1.829135e-02 | - | 8.101454e+00 |
| nw_64 | 2.500000e-02 | 4.627190e-03 | 1.982953e+00 | 8.101483e+00 |
| nw_128 | 1.250000e-02 | 1.160276e-03 | 1.995668e+00 | 8.101485e+00 |
| nw_256 | 6.250000e-03 | 2.902891e-04 | 1.998906e+00 | 8.101485e+00 |
| nw_512 | 3.125000e-03 | 7.258611e-05 | 1.999725e+00 | 8.101485e+00 |

Finest-grid error: 7.258611e-05 at `nw_512`.

Checks:
- observed_order: PASS
- final_error: PASS

## Manifests

- mms_quintic_matter_radial nr_32: `software/stage1_solver/runs/step2_manufactured_solutions/mms_quintic_matter_radial/nr_32/manifest.json`
- mms_quintic_matter_radial nr_64: `software/stage1_solver/runs/step2_manufactured_solutions/mms_quintic_matter_radial/nr_64/manifest.json`
- mms_quintic_matter_radial nr_128: `software/stage1_solver/runs/step2_manufactured_solutions/mms_quintic_matter_radial/nr_128/manifest.json`
- mms_quintic_matter_radial nr_256: `software/stage1_solver/runs/step2_manufactured_solutions/mms_quintic_matter_radial/nr_256/manifest.json`
- mms_tensor_laplacian_2d nr_24_nw_20: `software/stage1_solver/runs/step2_manufactured_solutions/mms_tensor_laplacian_2d/nr_24_nw_20/manifest.json`
- mms_tensor_laplacian_2d nr_48_nw_40: `software/stage1_solver/runs/step2_manufactured_solutions/mms_tensor_laplacian_2d/nr_48_nw_40/manifest.json`
- mms_tensor_laplacian_2d nr_96_nw_80: `software/stage1_solver/runs/step2_manufactured_solutions/mms_tensor_laplacian_2d/nr_96_nw_80/manifest.json`
- mms_tensor_laplacian_2d nr_192_nw_160: `software/stage1_solver/runs/step2_manufactured_solutions/mms_tensor_laplacian_2d/nr_192_nw_160/manifest.json`
- mms_tensor_laplacian_2d nr_384_nw_320: `software/stage1_solver/runs/step2_manufactured_solutions/mms_tensor_laplacian_2d/nr_384_nw_320/manifest.json`
- mms_tensor_laplacian_2d nr_768_nw_640: `software/stage1_solver/runs/step2_manufactured_solutions/mms_tensor_laplacian_2d/nr_768_nw_640/manifest.json`
- mms_complex_current_radial nr_32: `software/stage1_solver/runs/step2_manufactured_solutions/mms_complex_current_radial/nr_32/manifest.json`
- mms_complex_current_radial nr_64: `software/stage1_solver/runs/step2_manufactured_solutions/mms_complex_current_radial/nr_64/manifest.json`
- mms_complex_current_radial nr_128: `software/stage1_solver/runs/step2_manufactured_solutions/mms_complex_current_radial/nr_128/manifest.json`
- mms_complex_current_radial nr_256: `software/stage1_solver/runs/step2_manufactured_solutions/mms_complex_current_radial/nr_256/manifest.json`
- mms_complex_current_radial nr_512: `software/stage1_solver/runs/step2_manufactured_solutions/mms_complex_current_radial/nr_512/manifest.json`
- mms_localized_maxwell_h_equals_z nr_24_nw_20: `software/stage1_solver/runs/step2_manufactured_solutions/mms_localized_maxwell_h_equals_z/nr_24_nw_20/manifest.json`
- mms_localized_maxwell_h_equals_z nr_48_nw_40: `software/stage1_solver/runs/step2_manufactured_solutions/mms_localized_maxwell_h_equals_z/nr_48_nw_40/manifest.json`
- mms_localized_maxwell_h_equals_z nr_96_nw_80: `software/stage1_solver/runs/step2_manufactured_solutions/mms_localized_maxwell_h_equals_z/nr_96_nw_80/manifest.json`
- mms_localized_maxwell_h_equals_z nr_192_nw_160: `software/stage1_solver/runs/step2_manufactured_solutions/mms_localized_maxwell_h_equals_z/nr_192_nw_160/manifest.json`
- mms_localized_maxwell_h_equals_z nr_384_nw_320: `software/stage1_solver/runs/step2_manufactured_solutions/mms_localized_maxwell_h_equals_z/nr_384_nw_320/manifest.json`
- mms_wall_s_eta_2 nw_32: `software/stage1_solver/runs/step2_manufactured_solutions/mms_wall_s_eta_2/nw_32/manifest.json`
- mms_wall_s_eta_2 nw_64: `software/stage1_solver/runs/step2_manufactured_solutions/mms_wall_s_eta_2/nw_64/manifest.json`
- mms_wall_s_eta_2 nw_128: `software/stage1_solver/runs/step2_manufactured_solutions/mms_wall_s_eta_2/nw_128/manifest.json`
- mms_wall_s_eta_2 nw_256: `software/stage1_solver/runs/step2_manufactured_solutions/mms_wall_s_eta_2/nw_256/manifest.json`
- mms_wall_s_eta_2 nw_512: `software/stage1_solver/runs/step2_manufactured_solutions/mms_wall_s_eta_2/nw_512/manifest.json`

