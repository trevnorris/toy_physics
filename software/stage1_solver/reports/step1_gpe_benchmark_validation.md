# Step 1 GPE Benchmark Validation

Overall gate: PASS
Config hash: `06d1d0eca45b711e`

## Solver Controls

```yaml
dtype: float64
device: cpu
seed: 20260612
deterministic_algorithms: True
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
finite_difference_jvp_epsilon: 1e-06
```

## Discretization

Cell-centered finite volume on radial shells. Radial divergence is `(F_{i+1/2}-F_{i-1/2})/V_i` with face flux `F=4*pi*r_face^2*d_r psi`; the inner face has zero area, enforcing the regular r=0 flux condition. The reusable tensor-product `(r,w)` operator extends this with cell volume `shell_volume*dw`, radial flux `4*pi*r_face^2*dw*d_r psi`, and w-flux `shell_volume*d_w psi`. Dirichlet and Robin boundaries use the same ghost-cell boundary operator.

## Linear Benchmark

Linear radial harmonic oscillator ground state with exact Robin boundary.
Reference: closed form exp(-omega*r^2/2), mu=3*omega/2.

| nr | dr | computed_mu | eigenvalue_error | eigenvalue_order | field_l2_error | field_l2_order | discrete_eigen_residual_linf | mass_drift | origin_flux_abs | current_max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 128 | 6.250000e-02 | 1.499064e+00 | 9.357697e-04 | - | 5.108853e-04 | - | 8.044676e-13 | 1.110223e-16 | 0.000000e+00 | 0.000000e+00 |
| 256 | 3.125000e-02 | 1.499766e+00 | 2.339617e-04 | 1.999881e+00 | 1.276785e-04 | 2.000484e+00 | 5.300336e-12 | 2.220446e-16 | 0.000000e+00 | 0.000000e+00 |
| 512 | 1.562500e-02 | 1.499942e+00 | 5.849162e-05 | 1.999971e+00 | 3.191694e-05 | 2.000121e+00 | 1.932558e-11 | 3.330669e-16 | 0.000000e+00 | 0.000000e+00 |
| 1024 | 7.812500e-03 | 1.499985e+00 | 1.462299e-05 | 1.999992e+00 | 7.979067e-06 | 2.000030e+00 | 5.818346e-11 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 |

Linear checks:
- eigen_order: PASS
- field_order: PASS
- final_eigenvalue_error: PASS
- final_field_l2_error: PASS
- mass_drift: PASS
- current_abs: PASS

## Cubic GPE Benchmark

Mass-constrained stationary cubic GPE in a radial harmonic trap.
Reference: independent SciPy solve_bvp radial BVP with unknown chemical potential.
Reference details: mu=1.560971500350e+00, mass=9.999999999747e-01, nodes=1304, status=0.

| nr | dr | computed_mu | mu_error | mu_order | field_l2_error | field_l2_order | newton_iterations | newton_final_residual_linf | mass_drift | origin_flux_abs | current_max_abs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 64 | 1.250000e-01 | 1.557533e+00 | 3.438740e-03 | - | 2.233743e-03 | - | 4 | 4.157276e-14 | 0.000000e+00 | 0.000000e+00 | 0.000000e+00 |
| 128 | 6.250000e-02 | 1.560111e+00 | 8.601407e-04 | 1.999236e+00 | 5.583090e-04 | 2.000328e+00 | 4 | 3.697661e-14 | 2.220446e-16 | 0.000000e+00 | 0.000000e+00 |
| 256 | 3.125000e-02 | 1.560756e+00 | 2.150636e-04 | 1.999809e+00 | 1.398559e-04 | 1.997122e+00 | 4 | 3.957945e-14 | 2.220446e-16 | 0.000000e+00 | 0.000000e+00 |

Jacobian JVP cross-check: relative=9.849034e-11, absolute=1.317892e-07, epsilon=1.000000e-06.

Cubic checks:
- mu_order: PASS
- field_order: PASS
- final_mu_error: PASS
- final_field_l2_error: PASS
- mass_drift: PASS
- current_abs: PASS
- jacobian_relative: PASS
- jacobian_absolute: PASS

## Manifests

- linear_harmonic_robin nr=128: `software/stage1_solver/runs/step1_gpe_benchmark/linear_harmonic_robin/nr_128/manifest.json`
- linear_harmonic_robin nr=256: `software/stage1_solver/runs/step1_gpe_benchmark/linear_harmonic_robin/nr_256/manifest.json`
- linear_harmonic_robin nr=512: `software/stage1_solver/runs/step1_gpe_benchmark/linear_harmonic_robin/nr_512/manifest.json`
- linear_harmonic_robin nr=1024: `software/stage1_solver/runs/step1_gpe_benchmark/linear_harmonic_robin/nr_1024/manifest.json`
- cubic_gpe_harmonic_mass nr=64: `software/stage1_solver/runs/step1_gpe_benchmark/cubic_gpe_harmonic_mass/nr_64/manifest.json`
- cubic_gpe_harmonic_mass nr=128: `software/stage1_solver/runs/step1_gpe_benchmark/cubic_gpe_harmonic_mass/nr_128/manifest.json`
- cubic_gpe_harmonic_mass nr=256: `software/stage1_solver/runs/step1_gpe_benchmark/cubic_gpe_harmonic_mass/nr_256/manifest.json`

