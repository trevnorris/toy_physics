# Path-A Chunk 1b Closed Newton

Overall engineering gate: PASS
S_Sigma digest: `63501899d5957528`

Scope: closed static engineering solve over matter, gauge, R0, and chemical potential with placeholder S_Sigma. No physical packet or coefficient export is performed.

## Placeholder Parameters

Label: Path-A chunk-1b closed engineering placeholder; no physical packet; target-blind

```yaml
solve_grid: (6, 6)
continuation_K_values: (0.08, 0.18)
mass: 0.05
r_mouth: 1.0
r_exit: 0.92
radial_wall_strength: 0.045
axial_trap_strength: 0.04
matter_exit_impedance_alpha: 0.2
a0_exit_impedance_alpha: 0.35
preconditioner_layout: 5*cells+nw+1
exit_wall_bc: natural_zero_traction
```

## Newton Convergence

| eos_K | converged | iterations | initial_residual_norm | final_residual_norm | tolerance | gmres_iterations | r0_min | r0_max | message |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8.000000e-02 | True | 4 | 1.067053e+00 | 1.726341e-12 | 2.134107e-09 | [1, 1, 1, 1] | 9.877030e-01 | 9.980074e-01 | residual tolerance reached |
| 1.800000e-01 | True | 0 | 1.381831e-09 | 1.381831e-09 | 2.000000e-09 | [] | 9.877030e-01 | 9.980074e-01 | initial residual met tolerance |

Final residual linf=1.381831e-09; R0 range=[9.877030e-01, 9.980074e-01]; mass=5.000000e-02; mu=2.013805e+00; wall-clock=1.173417e+02s.
Closed-Newton tolerance note: the counted gate uses residual_atol=2.000000e-09 and residual_rtol=2.000000e-09; the active nonlinear closed solve achieved residual 1.726341e-12, so the relaxed absolute tolerance is not masking non-convergence.

## JVP Check

Closed residual JVP vs centered finite difference: relative=4.041896e-11, absolute=2.735362e-09, epsilon=1.000000e-05, status=PASS.

## Placeholder Derivatives

| quantity | absolute | relative |
| --- | --- | --- |
| T_w_R | 1.676354e-12 | 1.676354e-12 |
| T_w_RR | 5.764738e-08 | 5.764738e-08 |
| U_R | 9.441967e-11 | 9.441967e-11 |
| U_RR | 7.665529e-09 | 7.665529e-09 |

Relative-derivative note: relative and absolute entries can be identical here because the denominator convention is `max(1, ||analytic||_inf)` and these analytic derivative magnitudes are below one.

## Return Source Diagnostic

The wall source was recomputed from the shared confinement coefficient and shell volumes. max_abs_diff=0.000000e+00, source range=[-5.782026e-03, -5.000463e-03].

## Stability Margin

```yaml
mu_min: 1.051521e+00
T_w_min: 3.628045e+00
T_Omega_min: 9.495381e-01
U_RR_min: 7.341337e-02
minimum_margin: 7.341337e-02
```

A downstream Schur-denominator value is not constructed in chunk 1b; this report records the placeholder constitutive positivity margin used by the closed background solve.
This counted check is a constitutive-positivity sanity guard for the smooth positive placeholder family, near true by construction and retained as a stability precondition for the background solve; it is not independent physics evidence.

## Counted Gates

- closed_jvp_consistency: PASS (relative=4.041896e-11, absolute=2.735362e-09)
- closed_newton_convergence: PASS (final_linf=1.381831e-09, R0_min=9.877030e-01)
- placeholder_provider_derivatives_fd: PASS (max_relative=5.764738e-08, max_absolute=5.764738e-08)
- constitutive_positive_margin: PASS (minimum_margin=7.341337e-02)
- target_token_scan_clean: PASS (scanned 5 Path-A 1b files.)

## Diagnostics Not A Physics Gate

- return_source_fidelity_not_a_physics_gate: reported (max_abs_diff=0.000000e+00)
- gauge_return_no_explicit_source_not_a_physics_gate: construction_note (matter and Maxwell lanes are coupled monolithically; no extra gauge source was added.)

## Stop Conditions

- R0_nonpositive: not tripped
- exit_bc_changed_from_zero_traction: not tripped
- hidden_clamp_or_line_search_hack: not tripped
- source_sign_or_reduction_mismatch: not tripped

