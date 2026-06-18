# Path-A Chunk 1a Static Balance MMS Report

Overall gate: PASS
Config hash: `62db24d3532de346`
S_Sigma spec digest: `99a48c7c8163a213`

## Built Surface

Implemented a serializable `S_Sigma` spec/registry and a separate nonlinear flat-`dw` FV static-balance operator for `R0(w)` with Dirichlet face values.

## Genuine Gates

| check | pass | detail |
| --- | --- | --- |
| dual_engine_forcing_agreement | True | SymPy and Mathematica source samples agree to 8.438e-15. |
| second_order_static_balance_mms | True | Discrete LHS compared with independently derived continuum source on a refinement ladder. |
| flux_divergence_term_nonzero | True | final L_inf=3.792e+00. |
| gradient_square_term_nonzero | True | final L_inf=1.737e-02. |
| potential_gradient_term_nonzero | True | final L_inf=6.750e-01. |
| target_token_scan_clean | True | scanned 5 chunk-1a code/report files. |

## Construction Restatements

| check | pass | detail |
| --- | --- | --- |
| spec_roundtrip_not_a_physics_gate | True | Serializable registry selector round-trips without callables. |
| dirichlet_end_values_not_a_physics_gate | True | Chunk-1a MMS prescribes both manufactured end values; open-exit choice is deferred. |

## MMS Refinement

Path-A chunk-1a nonlinear static throat balance with R-dependent T_w, live gradient-square term, and nonzero U_R.
Continuum source: `-d_w(T_w(R0,w) d_w R0) + 0.5 T_w_R(R0,w)(d_w R0)^2 + U_R(R0,w)`
Manufactured field: `7.99666805497709*w**4*(1 - 0.714285714285714*w)**4 + 0.025*sin(0.714285714285714*pi*w)**4 + 1.08`
Forcing derivation: SymPy differentiates the continuum balance analytically; Mathematica exports an independent sample check of the same continuum expression.

| grid | spacing | error | observed_order | reference_norm | residual_linf | flux_divergence_linf | gradient_square_linf | potential_gradient_linf |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nw_32 | 4.375000e-02 | 1.924502e-02 | - | 2.539571e+00 | 2.709328e-02 | 3.764483e+00 | 1.677431e-02 | 6.746310e-01 |
| nw_64 | 2.187500e-02 | 5.039107e-03 | 1.933245e+00 | 2.539577e+00 | 8.076745e-03 | 3.784037e+00 | 1.723988e-02 | 6.749498e-01 |
| nw_128 | 1.093750e-02 | 1.309989e-03 | 1.943613e+00 | 2.539577e+00 | 2.556594e-03 | 3.790729e+00 | 1.734761e-02 | 6.750296e-01 |
| nw_256 | 5.468750e-03 | 3.352718e-04 | 1.966152e+00 | 2.539577e+00 | 7.100674e-04 | 3.791534e+00 | 1.736998e-02 | 6.750496e-01 |

Finest-grid summary: error=3.352718e-04, observed_order=1.966152e+00.

## Dual Engine

SymPy and Mathematica independently evaluated the continuum manufactured source samples; max absolute difference = 8.437695e-15.
Because the finest-grid gradient-square L_inf (1.736998e-02) is small relative to flux-divergence L_inf (3.791534e+00), that term's correctness is validated primarily by the dual-engine continuum-term agreement (`term_max_abs_diff.gradient_square` = 1.092876e-16) and its live residual contribution, not by the global second-order residual metric alone.
Mathematica diagnostics: `/var/projects/toy_physics/software/stage1_solver/mathematica/runs/pathA_02_chunk1a_static_balance_mms/pathA_02_chunk1a_static_balance_mms_diagnostics.json`

## Validation Scope

This MMS exercises only the `patha_static_mms_v1` provider; the `smooth_positive_placeholder_v1` family derivatives (`T_w_R`, `T_w_RR`, `U_R`, `U_RR`) are analytically asserted and hand-verified, with MMS validation of the actual solve family deferred to chunk 1b/1c.

## Under-Specified Items

None for chunk 1a. The background subtraction, radial reduction, and open-exit BC choices remain deferred to chunk 1b and are not used by this MMS.

