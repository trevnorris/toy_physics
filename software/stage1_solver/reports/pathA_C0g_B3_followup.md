# Path-A C0g B-3 Follow-Up v2

## Contract

- Battery: C0 read-offs, C1 tight two-start reconvergence, C2/C3 gated characterization, C4 gated resolution check.
- Tool selection is evidence-driven; pseudo-arclength is not pre-selected.
- Original `patha_closed_branch_residual` remains the single arbiter.

## C0 Read-Offs

- sigma trend call: `STILL_DROPPING_STEEPLY`
- deepest converged tau: `2.911320e-02`
- deepest sigma_min(JQ_perp): `7.026590e-06`
- fit zero crossing: `2.911390e-02`
- zero-crossing minus deepest tau: `7.000000e-07`

| tau | sigma_min_JQ_perp | cond_JQ_perp | status |
| --- | --- | --- | --- |
| 2.912250e-02 | 3.948475e-04 | - | - |
| 2.912200e-02 | 3.792773e-04 | - | - |
| 2.912000e-02 | 3.095487e-04 | - | - |
| 2.911813e-02 | 2.254878e-04 | - | - |
| 2.911625e-02 | 1.455266e-04 | - | - |
| 2.911443e-02 | 1.182048e-04 | - | - |
| 2.911380e-02 | 7.885724e-05 | - | - |
| 2.911320e-02 | 7.026590e-06 | - | - |

### Arclength Metric Balance

| id | previous_tau | tau | dx_norm_squared | tau_metric_squared | x_fraction | dominant |
| --- | --- | --- | --- | --- | --- | --- |
| B3_0_accepted_001 | 2.912201e-02 | 2.912134e-02 | 3.899405e-07 | 4.395602e-07 | 4.700906e-01 | tau |
| B3_0_accepted_002 | 2.912134e-02 | 2.912067e-02 | 4.581431e-07 | 4.559466e-07 | 5.012015e-01 | x |
| B3_0_accepted_003 | 2.912067e-02 | 2.912000e-02 | 5.150940e-07 | 4.443354e-07 | 5.368753e-01 | x |
| B3_0_accepted_004 | 2.912000e-02 | 2.911938e-02 | 5.313983e-07 | 3.905267e-07 | 5.764008e-01 | x |
| B3_0_accepted_005 | 2.911938e-02 | 2.911875e-02 | 6.389945e-07 | 3.903480e-07 | 6.207793e-01 | x |
| B3_0_accepted_006 | 2.911875e-02 | 2.911813e-02 | 8.017420e-07 | 3.900575e-07 | 6.727155e-01 | x |
| B3_0_accepted_007 | 2.911813e-02 | 2.911750e-02 | 1.076249e-06 | 3.892655e-07 | 7.343830e-01 | x |
| B3_0_accepted_008 | 2.911750e-02 | 2.911688e-02 | 1.636348e-06 | 3.863887e-07 | 8.089773e-01 | x |
| B3_0_accepted_009 | 2.911688e-02 | 2.911626e-02 | 3.658220e-06 | 3.895732e-07 | 9.037566e-01 | x |

## C1 Tight Reconvergence

- case: `Case 1`
- reading: `both_tight_same`
- recorded status: `TIGHT_CONVERGED` residual `2.052413e-13` step `2.142338e-11`
- continuation status: `TIGHT_CONVERGED` residual `2.142822e-13` step `7.320407e-12`

Gauge-invariant mismatch after phase alignment:

```yaml
{'rho_linf': 1.6214460329955216e-14, 'curl_A_linf': 7.228960273916914e-19, 'r0_linf': 5.024980431755921e-12, 'mu_abs': 9.72111280361787e-13}
```

## C2 Mode Characterization

status: `MEASURED`

### recorded

- sigma_min: `7.804140e-05`
- near-null count call: `ISOLATED_NEAR_NULL` gap `2.207995e+02`
- dominant sector fractions: `{'psi': 0.0011568954911383013, 'A': 3.8515697900584036e-08, 'r0': 0.9940375383256637, 'mu': 0.00480552766750007}`
- transversality: `{'wT_Ftau': 0.4859017832482557, 'normalized_abs_wT_Ftau': 0.9175629363132609, 'call': 'FOLD_TRANSVERSAL'}`
- gauge tests: `{'scaled_removed_gauge_projection_fraction': 2.5227656348386824e-33, 'independent_expanded_candidate_projection_fraction': 0.8109343699010971, 'curl_mode_l2': 3.368169681690008e-15, 'scaled_residual_response_l2': 7.804139752489914e-05, 'high_physical_overlap_alone_is_not_used': True}`
- localization: `{'throat_energy_fraction_first_quarter_w': 0.022155079346926575, 'classification': 'EXTENDED', 'peak_cell': {'r_index': 9, 'w_index': 15}}`

### continuation

- sigma_min: `7.804140e-05`
- near-null count call: `ISOLATED_NEAR_NULL` gap `2.207995e+02`
- dominant sector fractions: `{'psi': 0.00115689549092288, 'A': 3.8515697894105776e-08, 'r0': 0.9940375383262159, 'mu': 0.004805527667163325}`
- transversality: `{'wT_Ftau': 0.4859017832364471, 'normalized_abs_wT_Ftau': 0.9175629363168817, 'call': 'FOLD_TRANSVERSAL'}`
- gauge tests: `{'scaled_removed_gauge_projection_fraction': 2.1775715599688314e-33, 'independent_expanded_candidate_projection_fraction': 0.8109343699020704, 'curl_mode_l2': 3.173633255722171e-15, 'scaled_residual_response_l2': 7.804139703528706e-05, 'high_physical_overlap_alone_is_not_used': True}`
- localization: `{'throat_energy_fraction_first_quarter_w': 0.02215507934700143, 'classification': 'EXTENDED', 'peak_cell': {'r_index': 9, 'w_index': 15}}`

## C3 Line Scan

- status: `MEASURED`
- call: `FLAT_OR_NO_CLEAR_SPIKE`
- max normalized residual: `1.000001e+00`

## C4 Resolution Check

- status: `SKIPPED_BY_DEFAULT_GATE`
- reason: `resolution relocation is expensive; Case-1 stops for re-gate before any tool claim, so fold-vs-wall remains open`

## Tool Selection

- case: `Case 1`
- chosen tool: `STOP_RECONVERGE_DEEP_STATES_AND_GATE_REFINE_THEN_C2_C4`
- why: C1 found the same tight root at this tau; global fold-vs-wall remains open

## Scope Guard

```yaml
single_arbiter_residual: stage1_solver.coupled_branch.patha_closed_branch_residual
convergence_judged_only_on_original_residual: True
gauge_Q_perp_inside_every_C1_solve: True
patha_closed_branch_residual_touched: False
faithful_operators_touched: False
xi_or_grad_div_penalty_touched: False
physical_export_permitted_touched: False
frozen_physics_touched: False
no_LM_no_PTC_no_Sobolev_regularization: True
no_pseudoarclength_precommit: True
comparison_artifacts_diagnostics_only: True
```

## Git Diff Summary

```
No diff in src/stage1_solver/coupled_branch.py or src/stage1_solver/operators.py.
```

Machine JSON: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_followup/pathA_C0g_B3_followup.json`
Progress log: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_followup/progress.jsonl`
