# Path-A C0g B-3 Deep-Point Characterization

## Contract

- Deep-point battery only; no remedy tool was run.
- Fixed-tau convergence is judged only by the original `patha_closed_branch_residual` Linf norm and update norm.
- B-4 autodiff sparse Jacobian is used only for Newton directions and SVD diagnostics.

## Stage 1 - Deep Re-Converge

- reading: `BRACKETED_NEAR_SINGULARITY` / `monotone_tight_to_non_tight_as_tau_deepens`
- target tau: `2.911625e-02`
- first non-tight tau: `2.911594e-02`

| tau | label | status | final residual | final step | iterations |
|---:|---|---|---:|---:|---:|
| 2.911625e-02 | `stage1_B3_followup_C1_recorded_anchor_tau_0p02911625` | `TIGHT_CONVERGED` | 2.052413e-13 | 2.142338e-11 | 0 |
| 2.911594e-02 | `refine_02_tau_0p02911594375` | `STALLED` | 2.224632e-08 | 4.590617e-02 | 4 |
| 2.911564e-02 | `refine_01_tau_0p0291156375` | `STALLED` | 1.605095e-07 | 1.386098e-01 | 5 |
| 2.911502e-02 | `refine_00_tau_0p029115025` | `STALLED` | 3.484567e-07 | 6.909035e-01 | 2 |
| 2.911380e-02 | `stage1_deepcrawl_state_dir_attempt_024_tau_0p0291138` | `STALLED` | 8.053571e-07 | 3.700520e-01 | 2 |
| 2.911320e-02 | `stage1_deepcrawl_state_dir_attempt_025_tau_0p0291132` | `STALLED` | 9.795500e-07 | 4.426295e-01 | 1 |
| 2.911310e-02 | `stage1_deepcrawl_state_dir_attempt_028_tau_0p0291131` | `STALLED` | 1.026576e-06 | 4.202788e-01 | 2 |
| 2.911300e-02 | `stage1_deepcrawl_state_dir_attempt_027_tau_0p029113` | `STALLED` | 1.085899e-06 | 5.676729e-01 | 3 |
| 2.911280e-02 | `stage1_deepcrawl_state_dir_attempt_026_tau_0p0291128` | `STALLED` | 1.184823e-06 | 6.736036e-01 | 3 |
| 2.911265e-02 | `stage1_deepcrawl_state_dir_attempt_023_tau_0p0291126455688` | `STALLED` | 1.184761e-06 | 4.177653e-01 | 1 |
| 2.911086e-02 | `stage1_deepcrawl_state_dir_attempt_022_tau_0p0291108575439` | `STALLED` | 1.800293e-06 | 1.736098e+00 | 0 |

## Stage 2 - C2/C3 At Deep Point

- status: `MEASURED`
- smallest singular values: `[0.17269727105522853, 0.1557835137096228, 0.11134898043597927, 0.10298097396919434, 0.0852268495813875, 0.05998703863353096, 0.058762385126102136, 0.024105839877529244, 0.01723150165432018, 7.80413975247924e-05]`
- cluster call: `ISOLATED_NEAR_NULL`
- sector fractions: `{'psi': 0.0011568954911383013, 'A': 3.8515697900584036e-08, 'r0': 0.9940375383256637, 'mu': 0.00480552766750007}`
- fold transversality: `{'wT_Ftau': 0.4859017832482557, 'normalized_abs_wT_Ftau': 0.9175629363132609, 'call': 'FOLD_TRANSVERSAL'}`
- gauge tests: `{'scaled_removed_gauge_projection_fraction': 2.5227656348386824e-33, 'independent_expanded_candidate_projection_fraction': 0.8109343699010971, 'curl_mode_l2': 3.368169681690008e-15, 'scaled_residual_response_l2': 7.804139752489914e-05, 'high_physical_overlap_alone_is_not_used': True}`
- localization: `{'throat_energy_fraction_first_quarter_w': 0.022155079346926575, 'classification': 'EXTENDED', 'peak_cell': {'r_index': 9, 'w_index': 15}}`
- one-rung comparison: `FOLD_SIGNATURE_SHARPENS_OR_PERSISTS`
- C3 line scan: status `SKIPPED_NO_TIGHT_NEIGHBOR`, call `None`, max normalized `n/a`

## Stage 3 - C4 Resolution Relocation

- status: `LIMITED_C4_INSUFFICIENT_TIGHT_RELOCATION`
- omission log: `32x32 deferred by configuration; C4 verdict is candidate/limited only.`
- 16x16 reference: `{'tau': 0.02911625, 'raw_scaled_sigma_min_feature': 7.80413975247924e-05, 'sigma_min_over_sigma_max_feature': 7.946435710844274e-06, 'feature_mode_sector': {'psi': 0.0011568954911383013, 'A': 3.8515697900584036e-08, 'r0': 0.9940375383256637, 'mu': 0.00480552766750007}, 'feature_localization': {'throat_energy_fraction_first_quarter_w': 0.022155079346926575, 'classification': 'EXTENDED', 'peak_cell': {'r_index': 9, 'w_index': 15}}}`
- grid `[24, 24]` status `INSUFFICIENT_TIGHT_RELOCATION` feature tau `n/a`, sigma ratio `n/a`, collapse `n/a`, localization `None`

## Tool Selection

- case: `CANDIDATE_ONLY_C4_NOT_DECISIVE`
- tool: `STOP_NO_REMEDY_TOOL`
- why: C4 status LIMITED_C4_INSUFFICIENT_TIGHT_RELOCATION does not provide a tight same-grid relocation, so no physical/discretization case is finalized

## Scope Guard

- `{'single_arbiter_residual': 'stage1_solver.coupled_branch.patha_closed_branch_residual', 'convergence_judged_only_on_original_residual': True, 'gauge_Q_perp_rebuilt_inside_every_fixed_tau_solve': True, 'b4_autodiff_sparse_jacobian_direction_and_svd_only': True, 'patha_closed_branch_residual_touched': False, 'faithful_operators_touched': False, 'preconditioners_touched': False, 'xi_or_grad_div_penalty_touched': False, 'physical_export_permitted_touched': False, 'frozen_physics_touched': False, 'no_pseudoarclength_or_LM_or_deflation_remedy_run': True, 'comparison_artifacts_diagnostics_only': True}`
- git diff summary: `No diff in src/stage1_solver/coupled_branch.py, src/stage1_solver/operators.py, or src/stage1_solver/preconditioners.py.`

## Artifacts

- JSON: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_deeppoint/pathA_C0g_B3_deeppoint.json`
- progress: `/var/projects/toy_physics/software/stage1_solver/runs/pathA_C0g_B3_deeppoint/progress.jsonl`
