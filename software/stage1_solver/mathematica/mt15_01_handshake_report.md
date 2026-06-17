# MT15 M1a Geometry Stationary Handshake Report

Generated on 2026-06-17 from `software/stage1_solver/mathematica/mt15_01_geometry_stationary.wls`.

This is a pre-freeze, target-blind machinery handoff. The exported modal and port values are placeholder scaffolding only; they are not a Stage-1 physics result and are intended to be replaced in M1b by a genuine coupled eigensolve.

## Producer Output

Run command:

```bash
timeout 600 wolframscript -script software/stage1_solver/mathematica/mt15_01_geometry_stationary.wls
```

Valid packet:

`software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22b_placeholder_packet.json`

Bad target-leak teeth-test packet:

`software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22b_forbidden_R_norm_packet.json`

The run directory is gitignored by `.gitignore:89` via `software/**/runs/`.

## Geometry Diagnostics

- Coordinate: `w` on `[0, L]`.
- `L = 2.0`, mesh points = `129`.
- Grid endpoints: first `0.0`, last `2.0`; strictly increasing = `True`.
- Radius profile: smooth monotone cubic interpolation.
- `R0(0) = a = 1.0`.
- `R0(L) = R_exit = 1.65`, so the exit is open and strictly positive.
- Radius monotone nondecreasing = `True`; strictly increasing on the grid = `True`.
- Weight/measure range: min `1.0`, max `2.7225`.

Profile diagnostics from the WL export and V2 validator:

- Wall norm `integral W chi_eta^2 dw = 1.0000000000000004` in the validator report.
- BdG profile norm `= 0.9999999999999997` in the validator report.
- Mixed overlap `I_u_w = 0.9417803760099119`.
- Mixed determinant `Delta_eff = 169.49166281103697`, positive and above `1e-12`.

## Packet Structure

Top-level sections emitted:

`schema`, `branch_id`, `freeze`, `geometry`, `constants`, `grid`, `wall`, `profiles`, `bdg_modes`, `mixed_ports`, `solver_metadata`, `normalization_tolerances`.

Key values:

- `schema = "stage_v2_22b_solver_handoff/v1"`.
- `branch_id = "mt15_m1a_pre_freeze_placeholder_open_throat_machinery"`.
- Freeze flags: `pre_target_freeze = true`, `target_blind = true`.
- `gauge_convention = "localized_H_equals_Z"`.
- `boundary_protocol = "open_impedance_AC_reflecting_DC_leaking"`.
- `parent_action_status = "effective_wall_closure"`.
- Geometry: `boundary_class = "open_impedance"`, `exit_model = "impedance_mismatch_open_exit"`, `Y_L_limit = 0.0`.
- Constants: `G = 1.0`, `c_s = 1.0`, `c = 2.5`, `a = 1.0`, `mhat0 = 0.85`, `S_port = 1.2`, `theta_tail = 1.0`.
- Wall placeholders: `K = 2.25`, `M = 0.72`.
- One BdG mode: `name = "m1a_placeholder_bdg_halfwave"`, `lambda_B = 0.34`, `varpi = 2.6`.
- One mixed port: `name = "m1a_placeholder_open_port"`, `lambda_U = 0.24`, `lambda_W = 0.31`, `lambda_R = 0.18`, `Omega_U = 3.1`, `Omega_W = 4.2`.
- Metadata: `exporter = "mt15_01_geometry_stationary.wls"`, `coefficient_family = "v2_22b_stationary_open_throat_placeholder_machinery_v1"`, `source_commit = "2c6a9be"`, `mesh_points = 129`.

A recursive key scan of the valid packet found none of the forbidden target-output keys:

`P0_target`, `R_pole`, `R_norm`, `R_P2`, `R_P4`, `gamma_eff`, `gamma_GR`, `pass_flags`, `residuals`, `target_packet_pass`, `target_values`.

## External V2 Handshake

Run command:

```bash
timeout 600 python research/pde_audit/scripts/stage_v2_22c_end_to_end_smoke_pipeline.py --solver-output software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22b_placeholder_packet.json --tol 1e-9 --out-report software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22c_pipeline_report.json --out-profile-manifest software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22a_profile_manifest.json --out-v21-manifest software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_21_manifest.json --out-observable-packet software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22c_observable_packet.json --out-tolerance-budget software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22c_tolerance_budget.json
```

Result: exit code `0`.

V2-22B validation:

- `validation_pass = True`.
- `error_count = 0`, `warning_count = 0`.
- Open gates passed: `R_exit > 0`, `boundary_class = open_impedance`, hard-cap absent.
- Positivity gates passed: `wall.K > 0`, `wall.M > 0`, `varpi > 0`, `Omega_U > 0`, `Omega_W > 0`.
- Stability gate passed: `Delta_eff > 1e-12`.
- Grid/profile gates passed: monotone grid, endpoint match, wall norm and BdG norm within tolerance.

Pipeline extraction:

- `mechanical_pipeline_pass = True`.
- V2-22A profile manifest generated: hash `70da73b7f662026843f3d9574b21a746e11b57a1c36476b315adf5524a027e38`.
- V2-21 manifest generated: hash `e878bd1b1ba131ade955a9562bb05562c1585c7089baf9d81468280f27f007c1`.
- Observable packet generated: hash `9ede1584f04206c21542c12d61cb021cc6c9c34334d1d428a318ee045dd9cb81`.
- Coefficient bundle exists with grouped keys including `B0/B2/B4`, `Z0/Z2/Z4`, `N0/N2/N4`, `D0/D2/D4`, `u2/u4`, and `P0/P2/P4`.
- The pipeline report records `branch_target_realization_claimed = False`.

## Target-Leak Guard Teeth-Test

Run command:

```bash
timeout 600 python research/pde_audit/scripts/stage_v2_22b_solver_handoff_validator.py --solver-output software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_v2_22b_forbidden_R_norm_packet.json --out-report software/stage1_solver/mathematica/runs/mt15_01_geometry_stationary/mt15_01_bad_R_norm_validation_report.json
```

The validator rejected the deliberately bad packet:

- `validation_pass = False`.
- `error_count = 1`, `warning_count = 0`.
- Error: `R_norm: target residual/output fields are forbidden in frozen solver exports`.

This confirms the target-leak guard is active.

## Loud Placeholder Scoping

Any `R_pole`, `R_norm`, `R_P2`, or `R_P4` printed by the downstream V2 pipeline for this placeholder packet are machinery-proof only, physically meaningless, NOT a Stage-1 result, and were NOT tuned to any target. They are recorded only by the external judge after JSON handoff to prove the pipe validates, adapts, and extracts. The WL producer does not import the V2 Python chain and does not emit target residuals or target values.
