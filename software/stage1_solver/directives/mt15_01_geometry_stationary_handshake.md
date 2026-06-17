# Directive — M1a: geometry + stationary open-throat machinery + V2-22B handshake (pre-freeze, placeholder)

You are the **CODER** (Codex codes / Claude reviews). Build the first piece of the Mathematica 15
branch-PRODUCER: a Wolfram-Language stationary open-throat geometry/profile module that emits a
**schema-valid V2-22B branch packet** and prove the *handshake* by pushing that packet through the EXISTING
Python validator/extractor chain. M0 (env check) is done and green; this is M1a.

**This step proves the PIPE, not the physics.** All physical values are deliberately **placeholder** and
**target-blind**. The modal/port content here is acknowledged-posited machinery scaffolding — it will be
**replaced in M1b** by values from a genuine coupled eigensolve. Do **not** read, report, or tune to the
target residuals; do **not** choose any placeholder by reference to the target constant
`54 G c_s^5 / (5 a^5 c^5)`. The point is: a packet we produce validates, its coefficients extract, and the
interface's target-leak guard is confirmed.

## Context you need

- The producer→judge interface is the frozen schema **`stage_v2_22b_solver_handoff/v1`**, validated by
  `research/pde_audit/scripts/stage_v2_22b_solver_handoff_validator.py`, then adapted+extracted by
  `…/stage_v2_22a_profile_to_coefficient_adapter.py` and `…/stage_v2_21_branch_extraction_fixture.py`,
  orchestrated end-to-end by `…/stage_v2_22c_end_to_end_smoke_pipeline.py`
  (`--solver-output <packet.json>`, `--tol 1e-9`).
- **Firewall:** the V2 chain is an EXTERNAL judge. The only interface is the exported JSON packet file. Do
  NOT import any V2/solver Python into the WL module or vice-versa. Do NOT touch
  `research/pde_audit/simulation/` (separate firewalled dir — not needed here).
- The validator **forbids** target fields in the packet (`P0_target, R_pole, R_norm, R_P2, R_P4, gamma_eff,
  gamma_GR, pass_flags, residuals, target_packet_pass, target_values`) — never emit any of these.

## Deliverables (under `software/stage1_solver/mathematica/`)

1. `mt15_01_geometry_stationary.wls` — builds the geometry + profiles, assembles the V2-22B packet,
   exports it to a gitignored run dir (e.g. `software/stage1_solver/mathematica/runs/`).
2. `mt15_01_handshake_report.md` — **markdown** (not JSON) report (see the YAML/markdown-not-JSON rule).

## What to build

### A. Geometry + grid (1-D throat coordinate `w`)
- A monotone grid on `w ∈ [0, L]`: strictly increasing, `points[0] = 0`, `points[-1] = L`, ≥ 2 points
  (use a real resolution, e.g. 64–129 points).
- A stationary **open-throat** radius profile `R0(w)` with `R0(0) = a` (mouth) and `R0(L) = R_exit > 0`
  (open exit — **NOT a hard cap**; `R_exit` strictly positive). For M1a use an explicit smooth monotone
  profile (e.g. a smooth interpolation `a → R_exit`); a FEM-solved reduced-energy profile is a later
  refinement, not required now. (`NDSolve`FEM`` and curvilinear ops are available if you want them.)
- A smooth positive weight/measure `W(w)` for the overlap integrals.
- Report geometry diagnostics: `R0(0)=a`, `R0(L)=R_exit>0`, monotonicity, grid endpoints.

### B. Placeholder modal/port content (acknowledged scaffolding)
- `wall.K > 0`, `wall.M > 0` (placeholder positive).
- `profiles.weight` and `profiles.wall_chi_eta` (length = grid length); normalize `wall_chi_eta` so the wall
  norm `∫ W · χ_η² dw ≈ 1` (avoid the validator's norm warning, tol 5e-3).
- ≥ 1 `bdg_modes` entry: `name`, finite `lambda_B`, `varpi > 0`, `profile_values` (length = grid),
  normalized so `∫ W · φ² dw ≈ 1`.
- ≥ 1 `mixed_ports` entry: `name`, finite `lambda_U/lambda_W/lambda_R`, `Omega_U > 0`, `Omega_W > 0`,
  `u_values`/`w_values` (length = grid), chosen so the validator's stability gate
  `Δ_eff = Omega_U² Omega_W² − (lambda_R · I_uw)² > 1e-12` holds, where `I_uw = ∫ W · u · w dw`.

### C. Freeze block + constants + metadata (honest labels)
- `freeze`: `pre_target_freeze=true`, `target_blind=true`, a declared nonempty `gauge_convention`
  (e.g. `"localized_H_equals_Z"`), `boundary_protocol="open_impedance_AC_reflecting_DC_leaking"`.
  Also declare `parent_action_status="effective_wall_closure"` (honest: the wall is an effective interface
  closure, no `S_Σ[R]` parent throat action — matches the project's frozen decision). Mark the packet as a
  placeholder/pre-freeze machinery export in `branch_id` and a metadata note.
- `geometry`: `L>0`, `R_mouth=a`, `R_exit>0`, `boundary_class="open_impedance"`, finite `Y_L_limit`,
  `exit_model` (e.g. `"impedance_mismatch_open_exit"`).
- `constants`: finite **sane placeholder** `G, c_s, c, a, mhat0, S_port, theta_tail` — chosen for sanity,
  NOT reverse-engineered from any target.
- `solver_metadata`: nonempty `exporter`, `coefficient_family`, `source_commit`; `mesh_points` = grid
  length; optional residual norms finite ≥ 0.
- `schema = "stage_v2_22b_solver_handoff/v1"`. **No forbidden target fields.**

### D. Handshake dry-run (prove the pipe)
- Export the packet JSON, then run the real `stage_v2_22c_end_to_end_smoke_pipeline.py` on it
  (`timeout 600`). Confirm and record in the report:
  - V2-22B **validation_pass = True** (and which open/stability gates passed: `R_exit>0`,
    `boundary_class`, `K>0`, `M>0`, `varpi>0`, `Δ_eff>0`, monotone grid, profile norms).
  - The V2-22A→V2-21 coefficient bundle **extracts** (`B/Z/N/D, u2/u4, P0/P2/P4` produced).
  - **Target-leak guard teeth-test:** also export a deliberately-bad copy with one forbidden field injected
    (e.g. `R_norm`) and confirm the validator **REJECTS** it (proves the guard is live, not assumed).
- **LOUD scoping (mandatory):** state in the report that any `R_pole/R_norm/R_P2/R_P4` printed by the
  pipeline on this placeholder packet are **machinery-proof only, physically meaningless, NOT a Stage-1
  result, and were NOT tuned to any target.** Do not present them as findings.

## Working rules
- ≤ 2 concurrent `math`/`wolframscript` kernels (this step uses 1). `timeout 600` around every script RUN;
  iterate WL export + V2 run until both exit 0. Never wrap your own session in a shell `timeout`.
- Files only under `software/stage1_solver/mathematica/`; packet JSON → a gitignored `runs/` (or `_scratch/`)
  subdir. Do **not** `git add` or commit. No network, no GPU.

## Report back to the reviewer
- The geometry diagnostics; the exact packet structure you emitted (sections + key example values);
  confirmation V2-22B `validation_pass=True` with the gate list; the extracted coefficient bundle exists;
  the target-leak teeth-test (bad packet rejected, with the validator's error); and the LOUD placeholder
  scoping. Then Claude reviews (schema-fidelity + target-blindness/honesty) and an arbiter re-runs it.
