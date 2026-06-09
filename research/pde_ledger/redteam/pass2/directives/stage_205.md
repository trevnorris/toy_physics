---
unit_id: 205
batch: VI.1
created_at: 2026-06-09T21:25:48Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 205

There is exactly one finding, and it is a `paper_misalignment` pending user resolution.
**Codex applies nothing on this unit.** Do not edit `paper/`, `notes/`, or the scripts.
The orchestrator is holding for the user to choose a direction.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_205.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- The artifact in question is `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl` (9367 bytes), which exists and passes all eleven checks M1-M11 (see `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.txt`, "STAGE 205 MATHEMATICA AUDIT PASSED").

## Resolve before fix_loop

The Stage 205 card's `\stagefield{Verification}` line says "Mathematica audit: none yet," but a passing Mathematica audit `.wl` was added in the pass-1 dual-engine retrofit. The card under-reports the verification coverage that now exists. Which is correct?

Possible directions (the user picks one):
- (a) The `.wl` is intended → update `paper/stages/stage_205.tex:11` to cite `\StageFile{mathematica/moving_throat_pde_stage205_directional_hessian_and_quadratic_root_refinement_mathematica_audit.wl}` in place of "Mathematica audit: none yet." (This is the expected direction; it should also flow into the post-batch MATHEMATICA_MIRROR_POLICY tracker sync.) Paper-side edit is user-owned; Codex does not apply it.
- (b) The `.wl` should not exist for this stage → quarantine/remove the script and leave the card unchanged. (Unlikely given the dual-engine-required policy, but it is the user's call.)

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No math change is implied by either direction — the scripts are correct and dual-engine-independent.
