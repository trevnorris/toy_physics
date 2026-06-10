---
unit_id: 220
batch: VII.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 220

The only finding for this unit is a `paper_misalignment` (F1) requiring user resolution.
Codex must do **nothing** here: do not edit `paper/`, `notes/`, or any script. The
orchestrator is holding for the user to choose a direction. No Codex-applied changes
are specified. The scripts already pass in both engines (committed outputs all green);
there is no script-side fix to make.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (documentation understates verification: a passing second engine exists but the card/notes deny it)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_220.tex:11` quote:
  `\stagefield{Verification}{SymPy audit: \StageFile{...stage220...sympy_audit.py}.  Mathematica audit: none yet.}`
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_sympy_audit.md:561` quote:
  "The accompanying SymPy audit verifies:" (Section 10) — and the "Supporting file" list (line 592) names only the `.py`.

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl` exists, is committed, and passes all nine M-blocks. Committed output `mathematica/output/...stage220...mathematica_audit.txt:91` ends "All Stage 220 Mathematica checks passed."
- For convention, the neighbouring card already cites its `.wl` correctly: `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_221.tex:11`:
  `... Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage221_..._mathematica_audit.wl}.`

## Resolve before fix_loop

Stage 220's card (line 11) and notes (Section 10) state the stage is SymPy-only
("Mathematica audit: none yet"), but a committed, passing Mathematica audit `.wl`
exists and dual-engine-verifies every deliverable. Which direction is correct?

Possible directions (the user picks one):
- (a) The `.wl` is the intended second engine (it passes and matches stage 221's
  pattern) → update the **paper card** line 11 to cite the `.wl` via
  `\StageFile{mathematica/moving_throat_pde_stage220_dynamic_mixed_port_kernel_phase_lag_no_go_and_resonant_survival_gate_mathematica_audit.wl}` (replacing "none yet"),
  and refresh notes Section 10 to "SymPy and Mathematica audits verify …" with the
  `.wl` added to the supporting-file list (line 592). No script change. (This is a
  paper/notes edit that — per the file-ownership rule — Codex may apply only AFTER
  the user authorizes the direction.)
- (b) The `.wl` is not intended to be cited at this stage → leave the docs as-is and
  record why the engine is present-but-uncited. No edit.

The orchestrator will not invoke Codex on this unit until the user has chosen a
direction. This is a doc-side label correction only; the math is sound and
dual-engine green, so there is no `## F<n>` Codex-applied block.
