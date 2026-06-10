---
unit_id: 230
batch: VII.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 230

The only finding is a `paper_misalignment`. Codex applies NOTHING on this unit
until the user resolves direction. Do not edit `paper/`, `notes/`, or any script
to "fix" this discrepancy. The orchestrator is holding for user resolution.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_230.tex:11` quote: "\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_sympy_audit.py}.  Mathematica audit: none yet.}"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl` — exists, committed (1dfc3fe, batch VII.1), and PASSES; its saved output ends "All Stage 230 Mathematica audits passed."

## Resolve before fix_loop

The card says "Mathematica audit: none yet," but a passing dual-engine `.wl` now
exists for Stage 230. The card understates verification coverage. Which side is correct?

Possible directions (the user picks one):
- (a) Paper is stale → update `stage_230.tex:11` so the Verification line cites the
  `.wl` (e.g. "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage230_selected_branch_classifier_to_dynamic_window_compiler_and_static_first_theorem_mathematica_audit.wl}."), no script change. (This is the expected direction under the dual-engine-required policy.)
- (b) The `.wl` should not count as the stage's Mathematica audit → leave the card and
  justify why the present `.wl` is not the canonical second engine.

The orchestrator will not invoke Codex on this unit until the user has chosen a
direction. Because the resolution is a paper-side prose edit, the follow-up is a
paper-side edit by Codex (per file-ownership rules) reviewed by Claude — not a
script change.
