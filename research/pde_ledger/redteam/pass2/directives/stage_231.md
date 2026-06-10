---
unit_id: 231
batch: VII.2
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 231

The only finding is a `paper_misalignment`. Codex applies NOTHING on this unit
until the user resolves direction. Do not edit `paper/`, `notes/`, or any script
to "fix" this discrepancy. The orchestrator is holding for user resolution.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_231.tex:11` quote: "\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_sympy_audit.py}.  Mathematica audit: none yet.}"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl` — exists, committed (Jun 2), and PASSES; its saved output ends "All Stage 231 Mathematica audits passed."

## Resolve before fix_loop

The card says "Mathematica audit: none yet," but a passing dual-engine `.wl` now
exists for Stage 231 and independently confirms every deliverable (it is stronger
than the SymPy: it proves the monotonicity/positivity claims via `Resolve[ForAll]`
quantifier elimination). The card understates verification coverage. Which side is correct?

Possible directions (the user picks one):
- (a) Paper is stale → update `stage_231.tex:11` so the Verification line cites the
  `.wl` (e.g. "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage231_continuum_placement_pullback_of_the_selected_branch_dynamic_class_map_mathematica_audit.wl}."), no script change. (This is the expected direction under the dual-engine-required policy, matching the Stage-230 resolution.)
- (b) The `.wl` should not count as the stage's Mathematica audit → leave the card and
  justify why the present `.wl` is not the canonical second engine.

The orchestrator will not invoke Codex on this unit until the user has chosen a
direction. Because the resolution is a paper-side prose edit, the follow-up is a
paper-side edit by Codex (per file-ownership rules) reviewed by Claude — not a
script change.
</content>
