---
unit_id: 241
batch: VII.2
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 241

This unit's only finding is a `paper_misalignment`. Codex applies NOTHING. The
orchestrator holds for user resolution. Do not edit paper.tex, notes/, or scripts.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (card understates verification coverage)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_241.tex:11` quote:
  "SymPy audit: \StageFile{scripts/moving_throat_pde_stage241_..._sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage241_exact_primitive_ranking_on_the_selected_twin_support_branch_mathematica_audit.wl` exists, is independently derived (Solve-based, not a port), and its committed output passes all M1-M8 checks.

## Resolve before fix_loop

The Stage 241 card says "Mathematica audit: none yet," but a passing, genuinely
independent Mathematica audit `.wl` (with committed output) already exists. The
card understates verification coverage.

Possible directions (the user picks one):
- (a) Card is stale -> update `stage_241.tex:11` to cite the `.wl` audit and its
  output (paper-side edit, no script change). RECOMMENDED.
- (b) The `.wl` is not intended to count as the stage's Mathematica audit -> leave
  the card and note why the `.wl` is provisional.

The orchestrator will not invoke Codex on this unit until the user has chosen a
direction. No script-side change is required either way.
