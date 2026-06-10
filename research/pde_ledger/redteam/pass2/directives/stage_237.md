---
unit_id: 237
batch: VII.2
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 237

This unit has exactly one finding, and it is a `paper_misalignment` requiring user resolution. Codex applies NOTHING until the user chooses a direction. Do not edit paper.tex, notes/, or scripts.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale Verification status)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_237.tex:11` quote: "SymPy audit: \StageFile{scripts/...stage237..._sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage237_actual_branch_dressing_compiler_finite_static_blind_curve_and_support_blind_post_static_orbit_lock_theorem_mathematica_audit.wl` — a complete, independent, all-passing Mathematica audit (304 lines; saved output reports every M1-M7 check PASS, output mtime 2026-06-03 08:16).

## Resolve before fix_loop

The card's `\stagefield{Verification}` says "Mathematica audit: none yet," but a full, independent, passing Mathematica `.wl` audit exists for this stage. The card understates the dual-engine coverage of an `\StatusExactClosure` stage. Which way should this be reconciled?

Possible directions (the user picks one):
- (a) The `.wl` is canonical and intended → update the card line 11 to name the existing `mathematica/...stage237..._mathematica_audit.wl` (paper-side edit, no script change). This is the expected resolution.
- (b) The `.wl` should not count (e.g. it is provisional / to be removed) → leave the card as-is and the orchestrator notes the `.wl` as out-of-band.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script-side change is warranted — the math (both engines) is sound and fully aligned.
</content>
