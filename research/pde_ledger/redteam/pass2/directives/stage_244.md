---
unit_id: 244
batch: VIII.1
created_at: 2026-06-10T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 244

The only finding is a `paper_misalignment`. Codex applies nothing until the user
resolves direction. Do not edit paper.tex, notes/, or scripts to "fix" this.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_244.tex:4` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl` — exists; saved output shows all M1–M9 PASS.

## Resolve before fix_loop

The card's `\stagefield{Verification}` says "Mathematica audit: none yet," but pass-1
added an independent `.wl` that passes. Should the card be updated to cite the
Mathematica audit?

Possible directions (the user picks one):
- (a) Card is stale → update line 4 to cite `mathematica/moving_throat_pde_stage244_selected_branch_leakage_and_scalar_photon_work_compiler_mathematica_audit.wl` (paper-side edit, no script change).
- (b) Leave as-is → no change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
