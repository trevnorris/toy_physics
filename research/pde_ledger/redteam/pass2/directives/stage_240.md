---
unit_id: 240
batch: VII.2
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 240

All findings on this unit are `paper_misalignment` pending user resolution. Codex applies nothing. Do not edit paper.tex, notes/, or scripts until the user chooses a direction.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (card understates verification coverage)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_240.tex:11` quote: `\stagefield{Verification}{SymPy audit: ... Mathematica audit: none yet.}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage240_selected_branch_loading_ratio_from_the_minimal_isotropic_quadrupole_precursor_mathematica_audit.wl` exists, is committed, and passes all six modules (M1-M6) per `mathematica/output/...stage240..._mathematica_audit.txt`.

## Resolve before fix_loop

The card's `\stagefield{Verification}` line says "Mathematica audit: none yet," but a passing `.wl` exists. Should the card be updated to cite the existing Mathematica audit?

Possible directions (the user picks one):
- (a) Card is stale → update line 11 to cite the `.wl` filename (paper-side edit, no script change). RECOMMENDED.
- (b) Intentional → leave as-is (e.g. card text deliberately reserved); no change.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script-side change is required either way; this is a documentation reconciliation only.
