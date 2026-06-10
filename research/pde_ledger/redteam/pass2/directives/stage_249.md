---
unit_id: 249
batch: VIII.1
created_at: 2026-06-10T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 249

The only finding is a `paper_misalignment`. Do nothing — the orchestrator is holding
for user resolution. Do not edit paper.tex, notes/, or scripts to "fix" this unless
a follow-up directive explicitly authorizes a direction.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale Verification line; resolution is paper-side)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_249.tex:4` quote: "\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage249_..._sympy_audit.py}.  Mathematica audit: none yet.}"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage249_conditional_helicity_export_diagnostic_compiler_on_the_dynamic_event_chain_and_aligned_vs_anti_aligned_mixed_sector_closure_mathematica_audit.wl` exists and passes (M1-M5 all PASS; see `mathematica/output/..._mathematica_audit.txt`).

## Resolve before fix_loop

The stage card's Verification line states "Mathematica audit: none yet," but a passing,
independent Mathematica `.wl` audit now exists (added in pass-1). The card text is stale.

Possible directions (the user picks one):
- (a) Update the card's `\stagefield{Verification}` line at `stage_249.tex:4` to name the
  `.wl` file (e.g. "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage249_..._mathematica_audit.wl}.")
  instead of "none yet". No script change. (Recommended — the `.wl` is real, independent,
  and passing.)
- (b) If the `.wl` is intended to be removed/not counted, drop it and leave the card as-is.
  (Not recommended — it is a genuine independent second engine.)

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
Codex makes no edits in the meantime.
