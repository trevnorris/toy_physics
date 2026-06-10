---
unit_id: 226
batch: VII.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 226

The only finding is a `paper_misalignment` (F1) held for user resolution. Codex applies NOTHING on this unit until the user chooses a direction. Do not edit `paper/`, `notes/`, or scripts to "fix" this.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale verification pointer — card understates what is verified)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_226.tex:11` quote: `\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_sympy_audit.py}.  Mathematica audit: none yet.}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl` (full file) — a complete second-engine audit that passes; committed output `mathematica/output/...stage226...mathematica_audit.txt:88` quote: `All Stage 226 Mathematica checks passed.`

## Resolve before fix_loop

The stage 226 card says "Mathematica audit: none yet," but a full Mathematica `.wl` second-engine audit now exists and passes (engine values agree with SymPy to ~20 digits). The card's verification pointer is stale. Which side is corrected?

Possible directions (the user picks one):
- (a) Card is stale (expected) → update `paper/stages/stage_226.tex:11` to cite the `.wl` at `mathematica/moving_throat_pde_stage226_strict_5pn_even_gate_package_surviving_mixed_corridor_and_pure_transfer_subcorridor_mathematica_audit.wl`. No script change. (This is a paper-side edit; per file-ownership policy Codex applies paper/notes edits only after the user authorizes a direction.)
- (b) The `.wl` should not be advertised yet → leave the card as-is; no edit. (Math is unaffected either way.)

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
