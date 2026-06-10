---
unit_id: 246
batch: VIII.1
created_at: 2026-06-10T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 246

The only finding is a paper-side `paper_misalignment`. Do nothing — the orchestrator is holding for user resolution. Do not edit paper.tex, notes/, or scripts until the user chooses a direction in a follow-up directive.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale Verification field)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_246.tex:4` quote: "SymPy audit: \StageFile{...stage246..._sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage246_compensated_multimode_source_compiler_beyond_positive_family1_mathematica_audit.wl` — a full independent second-engine audit (M1–M9) now exists and passes (output committed 2026-06-03, all PASS).

## Resolve before fix_loop

The card's `\stagefield{Verification}` line still says "Mathematica audit: none yet," but a complete, independent Mathematica audit was added in pass-1 and passes. The card is stale.

Possible directions (the user picks one):
- (a) Update the card's Verification line to cite the `.wl` file (paper-side prose edit, user-authorized; Codex/Claude per file-ownership policy). No script change.
- (b) Leave the card as-is intentionally (e.g., the `.wl` is provisional) — then close as wontfix.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script edits are warranted; the math is clean in both engines.
