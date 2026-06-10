---
unit_id: 250
batch: VIII.1
created_at: 2026-06-10T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 250

The only finding is a `paper_misalignment` pending user resolution. Codex applies nothing on this unit until the user chooses a direction. Do NOT edit paper.tex, notes/, or scripts to "fix" this.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_250.tex:4` quote: "SymPy audit: \StageFile{...sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage250_crossing_vs_collapse_goldilocks_window_compiler_from_the_stage248_event_chain_and_relaxed_wall_timescale_closure_mathematica_audit.wl` — a complete, independent, passing Mathematica audit (M1–M7 all PASS in the committed saved output `mathematica/output/...mathematica_audit.txt`).

## Resolve before fix_loop

The stage card's `\stagefield{Verification}` line says "Mathematica audit: none yet," but a full passing dual-engine `.wl` exists and verifies the load-bearing claims (including the global `Resolve[ForAll]` monotonicity / one-sided-window certificates). The card under-reports verification coverage.

Possible directions (the user picks one):
- (a) Update the card line 4 to cite the existing `.wl` (e.g., "Mathematica audit: \StageFile{mathematica/...mathematica_audit.wl}.") — paper-side text edit, no script change. (Recommended: the `.wl` is real and passing.)
- (b) Leave the card text as a known doc-lag to be swept by the dedicated numbering/coverage pass — no edit now.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No math/script change is implied either way.
