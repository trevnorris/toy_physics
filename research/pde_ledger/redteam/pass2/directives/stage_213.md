---
unit_id: 213
batch: VI.1
created_at: 2026-06-09T00:00:00-06:00
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 213

There is no Codex-applicable script-side change in this unit. Do nothing until the user resolves F1.
The orchestrator is holding for user resolution on the paper_misalignment below; F2 is an
informational stale-output flag the verifier's independent re-run refreshes automatically.

Do NOT edit paper.tex, notes/, or scripts to "fix" F1 unless the user explicitly chooses a
direction in a follow-up directive.

## F1 — paper_misalignment

**Subtype:** paper-prose stale relative to present script (verification-coverage statement)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_213.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl` — present, runs clean (output dated 2026-06-02), independently derives the load-bearing objects (M3 `Solve` for the gradient ray, M5 `Maximize` for the leverage bound, M7 `Coefficient`-extraction for the discriminant coefficients).

## Resolve before fix_loop

The card's `\stagefield{Verification}` line says "Mathematica audit: none yet", but a substantive, independently-derived Mathematica audit `.wl` now exists and passes. The card prose understates verification coverage.

Possible directions (the user picks one):
- (a) Card is stale → update line 11 to cite the present `.wl`, e.g. "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_mathematica_audit.wl}." (paper-side edit, no script change). This is the expected resolution given the first-pass dual-engine retrofit.
- (b) The `.wl` is considered not-yet-blessed for this card → leave the card and track separately.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.

## F2 — stale_output (informational, no Codex action)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage213_four_coordinate_mixed_simplex_and_support_cardinality_4_gate_sympy_audit.txt`

**Issue:** The saved SymPy transcript (mtime 2026-05-11) predates the current `.py` (mtime 2026-06-03) and carries the pre-renumber banner "STAGE 196" (the +17 offset; 196+17=213). The current `.py` prints "STAGE 213". Mathematical content of the checks is unchanged (all reduce to 0/True).

**Required change:** None by Codex. The orchestrator's independent exec re-run regenerates the transcript with the correct "STAGE 213" banner.
