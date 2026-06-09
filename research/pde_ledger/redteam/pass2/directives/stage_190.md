---
unit_id: 190
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 190

This unit has NO Codex-applied script edits. Both findings are documentation-side:

- **F1 (stale_output)** is informational only — the orchestrator's independent
  `exec-sympy 190` re-run refreshes the saved `.txt` (the only stale content is the
  `STAGE 173 → STAGE 190` banner; all math lines are byte-identical). Codex applies nothing.
- **F2 (paper_misalignment)** is held for USER RESOLUTION. Do not edit `paper/`,
  `notes/`, or scripts to "fix" it. The orchestrator halts until the user picks a direction.

Codex: take no action on this unit until a follow-up directive authorizes a specific edit.

## F2 — paper_misalignment

**Subtype:** paper_missing_script_claim

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_190.tex:11` quote:
  "SymPy audit: \StageFile{scripts/moving_throat_pde_stage190_direct_defect_vs_dressing_split_sympy_audit.py}.  Mathematica audit: none yet."

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl` (283 lines, exists, all checks PASS in `mathematica/output/...mathematica_audit.txt`)

## Resolve before fix_loop

The card's `\stagefield{Verification}` says "Mathematica audit: none yet", but a complete, independent Mathematica audit `.wl` (created during the pass-1 dual-engine retrofit) exists and passes. Which is correct?

Possible directions (the user picks one):
- (a) The `.wl` is the intended second engine → update `stage_190.tex:11` to read
  "Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage190_direct_defect_vs_dressing_split_mathematica_audit.wl}." (paper-side edit, applied by Codex only after this resolution; no script change).
- (b) The `.wl` should be retired → remove it and leave the card as-is (script-side change, user-authorized).
- (c) Defer the card-prose sync to the dedicated post-V.3 tracker pass → no edit now.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
