---
unit_id: 235
batch: VII.2
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 235

The only finding is a `paper_misalignment` pending user resolution. Codex applies
NOTHING on this unit until the user chooses a direction. Do not edit
`paper/stages/stage_235.tex`, the notes file, or the scripts.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (stale verification-coverage statement)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_235.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage235_..._sympy_audit.py}.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage235_rigid_mouth_packet_projectors_static_blind_dressing_line_and_codimension_two_orbit_lock_point_sympy_audit.md:395` quote: "## 8. SymPy-backed status" (lists only the `.py` as supporting file)

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage235_..._mathematica_audit.wl` exists, runs clean, fresh output ends "All Stage 235 Mathematica checks passed." (dual-engine coverage now present)

## Resolve before fix_loop

The card line 11 says the Mathematica audit does not exist ("none yet") and notes §8 is
titled "SymPy-backed status," but a verified Mathematica audit `.wl` (with fresh saved
output, all checks passing) is in fact present. The paper-side coverage statement is stale.
Should the card and notes be updated to record the now-present second engine?

Possible directions (the user picks one):
- (a) Documentation should reflect reality (recommended) → update card line 11 to cite the
  `.wl` ("Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage235_..._mathematica_audit.wl}.")
  and update notes §8 to reference the `.wl` alongside the `.py`. No script change.
- (b) The `.wl` was added in error / should not be advertised → leave the card/notes as-is and
  flag the `.wl` for review. (Conflicts with the standing dual-engine-required policy; unlikely.)
- (c) Defer to the broader numbering/coverage doc-sync pass → no edit now, log only.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction.
Per project policy, any resulting card/notes edit is a paper-side change Codex applies only
under a follow-up directive that explicitly authorizes it.
