---
unit_id: 224
batch: VII.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 224

This directive contains ONLY a `paper_misalignment` finding requiring user resolution.
Codex does nothing here: do not edit paper.tex, notes/, or scripts. The orchestrator
holds for the user to choose a direction.

## F1 — paper_misalignment

**Subtype:** notes_contradicts_script (paper-card verification-status staleness)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_224.tex:11` quote: "\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_sympy_audit.py}.  Mathematica audit: none yet.}"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage224_pde_branch_packet_compiler_weak_axisymmetric_ceiling_transport_and_first_actual_branch_kill_test_mathematica_audit.wl:1-218` — a complete, independent Mathematica audit (sections M1-M6) exists and all checks PASS in the committed output `mathematica/output/...mathematica_audit.txt:87` ("All Stage 224 Mathematica checks passed.").

## Resolve before fix_loop

The card says "Mathematica audit: none yet", but a complete passing `.wl` audit now exists for this stage. The card understates the dual-engine verification coverage. Which direction is correct?

Possible directions (the user picks one):
- (a) Card is stale → update `paper/stages/stage_224.tex:11` to reference the existing `.wl` (e.g. `Mathematica audit: \StageFile{mathematica/moving_throat_pde_stage224_..._mathematica_audit.wl}.`), no script change. (This is the expected resolution — the `.wl` is genuinely independent and passes.)
- (b) The `.wl` should not count (e.g. policy reason) → leave the card and note why the `.wl` is non-authoritative.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script-side changes are required; the SymPy and Mathematica scripts are both clean and load-bearing.
