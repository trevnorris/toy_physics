---
unit_id: 219
batch: VII.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 219

The only finding is a `paper_misalignment` pending user resolution. Codex applies
NOTHING for this unit. Do not edit `paper/`, `notes/`, or any script to "fix" this
until the user has explicitly chosen a direction in a follow-up directive. The
orchestrator holds before invoking Codex.

## F1 — paper_misalignment

**Subtype:** paper_missing_script_claim (engine-coverage statement, not a value)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_219.tex:11` quote: "SymPy audit: \StageFile{scripts/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.py}.  Mathematica audit: none yet."
- `/var/projects/toy_physics/research/pde_ledger/notes/stages/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_sympy_audit.md:462-486` (section 9 "Script-backed status" references only the SymPy file).

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage219_one_port_mixed_bundle_static_kernel_and_square_law_suppression_test_mathematica_audit.wl` — a Mathematica audit now exists (added in VII.1 forward pass, commit `1dfc3fe`) and passes all of M1–M7 (`mathematica/output/...stage219..._mathematica_audit.txt`).

## Resolve before fix_loop

The card's `\stagefield{Verification}` field (and notes §9) say "Mathematica audit: none yet," but a passing dual-engine `.wl` now exists. The card understates the actual verification coverage. Which side should change?

Possible directions (the user picks one):
- (a) Paper is stale → update `stage_219.tex:11` so the `Verification` field references the existing `.wl` (and optionally add the `.wl` to notes §9 / §3 "Mathematica audit"). Paper-side edit; no script change. (Most likely correct: the `.wl` is genuine, independent, and passing.)
- (b) The `.wl` should not count (e.g., a policy reason to keep this stage single-engine) → leave the card and remove/retire the `.wl`. (Unlikely; the dual-engine policy requires the `.wl` where Mathematica can verify, and it can here.)
- (c) Defer as a known doc-sync lag to be batched with other VII.1 card `Verification`-field updates.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No mathematical defect exists in either script; all seven deliverables reconcile and both engines pass.
