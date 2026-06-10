---
unit_id: 252
batch: VIII.1
created_at: 2026-06-10T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: true
---

# Codex directive — unit 252

Apply each non-paper_misalignment finding below in order. There are NONE — the
only finding is a `paper_misalignment` held for user resolution. Codex applies
nothing on this unit until the user picks a direction. Do not edit paper.tex,
notes/, or scripts to "fix" this paper_misalignment.

## F1 — paper_misalignment (subtype: paper_missing_script_claim)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_252.tex:4` quote:
  `\stagefield{Verification}{SymPy audit: \StageFile{scripts/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_sympy_audit.py}.  Mathematica audit: none yet.}`

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage252_vacuum_vs_lattice_heat_partition_and_cold_survival_compiler_from_the_microscopic_export_kernel_mathematica_audit.wl` — a complete passing Mathematica audit (sections M1-M9, all PASS in the saved output) now exists; it was created in pass-1.

## Resolve before fix_loop

The stage-252 card's verification line still says "Mathematica audit: none yet," but a working Mathematica audit `.wl` now exists and all of its checks (M1-M9) pass. This is a card-text-lag, not a math defect.

Possible directions (the user picks one):
- (a) Update the card → replace the verification line so it cites the Mathematica audit file (paper-side edit, no script change). This is the expected resolution.
- (b) Leave the card as-is (defer to a later doc-sync sweep) → no edit now; close F1 as a known card-text-lag.

The orchestrator will not invoke Codex on this unit until the user has chosen a direction. No script edits are warranted; the math is sound on both engines.
