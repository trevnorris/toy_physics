---
unit_id: 214
batch: VI.1
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 214

Apply the finding below. After applying, append an `## Applied: F1` block with:
`files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Do NOT edit the
`.py` or `.wl` script logic — the math is correct and verified on both engines.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`

**Issue:** The committed SymPy output is stale relative to the script. Its mtime
(2026-06-02 11:59:40) predates the `.py` mtime (2026-06-03 15:59:11), and its
content carries the pre-renumber banner `STAGE 197` on line 3
(`STAGE 197 — FULL INTERIOR FOUR-COORDINATE SIMPLEX OPTIMIZER`) and line 1209
(`STAGE 197 SYMPY AUDIT COMPLETED SUCCESSFULLY`), whereas the current script
banner (`.py:41`) reads `STAGE 214`. All numerical/symbolic content (degrees,
`54`, `150`, all `=0` residuals, the `924`-sample sweep counts) already matches
the current script; only the stage label is stale.

**Required change:**
Re-run the SymPy script and recommit its captured output so the banner reads
`STAGE 214`. No source edit. Concretely:
1. `python3 scripts/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.py` and capture stdout to
   `scripts/output/moving_throat_pde_stage214_full_interior_four_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`.
2. Confirm the script exits 0 with every in-file check passing (no `AssertionError`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 214`, confirm the
script exits 0, and confirm the regenerated `.txt` line 3 reads
`STAGE 214 — FULL INTERIOR FOUR-COORDINATE SIMPLEX OPTIMIZER` and the final line
reads `STAGE 214 SYMPY AUDIT COMPLETED SUCCESSFULLY`, with all residuals `0` and
the two sweep counts equal to `924`.
