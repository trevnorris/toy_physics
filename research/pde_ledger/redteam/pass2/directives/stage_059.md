---
unit_id: 059
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T18:00:17Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 059

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the lines named — number-only label corrections, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py:124`

**Issue:** This is stage 059 but two SELF-labels carry the pre-renumber "Stage 42" (42+17=59): the docstring header and the closing print. Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Labels only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3` docstring line: `Moving-Throat PDE — Stage 42 SymPy audit.` → `Moving-Throat PDE — Stage 59 SymPy audit.`
2. `py:124`: `print("\nStage 42 audit passed.")` → `print("\nStage 59 audit passed.")`
3. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

**DO NOT TOUCH (deferred to the dedicated SCRIPT/OUTPUT-band numbering plan — CROSS-references to other stages, not self-labels):**
- `py:6` `brackets on the Stage-41 branch interval` — cross-ref to stage 058 (41+17=58).
- `py:9` and `py:75` `the Stage-39 Omega_Pe series` — cross-ref to stage 056 (39+17=56).
Leave all three exactly as-is.

Do NOT pad the already-correct `STAGE 59` banner (py:45). Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 059` exits 0; `py:3` reads `Stage 59 ...`, `py:124` reads `Stage 59 audit passed.`; `py:6/9/75` still read `Stage-41`/`Stage-39`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage059_operator_branch_residual_bounds_sympy_audit.txt`
- summary: Corrected the two stale Stage 42 self-labels to Stage 59 and regenerated the SymPy audit output.
- deviation: none
