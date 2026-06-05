---
unit_id: 058
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T12:00:10-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 058

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the lines named — number-only label corrections, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py:229`

**Issue:** This is stage 058 but two SELF-labels carry the pre-renumber "Stage 41" (41+17=58): the docstring header and the closing print. Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Labels only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3` docstring line: `Moving-Throat PDE — Stage 41 SymPy audit.` → `Moving-Throat PDE — Stage 58 SymPy audit.`
2. `py:229`: `print("\nStage 41 audit passed.")` → `print("\nStage 58 audit passed.")`
3. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

Do NOT pad the already-correct `STAGE 58` banner (py:33). Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 058` exits 0; `py:3` reads `Stage 58 ...`, `py:229` reads `Stage 58 audit passed.`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage058_coupled_support_source_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage058_coupled_support_source_sympy_audit.txt`
- summary: Corrected the two stale SymPy self-labels from Stage 41 to Stage 58 and refreshed the saved SymPy transcript.
- deviation: none
