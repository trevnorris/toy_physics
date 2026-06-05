---
unit_id: 057
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T11:57:42-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 057

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the lines named — number-only label corrections, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py:134`

**Issue:** This is stage 057 but two SELF-labels carry the pre-renumber "Stage 40" (40+17=57): the docstring header and the closing print. Committed `.txt` outputs are stale and refresh on re-run. The `.wl` self-labels are already canonical. Labels only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3` docstring line: `Moving-Throat PDE — Stage 40 SymPy audit.` → `Moving-Throat PDE — Stage 57 SymPy audit.`
2. `py:134`: `print("\nStage 40 audit passed.")` → `print("\nStage 57 audit passed.")`
3. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

`py:83-84` reference `Stage 056` (carry-forward cross-ref) — already correct, leave it. Do NOT pad the already-correct `STAGE 57` banner (py:33). Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 057` exits 0; `py:3` reads `Stage 57 ...`, `py:134` reads `Stage 57 audit passed.`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage057_physical_parameter_map_sympy_audit.txt`
- summary: Updated the two stale SymPy self-labels from Stage 40 to Stage 57 and refreshed the saved SymPy transcript.
- deviation: none
