---
unit_id: 056
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T11:58:23-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 056

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the lines named — number-only label corrections, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py:97`

**Issue:** This is stage 056 but two SELF-labels carry the pre-renumber "Stage 39" (39+17=56): the docstring header and the closing print. Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Labels only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3` docstring line: `Moving-Throat PDE — Stage 39 SymPy audit.` → `Moving-Throat PDE — Stage 56 SymPy audit.`
2. `py:97`: `print("\nStage 39 audit passed.")` → `print("\nStage 56 audit passed.")`
3. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

**DO NOT TOUCH (deferred to the dedicated SCRIPT/OUTPUT-band numbering plan — CROSS-reference to another stage, not a self-label):**
- `py:7` `reproduces the Stage-36 boost formula` — "Stage-36" is a cross-reference to stage 053 (36+17=53), not a self-label. Leave it exactly as-is.

Do NOT pad the already-correct `STAGE 56` banner (py:31). Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 056` exits 0; `py:3` reads `Stage 56 ...`, `py:97` reads `Stage 56 audit passed.`; `py:7` still reads `Stage-36`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage056_transport_source_asymmetry_sympy_audit.txt`
- summary: Updated the two stale Stage 39 self-labels to Stage 56 and refreshed the SymPy audit output.
- deviation: none
