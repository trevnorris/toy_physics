---
unit_id: 050
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T17:49:14Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 050

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the lines named — number-only label corrections, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py:118`

**Issue:** This is stage 050 but two SELF-labels carry the pre-renumber "Stage 33" (33+17=50). Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Labels only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3` docstring first line: `Stage 33 SymPy audit.` → `Stage 50 SymPy audit.`
2. `py:118`: `print("\nAll Stage-33 symbolic checks passed.")` → `print("\nAll Stage-50 symbolic checks passed.")`
3. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

**DO NOT TOUCH (deferred to the dedicated SCRIPT/OUTPUT-band numbering plan — these are CROSS-references to other stages, not self-labels):**
- `py:61-62` comment `Imported from Stage 32's explicit D/N overlap extraction` — "Stage 32" here is a cross-reference to the IMPORTED stage 049, not a self-label. Leave it exactly as-is.
- `py:17` import statement and `py:48` "stage 049 import" — already correct; leave.

Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 050` exits 0; refreshed output closing reads `All Stage-50 symbolic checks passed.`; `grep -nE "Stage 33|Stage-33" <py>` returns nothing; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage050_zeta_threshold_comparison_sympy_audit.txt`
- summary: Updated the two stale Stage 33 self-labels to Stage 50 and refreshed the SymPy saved output.
- deviation: none
