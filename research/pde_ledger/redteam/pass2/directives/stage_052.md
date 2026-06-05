---
unit_id: 052
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T17:53:13Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 052

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the line named — number-only label correction, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous filename self-label; number-only)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py:4`

**Issue:** This is stage 052 but the docstring restates the module's own filename with the pre-renumber stem `stage35` (35+17=52). The real on-disk filename is `...stage052...` (3-digit). Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Label only.

**Required change (number-only; filename form is canonically 3-digit):**
1. `py:4`: `moving_throat_pde_stage35_nontwin_asymmetry_threshold_sympy_audit.py` → `moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`
2. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

Do NOT pad the already-correct `STAGE 52` banners (py:42, py:124). Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 052` exits 0; `py:4` reads `...stage052...`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage052_nontwin_asymmetry_threshold_sympy_audit.txt`
- summary: Corrected the SymPy audit docstring filename label to the canonical stage052 stem and refreshed the saved SymPy output.
- deviation: none
