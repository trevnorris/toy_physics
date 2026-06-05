---
unit_id: 051
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T17:49:41Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 051

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the line named — number-only label correction, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous filename self-label; number-only)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py:3`

**Issue:** This is stage 051 but the docstring restates the module's own filename with the pre-renumber stem `stage34` (34+17=51). The real on-disk filename is `...stage051...` (3-digit). Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Label only.

**Required change (number-only; filename form is canonically 3-digit):**
1. `py:3`: `moving_throat_pde_stage34_lowest_twin_criterion_sympy_audit.py` → `moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`
2. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

**DO NOT TOUCH (deferred to the dedicated SCRIPT/OUTPUT-band numbering plan — compound/cross-references, not unambiguous self-labels):**
- `py:20-21` `the Stage 050/034 product law` (compound dual-epoch label) and `the Stage 047 coherent map` (cross-ref).
- `py:126` and `wl:87` `Stage 047/030 coherent forward map` (compound dual-epoch cross-ref).
- `STAGE 51` banners (py:63, py:149) already carry the correct number — leave their format.

Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 051` exits 0; `py:3` reads `...stage051...`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage051_lowest_twin_criterion_sympy_audit.txt`
- summary: Corrected the stale stage34 filename self-label to the canonical stage051 form and refreshed the SymPy output by rerunning the audit.
- deviation: none
