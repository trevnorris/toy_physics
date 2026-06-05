---
unit_id: 053
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T17:52:36Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 053

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the line named — number-only label correction, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-label; number-only, format preserved)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py:3`

**Issue:** This is stage 053 but the docstring carries the pre-renumber "Stage 36" (36+17=53). Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Label only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3`: `Stage 36 SymPy audit: exact overlap-boost window for the lowest support lane.` → `Stage 53 SymPy audit: exact overlap-boost window for the lowest support lane.`
2. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

Do NOT pad the already-correct `STAGE 53` banner (py:25) — leave its 2-digit format. Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 053` exits 0; `py:3` reads `Stage 53 ...`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage053_overlap_boost_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage053_overlap_boost_sympy_audit.txt`
- summary: Corrected the SymPy audit docstring from Stage 36 to Stage 53 and refreshed the SymPy output transcript.
- deviation: none
