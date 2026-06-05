---
unit_id: 055
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T17:55:15Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 055

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

Edit exactly the line named — number-only label correction, nothing else. No refactors.
After editing, RUN `python3 <path>` and iterate until it exits 0 with all checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-label; number-only, format preserved)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py:3`

**Issue:** This is stage 055 but the docstring carries the pre-renumber "Stage 38" (38+17=55). Committed `.txt` outputs are stale and refresh on re-run. The `.wl` is already canonical. Label only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3`: `Stage 38 SymPy audit: explicit reachability of the non-twin lowest support lane.` → `Stage 55 SymPy audit: explicit reachability of the non-twin lowest support lane.`
2. Re-run SymPy (and Mathematica via the orchestrator) to refresh outputs.

**DO NOT TOUCH (deferred to the dedicated SCRIPT/OUTPUT-band numbering plan — CROSS-reference to another stage, not a self-label):**
- `py:73` `It reaches the Stage-35 threshold iff ...` — "Stage-35" is a cross-reference to the threshold stage 052 (35+17=52), not a self-label. Leave it exactly as-is.
- `py:27` `(Stage 054 softening)` — already-correct cross-ref; leave.

Do NOT pad the already-correct `STAGE 55` banner (py:25). Do not alter any assertion, symbol, or numeric expression.

**Verification:** `redteam exec-sympy 055` exits 0; `py:3` reads `Stage 55 ...`; `py:73` still reads `Stage-35`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage055_explicit_reachability_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage055_explicit_reachability_sympy_audit.txt`
- summary: Corrected the SymPy audit self-label from Stage 38 to Stage 55 and regenerated the matching SymPy output.
- deviation: none
