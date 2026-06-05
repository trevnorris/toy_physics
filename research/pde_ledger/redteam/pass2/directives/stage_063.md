---
unit_id: 063
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:43:30-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 063

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edits; preserve format exactly.

⚠️ SCOPE GUARD: Edit ONLY the two SELF-labels named below. Do NOT touch line 76 (`# Insert Stage-44 threshold formulas`) — that is a CROSS-reference to upstream stage 061 and is DEFERRED to the dedicated numbering plan. Do NOT pad the `STAGE 63` banner at line 31 (correct number; padding is reserved for the dedicated plan). Do not touch any other line or variable name.

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py` and confirm exit 0. No `.wl` edit (the `.wl` self-labels are already canonical).

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-labels (numbering)

**Target:** `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py:3` and `:124`

**Issue:** Stale +17 pre-renumber self-labels "Stage 46" in the docstring and closing summary print. Canonical is 063.

**Required change:**
- Line 3: `Stage 46 SymPy audit.` → `Stage 63 SymPy audit.`
- Line 124: `print("\nAll Stage 46 symbolic checks passed.")` → `print("\nAll Stage 63 symbolic checks passed.")`
- Change ONLY `46`→`63` on each line; preserve everything else. Do not pad.

**Verification:** `git diff` shows exactly two changed lines (lines 31 and 76 untouched), each identical except `46`→`63`. Script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage063_parent_thresholds_sympy_audit.py`
- summary: Updated the two stale Stage 46 self-labels to Stage 63.
- deviation: none
