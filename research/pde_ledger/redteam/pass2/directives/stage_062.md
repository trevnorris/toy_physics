---
unit_id: 062
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T19:44:05Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 062

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edits; preserve format exactly. Do NOT touch any other line, cross-reference, variable name, or banner. After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py` and confirm exit 0. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-labels (numbering)

**Target:** `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:3` and `:112`

**Issue:** Stale +17 pre-renumber self-labels "Stage 45" in the docstring and the closing summary print. Canonical is 062. Both are unambiguous SELF-labels (the closing print writes into the committed transcript).

**Required change:**
- Line 3: `Stage 45 SymPy audit.` → `Stage 62 SymPy audit.`
- Line 112: `print("\nAll Stage 45 symbolic checks passed.")` → `print("\nAll Stage 62 symbolic checks passed.")`
- Change ONLY `45`→`62` on each line; preserve everything else byte-identical. Do not pad.

**Verification:** `git diff` shows exactly two changed lines, each identical except `45`→`62`. Script exits 0. Refreshed transcript closing line reads "All Stage 62 symbolic checks passed."

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py`
  - `scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt`
- summary: Updated the Stage 45 self-labels to Stage 62 in the SymPy audit and refreshed its saved transcript.
- deviation: Regenerated the saved transcript so its closing line reflects the updated print.
