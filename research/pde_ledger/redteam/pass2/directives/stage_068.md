---
unit_id: 068
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:46:15-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 068

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edit; preserve format exactly.

⚠️ SCOPE GUARD: Edit ONLY the filename-style self-label in the docstring (line 4). Do NOT pad or alter the `STAGE 68` ledger banner at line 37 — it is already the correct canonical number; banner padding is reserved for the dedicated plan. Do not touch any other line or variable name.

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py` and confirm exit 0. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-label (numbering)

**Target:** `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py:4`

**Issue:** The docstring's filename-style self-label still reads the stale +17 pre-renumber `stage51`. The on-disk filename is `...stage068_...`. Unambiguous SELF-label; filename-style, so the 3-digit `stage068` form is correct (matches the actual filename).

**Required change:**
- Line 4: `moving_throat_pde_stage51_resonance_thresholds_sympy_audit.py` → `moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
- Change ONLY `stage51`→`stage068` (3-digit, matching the on-disk filename); preserve everything else byte-identical.

**Verification:** `git diff` shows exactly one changed line, identical except `stage51`→`stage068`. Script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage068_resonance_thresholds_sympy_audit.py`
- summary: Updated the docstring self-label from `stage51` to `stage068` to match the on-disk filename.
- deviation: none
