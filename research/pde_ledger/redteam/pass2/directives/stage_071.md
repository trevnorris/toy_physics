---
unit_id: 071
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T19:46:19Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 071

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edits; preserve each label's existing format exactly. Do NOT pad the `STAGE 71` ledger banner (correct number, reserved for the dedicated plan), and do not touch any cross-reference or variable name.

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py` and confirm exit 0. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-labels (numbering)

**Target:** `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py:3` and `:5`

**Issue:** The docstring carries two stale +17 pre-renumber SELF-labels: a filename-style label (line 3, `stage54`) and a prose label (line 5, `Stage 54`). Canonical is 071. The filename-style label keeps the 3-digit `stage071` form (matches the on-disk filename); the prose label keeps its 2-digit `Stage 71` form.

**Required change:**
- Line 3: `moving_throat_pde_stage54_tanh_wall_branch_sympy_audit.py` → `moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py` (`stage54`→`stage071`, 3-digit)
- Line 5: `SymPy audit for Stage 54:` → `SymPy audit for Stage 71:` (`54`→`71`, 2-digit)
- Change ONLY the number tokens; preserve everything else byte-identical.

**Verification:** `git diff` shows exactly two changed lines, identical except the number tokens. Script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage071_tanh_wall_branch_sympy_audit.py`
- summary: Updated the two stale Stage 54 self-label numbers to canonical stage071/Stage 71.
- deviation: none
