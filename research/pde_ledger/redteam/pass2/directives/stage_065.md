---
unit_id: 065
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T19:46:42Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 065

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edit; preserve format exactly.

⚠️ SCOPE GUARD: Edit ONLY the one SELF-label on line 3. Do NOT touch line 22 (`5. Inserting the Stage-44 thresholds ...`) — that is a CROSS-reference to upstream stage 061, DEFERRED to the dedicated numbering plan. Do NOT pad the `STAGE 65` banner. Do not touch any other line or variable name.

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py` and confirm exit 0. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-label (numbering)

**Target:** `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py:3`

**Issue:** Stale +17 pre-renumber self-label "Stage 48" in the docstring. Canonical is 065.

**Required change:**
- Line 3: `Moving-Throat PDE — Stage 48 SymPy audit` → `Moving-Throat PDE — Stage 65 SymPy audit`
- Change ONLY `48`→`65`; preserve everything else. Do not pad.

**Verification:** `git diff` shows exactly one changed line (line 22 untouched), identical except `48`→`65`. Script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage065_thin_wall_confinement_sympy_audit.py`
- summary: Updated the stale Stage 48 self-label in the script docstring to Stage 65.
- deviation: none
