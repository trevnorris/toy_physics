---
unit_id: 066
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T19:46:30Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 066

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edit; preserve format exactly.

⚠️ SCOPE GUARD: Edit ONLY the one SELF-label on line 3. Do NOT touch line 14 (`2. Then the Stage-48 thin-wall thresholds become`) or line 59 (`# Stage-48 thresholds`) — those are CROSS-references to upstream stage 065, DEFERRED to the dedicated numbering plan. Do NOT pad the `STAGE 66` banner. Do not touch any other line or variable name.

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py` and confirm exit 0. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-label (numbering)

**Target:** `scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py:3`

**Issue:** Stale +17 pre-renumber self-label "Stage 49" in the docstring. Canonical is 066.

**Required change:**
- Line 3: `Moving-Throat PDE — Stage 49 SymPy audit` → `Moving-Throat PDE — Stage 66 SymPy audit`
- Change ONLY `49`→`66`; preserve everything else. Do not pad.

**Verification:** `git diff` shows exactly one changed line (lines 14/59 untouched), identical except `49`→`66`. Script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage066_wall_figure_of_merit_sympy_audit.py`
- summary: Changed the SymPy audit docstring self-label from Stage 49 to Stage 66.
- deviation: none
