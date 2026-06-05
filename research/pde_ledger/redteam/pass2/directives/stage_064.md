---
unit_id: 064
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
verification_status: pending
needs_user_resolution: false
applied_at: 2026-06-05T19:43:25Z
findings_applied: 1
findings_blocked: 0
---

# Codex directive — unit 064

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edit; preserve format exactly.

⚠️ SCOPE GUARD: Edit ONLY the one SELF-label on line 3. Do NOT touch lines 25, 122, 180 (`.py`) or line 104 (`.wl`) — those reference "Stage-45/46", CROSS-references to stages 062/063 (one is inside an `expect_zero`/`expectZero` assertion label), DEFERRED to the dedicated numbering plan. Do NOT pad the `STAGE 64` banners (lines 49, 176). Do not touch variable names or any other line.

After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py` and confirm exit 0. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-label (numbering)

**Target:** `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py:3`

**Issue:** Stale +17 pre-renumber self-label "Stage 47" in the docstring. Canonical is 064.

**Required change:**
- Line 3: `Moving-Throat PDE — Stage 47 SymPy audit` → `Moving-Throat PDE — Stage 64 SymPy audit`
- Change ONLY `47`→`64`; preserve everything else. Do not pad.

**Verification:** `git diff` shows exactly one changed `.py` line (lines 25/49/122/176/180 and `.wl:104` untouched), identical except `47`→`64`. Script exits 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage064_equilibrium_alignment_sympy_audit.py`
- summary: Updated the stage self-label in the SymPy audit docstring from Stage 47 to Stage 64.
- deviation: none
