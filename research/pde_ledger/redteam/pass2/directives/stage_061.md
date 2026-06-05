---
unit_id: 061
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T19:43:38Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 061

Label-only self-label canonicalization (pass-2 numbering, in-loop Reading-2 scope). NUMBER-only edit; preserve the existing format exactly. Do NOT touch any other line, any cross-reference to another stage, any variable name, or any banner (the banner is already canonical). After editing, RUN `timeout 600 python3 scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py` and confirm exit 0 with all checks passing. No `.wl` edit.

After applying, append an `## Applied: F1` block with `files_changed`, `summary`, `deviation` (or "none").

## F1 — stale self-label (numbering)

**Target:** `scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py:3`

**Issue:** The module docstring carries the stale +17 pre-renumber self-label "Stage 44". Canonical stage number is 061; the filename and banner are already canonical. Unambiguous SELF-label.

**Required change:**
- Line 3: `Stage 44 SymPy audit — microscopic gain thresholds and operator phase diagram.` → `Stage 61 SymPy audit — microscopic gain thresholds and operator phase diagram.`
- Change ONLY `44`→`61`. Keep "Stage" and the rest byte-identical. Do not pad to 3 digits.

**Verification:** `git diff` shows exactly one changed line, identical except `44`→`61`. Script exits 0. The orchestrator re-runs and refreshes the committed `.txt` (banner already prints `STAGE 61`/`STAGE 061`).

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage061_microscopic_gain_thresholds_sympy_audit.py`
- summary: Updated the self-label in the module docstring from Stage 44 to Stage 61.
- deviation: none
