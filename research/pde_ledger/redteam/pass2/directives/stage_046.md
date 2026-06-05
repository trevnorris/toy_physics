---
unit_id: 046
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T07:52:22-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 046

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `29`; canonical = `046`) in the SymPy source, matching the file's canonical banner
(`STAGE 46`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 046): docstring header + the closing pass-line. Banner already canonical.

**Required change (label-only — change ONLY the indicated numeric token, preserve the format):**
- line 3: `Stage 29 SymPy audit.` → `Stage 46 SymPy audit.`
- line 193: `All Stage-29 symbolic checks passed.` → `All Stage-46 symbolic checks passed.`

**DO NOT TOUCH:**
- the already-canonical `STAGE 46` banner (line 39) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring + closing pass-line read `46`; re-run all PASS; strip-the-number diff identical.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage046_tracking_branch_bounds_sympy_audit.py`
- summary: Updated the stale self-label stage number in the docstring header and closing pass-line from 29 to 46.
- deviation: none
