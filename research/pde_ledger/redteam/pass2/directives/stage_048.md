---
unit_id: 048
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T07:54:39-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 048

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-label (pre-renumber
self-number `31`; canonical = `048`) in the SymPy source, matching the file's canonical banner
(`STAGE 48`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`

**Issue:** stale pre-renumber SELF-label for THIS stage (canonical 048): the docstring header. Banner already canonical; the closing pass-line (`All Stage-048 ... passed.`) is already canonical. (Both committed `.txt` are also stale — the orchestrator's re-run refreshes them.)

**Required change (label-only — change ONLY the indicated token):**
- line 3: `Stage 31 SymPy audit.` → `Stage 48 SymPy audit.`

**DO NOT TOUCH:**
- line 127 `All Stage-048 symbolic checks passed.` — already canonical, LEAVE.
- the already-canonical `STAGE 48` banner (line 31) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring reads `Stage 48`; re-run all PASS; strip-the-number diff identical.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage048_support_compensation_sympy_audit.py`
- summary: Updated the SymPy audit docstring self-label from Stage 31 to Stage 48.
- deviation: none
