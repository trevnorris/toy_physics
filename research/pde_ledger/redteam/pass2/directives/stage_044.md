---
unit_id: 044
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:52:43Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 044

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `27`; canonical = `044`) in the SymPy source, matching the file's canonical banner
(`STAGE 44`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 044): docstring header + `27.k` subbanner indices. Banner already canonical.

**Required change (label-only — change ONLY the leading numeric token; preserve the rest of each subbanner string verbatim):**
- line 3: `Moving-throat PDE — Stage 27 SymPy audit.` → `Moving-throat PDE — Stage 44 SymPy audit.`
- line 53: subbanner `27.1` → `44.1`
- line 80: subbanner `27.2` → `44.2`
- line 111: subbanner `27.3` → `44.3`
- line 134: subbanner `27.4` → `44.4`
- line 148: subbanner `27.5` → `44.5`

**DO NOT TOUCH (deferred CROSS-refs — dedicated pass owns these):**
- lines 8, 11, 55, 82, 121 (`Stage-24` / `Stage-25` / `Stage 24/25` cross-refs to upstream stages 041/042) — LEAVE.
- the already-canonical `STAGE 44` banners (lines 51, 156) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring + subbanners read `44` / `44.1`–`44.5`; re-run all PASS; strip-the-number diff identical.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py`
- summary: Updated the stale self-stage docstring and subbanner labels from 27/27.k to 44/44.k.
- deviation: none
