---
unit_id: 041
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T14:57:23Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 041

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `24`; canonical = `041`) in the SymPy source, matching the file's canonical banner
(`STAGE 41`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 041): docstring header + `24.k` subbanner indices. Banner already canonical.

**Required change (label-only — change ONLY the leading numeric token; preserve the rest of each subbanner string verbatim):**
- line 3: `Moving-throat PDE — Stage 24 SymPy audit.` → `Moving-throat PDE — Stage 41 SymPy audit.`
- line 53: subbanner `24.1` → `41.1`
- line 84: subbanner `24.2` → `41.2`
- line 96: subbanner `24.3` → `41.3`
- line 104: subbanner `24.4` → `41.4`

**DO NOT TOUCH (deferred CROSS-refs — dedicated pass owns these):**
- lines 13, 151 (`Stage-23 one-direction geometry` cross-refs to upstream stage 040) — LEAVE.
- the already-canonical `STAGE 41` banners (lines 51, 142) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring + subbanners read `41` / `41.1`–`41.4`; re-run all PASS; strip-the-number diff identical.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`
- summary: Updated only the stale self-label stage tokens in the SymPy audit docstring and four subbanner strings from 24/24.k to 41/41.k.
- deviation: none

## F2 — stale_output (self-label, missed comment) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py:107`

**Issue:** an inline comment references THIS stage's own subsection by the stale pre-renumber number: `from section 24.1`. Section `24.1` is this stage's own section (its subbanner is now `41.1`). Unambiguous SELF-label the first pass missed (a comment, not a print — does not reach the transcript — but should be canonical for internal consistency).

**Required change (label-only — change ONLY the numeric token):**
- line 107: `from section 24.1` → `from section 41.1`

**DO NOT TOUCH:** any other line; `n_src`/`n_expected` are variable names — LEAVE.

After editing, RUN `python3 <the .py>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage041_rank2_support_sympy_audit.py`
- summary: Updated the stale self-section comment reference from section 24.1 to section 41.1.
- deviation: none
