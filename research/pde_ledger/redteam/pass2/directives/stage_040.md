---
unit_id: 040
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-06-05T14:57:32Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 040

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `23`; canonical = `040`) in the SymPy source, matching the file's canonical banner
(`STAGE 40`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 040): docstring header + the `23.k` subbanner section indices. Banner already canonical.

**Required change (label-only — change ONLY the leading numeric token; preserve the rest of each subbanner string verbatim):**
- line 3: `Moving-throat PDE — Stage 23 SymPy audit.` → `Moving-throat PDE — Stage 40 SymPy audit.`
- line 55: subbanner `23.1` → `40.1`
- line 81: subbanner `23.2` → `40.2`
- line 101: subbanner `23.3` → `40.3` (leave the trailing `... of Stage 22` text unchanged — that is a cross-ref)
- line 124: subbanner `23.4` → `40.4`

**DO NOT TOUCH (deferred CROSS-refs / variable names — dedicated pass owns these):**
- lines 10, 12, 14, 152, 160 (`Stage-18` / `Stage 22` / `Stage-18/19` cross-refs), and the `... of Stage 22` tail on line 101 — LEAVE.
- variable names `F_stage18`, `G_stage19` and their provenance comments (lines 113–141, which cite the canonical `stage035`/`stage036` paths) — CODE identifiers / cross-refs, LEAVE.
- the already-canonical `STAGE 40` banners (lines 53, 151) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring + subbanners read `40` / `40.1`–`40.4`; re-run all PASS; strip-the-number diff identical; no `F_stage18`/`G_stage19` rename.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
- summary: Updated only the stale self-label tokens from Stage 23/23.k to Stage 40/40.k in the SymPy audit source.
- deviation: none

## F2 — stale_output (self-label, missed comment) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py:136`

**Issue:** an inline comment references THIS stage's own subsection by the stale pre-renumber number: `(section 23.2)`. Section `23.2` is this stage's own section (its subbanner is now `40.2`). This is an unambiguous SELF-label the first pass missed (it is a comment, not a subbanner/print, so it does not reach the transcript — but it should still be canonical for internal consistency).

**Required change (label-only — change ONLY the numeric token):**
- line 136: `(section 23.2)` → `(section 40.2)`

**DO NOT TOUCH:** any other line; `F_general`/`G_general` are variable names — LEAVE.

After editing, RUN `python3 <the .py>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

## Applied: F2

- files_changed:
  - `scripts/moving_throat_pde_stage040_generalized_selected_branch_sympy_audit.py`
- summary: Updated the missed inline self-section reference from `23.2` to `40.2`.
- deviation: none
