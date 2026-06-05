---
unit_id: 042
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:45:00Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 042

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `25`; canonical = `042`) in the SymPy source, matching the file's canonical banner
(`STAGE 42`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 042): docstring header + `25.k` subbanner indices. Banner already canonical.

**Required change (label-only — change ONLY the leading numeric token; preserve the rest of each subbanner string verbatim):**
- line 3: `Moving-throat PDE — Stage 25 SymPy audit.` → `Moving-throat PDE — Stage 42 SymPy audit.`
- line 53: subbanner `25.1` → `42.1` (leave the trailing `... the Stage-24 support loading` text — cross-ref)
- line 75: subbanner `25.2` → `42.2`
- line 103: subbanner `25.3` → `42.3` (leave the trailing `... back to Stage 23` text — cross-ref)
- line 114: subbanner `25.4` → `42.4`
- line 132: subbanner `25.5` → `42.5`
- line 140: subbanner `25.6` → `42.6`

**DO NOT TOUCH (deferred CROSS-refs / variable names — dedicated pass owns these):**
- lines 7, 12, 184 (`Stage-24` / `Stage-23 two-vector law` cross-refs), and the `Stage-24`/`Stage 23` tails on lines 53/103 — LEAVE.
- variable name `F_stage23` (lines 106, 112) — CODE identifier, LEAVE.
- the already-canonical `STAGE 42` banners (lines 51, 180) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring + subbanners read `42` / `42.1`–`42.6`; re-run all PASS; strip-the-number diff identical; no `F_stage23` rename.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage042_rank2_selected_mode_sympy_audit.py`
- summary: Updated the stale self-stage docstring and subbanner numeric prefixes from 25/25.k to 42/42.k.
- deviation: none
