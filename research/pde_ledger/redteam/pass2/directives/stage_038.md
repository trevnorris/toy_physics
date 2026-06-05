---
unit_id: 038
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:43:27Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 038

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `21`; canonical = `038`) in the SymPy source, matching the file's canonical banner
(`STAGE 38`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 038); the banner was already fixed, the docstring filename + header were missed. (Both committed `.txt` are also stale — the orchestrator's re-run refreshes them.)

**Required change (label-only — change ONLY the indicated token):**
- line 3: `moving_throat_pde_stage21_dimensionless_continuum_placement_sympy_audit.py` → `moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py`
- line 5: `Stage 21 SymPy audit:` → `Stage 38 SymPy audit:`

**DO NOT TOUCH (deferred CROSS-refs to upstream stage 037 (pre-renumber "20") — dedicated pass owns these):**
- lines 6, 39, 42, 137, 225 — every `Stage-20` / `Stage-20 continuum formulas` is a cross-ref to upstream stage 037, LEAVE.
- the already-canonical `STAGE 38` banners (lines 36, 224) and unprefixed subbanners (`1.`–`4.`) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring reads `stage038`/`Stage 38`; re-run all PASS; strip-the-number diff identical.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage038_dimensionless_continuum_placement_sympy_audit.py`
- summary: Updated the stale self-labels in the SymPy audit docstring from stage21/Stage 21 to stage038/Stage 38.
- deviation: none
