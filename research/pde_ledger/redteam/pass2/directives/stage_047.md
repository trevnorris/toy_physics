---
unit_id: 047
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:54:52Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 047

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `30`; canonical = `047`) in the SymPy source, matching the file's canonical banner
(`STAGE 47`). **Label-only** — change ONLY the stage-number token, preserve surrounding text. Re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary`, `deviation` (or "none").

Do NOT touch paper.tex, notes/, prose docs, or the `.wl` (its labels are already canonical).
Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`

**Issue:** stale pre-renumber SELF-labels for THIS stage (canonical 047): docstring header + the closing pass-line. Banner already canonical.

**Required change (label-only — change ONLY the indicated numeric token, preserve the format):**
- line 3: `Stage 30 SymPy audit.` → `Stage 47 SymPy audit.`
- line 158: `All Stage-30 symbolic checks passed.` → `All Stage-47 symbolic checks passed.`

**DO NOT TOUCH (deferred CROSS-ref — dedicated pass owns it):**
- lines 44–45 `That saturation is established upstream at Stage 28 ...` — cross-ref to upstream stage 045, LEAVE.
- the already-canonical `STAGE 47` banner (line 32) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm exit 0, all checks PASS. Do NOT run or edit the `.wl`.

**Verification:** docstring + closing pass-line read `47`; re-run all PASS; strip-the-number diff identical; the upstream `Stage 28` cross-ref is untouched.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage047_coherent_kernel_map_sympy_audit.py`
- summary: Updated the two stale self-label stage-number tokens from 30 to 47.
- deviation: none
