---
unit_id: 037
batch: III.1
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T13:40:32Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 037

Single finding: in-loop fix of the audit-flagged UNAMBIGUOUS stale SELF-labels (pre-renumber
self-number `20`, which identifies THIS stage; canonical = `037`) in the SymPy source, matching
the file's canonical banner (`STAGE 37`). **Label-only** — change ONLY the stage-number token,
preserve all surrounding text/format. Then re-run to confirm exit 0.

After applying, append `## Applied: F1` with `files_changed`, `summary` (one sentence), `deviation` (or "none").

Do NOT touch paper.tex, notes/, or any prose documents. Do NOT touch the Mathematica `.wl` (its
source labels are already canonical). Do NOT change any equation, value, variable name, or assertion.

## F1 — stale_output (self-label) [label-only]

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`

**Issue:** the SymPy source carries stale pre-renumber SELF-labels for THIS stage (canonical 037).
The `e2a4780` numbering commit fixed the main `banner(...)` but missed the docstring filename + header.

**Required change (label-only — change ONLY the indicated token):**
- line 3: `moving_throat_pde_stage20_continuum_kernel_sympy_audit.py` → `moving_throat_pde_stage037_continuum_kernel_sympy_audit.py` (match the real filename, 3-digit)
- line 5: `Stage 20 SymPy audit:` → `Stage 37 SymPy audit:`

**DO NOT TOUCH (deferred cross-refs — the content-keyed dedicated numbering pass owns these):**
- line 6 `the Stage-17/19 reduced branch data` — cross-ref to upstream sources, LEAVE.
- line 224 `Stage-17/19 branch data exactly.` — cross-ref, LEAVE.
- the unprefixed subbanners (`1.`–`5.`) and the already-canonical `STAGE 37` banners (lines 51, 222) — LEAVE.

After editing, RUN `python3 <the .py above>` and confirm it exits 0 with all checks passing. Do NOT run or edit the `.wl`.

**Verification:** the `.py` docstring reads `stage037`/`Stage 37`; re-run shows all checks PASS; no math/value/variable/assertion byte changed (strip-the-number diff is identical).

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage037_continuum_kernel_sympy_audit.py`
- summary: Updated the SymPy audit docstring self-labels from stage20/Stage 20 to stage037/Stage 37.
- deviation: none
