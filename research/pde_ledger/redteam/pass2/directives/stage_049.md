---
unit_id: 049
batch: III.2
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-05T18:03:17Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 049

Apply the finding below. After applying, append an `## Applied: F1` block with `files_changed`, `summary` (one sentence), `deviation` (or "none").

Edit exactly the lines named — number-only label corrections, nothing else. No refactors, no reformatting.
After editing, RUN both engines (`python3 <path>`, `math -script <path>`) and iterate until they exit 0 with all in-file checks passing. The orchestrator re-runs independently afterward.
Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (unambiguous self-labels; number-only, format preserved)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py:3`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl:93`

**Issue:** This is stage 049 but two live SELF-labels carry the pre-renumber "Stage 32" (known +17 drift; 32+17=49). Committed `.txt` outputs are stale (print the old `STAGE 32` banner) and refresh on re-run. Math is verified correct — labels only.

**Required change (number-only, keep existing 2-digit format):**
1. `py:3` docstring first line: `Stage 32 SymPy audit.` → `Stage 49 SymPy audit.`
2. `wl:93`: `Print["Stage 32 Mathematica audit passed."];` → `Print["Stage 49 Mathematica audit passed."];`
3. Re-run both engines to refresh the committed `.txt` outputs.

Do NOT pad the already-correct `STAGE 49` banners (py:51, wl:36). Do not alter any assertion, helper, formula, or assumption.

**Verification:** `redteam exec-sympy 049` and `redteam exec-mathematica 049` both exit 0; refreshed outputs print `STAGE 49 — ...` and `Stage 49 Mathematica audit passed.`; no result value changes.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.wl`
  - `scripts/output/moving_throat_pde_stage049_dn_overlap_zeta_sympy_audit.txt`
  - `mathematica/output/moving_throat_pde_stage049_dn_overlap_zeta_mathematica_audit.txt`
- summary: Corrected the stale Stage 32 self-labels to Stage 49 and refreshed both saved engine outputs.
- deviation: none
