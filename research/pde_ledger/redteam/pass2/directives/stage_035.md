---
unit_id: 035
batch: II.1
created_at: 2026-06-04T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-04T23:09:41-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 035

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes beyond what is named. Edit exactly the file:line ranges named.

After editing, RUN the affected script (`python3 <path>`) and confirm it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output (low severity, informational)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py:3` and `:112`

**Issue:** The saved `.txt` transcripts (mtime 2026-05-21) are older than the scripts (mtime 2026-06-03), but the only intervening commit (`e2a4780`) was a label-only banner relabel — no math changed, so the captured numeric/symbolic results are still correct. The remaining defect is purely cosmetic: the sympy module docstring at `.py:3` still reads `Moving-throat PDE Stage 18 SymPy audit.` and the final print at `.py:112` still reads `All Stage 18 checks passed.`, both stale "Stage 18" labels left over from the relabel pass. The Mathematica script's banners were already corrected to "STAGE 035".

**Required change** — match the `.py`'s OWN main banner format (2-digit `STAGE 35.x`); change only the number `18` → `35`:
- `.py:3`: change `Moving-throat PDE Stage 18 SymPy audit.` → `Moving-throat PDE Stage 35 SymPy audit.`
- `.py:112`: change `print("All Stage 18 checks passed.")` → `print("All Stage 35 checks passed.")`

These are label-only edits to comment/print strings; do NOT touch any equation, value, variable, assertion, or `expect_zero` target. If you judge any of these strings load-bearing for a downstream parser, append `## Blocked: F1` with the concern instead of editing.

**Verification:**
After edit, `python3 <path>` exits 0; the regenerated sympy transcript header reads `STAGE 35.1 — ...` and the final line reads `All Stage 35 checks passed.`; every residual line remains `0`. (The output-staleness itself resolves automatically on re-run; the orchestrator re-runs both engines regardless.)

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage035_dimensionless_normalization_locus_sympy_audit.py`
- summary: Updated the stale SymPy stage labels from Stage 18 to Stage 35 in the module docstring and final success print.
- deviation: none
