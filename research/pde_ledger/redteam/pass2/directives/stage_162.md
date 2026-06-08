---
unit_id: 162
batch: IV.6
created_at: 2026-06-08T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-06-08T11:48:27-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 162

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line named.

After editing, RUN the affected script (`python3 scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py`) and confirm it exits 0 with all three `expect_zero` checks reporting `= 0`. (The change is a comment only and cannot alter execution; this is a regression check.)

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output / numbering (comment-only)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py:39`

**Issue:** The source comment misattributes the parent-family formulas to "Stages 99 and 102". This stage's notes attribute them to Stage 119 (parent compensation family, balance law, and `L_W/a=(π/2)√((1+𝔯²)/3)`) and to Stages 115–116 (`γ₀=(1+𝔯²)/9`). The wrong stage numbers are residue of the known +17-era numbering drift (the same garble the Step-0 fix already corrected to `Stage 119` in the notes). The comment is non-load-bearing, so the math and all assertions are unaffected; this is provenance hygiene only.

**Required change:**
Replace the comment at line 39:

- Before:
  `# Exact parent family formulas from Stages 99 and 102.`
- After:
  `# Exact parent family formulas from Stage 119 (with gamma0 from Stages 115-116).`

Do not change any code, any expression, or the Mathematica script (it carries no such comment).

**Verification command:**
The verifier runs `redteam exec-sympy 162` and confirms (a) the script exits 0 with `similarity identity = 0`, `lower-branch differential law = 0`, `positive slope decomposition = 0`, and (b) `grep -n "Stages 99 and 102"` on the `.py` returns nothing while the new comment cites Stage 119 (and 115–116 for gamma0).

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage162_parent_compensation_rigidity_sympy_audit.py`
- summary: Updated the parent-family provenance comment to cite Stage 119 and Stages 115-116 for gamma0.
- deviation: none
