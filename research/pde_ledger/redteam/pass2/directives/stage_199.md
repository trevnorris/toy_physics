---
unit_id: 199
batch: V.3
created_at: 2026-06-09T18:51:04Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 199

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

There are **no code edits** required in this directive. Both findings are informational. Do NOT touch paper.tex, notes/, or any prose document.

## F1 — stale_output (no code change; refresh only)

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.txt`

**Issue:** The committed SymPy transcript is older than the `.py` (output mtime 2026-06-01 12:17:54 < script mtime 2026-06-03 15:59:11) and carries the stale pre-renumber banner `STAGE 182 — …` (line 3) / `STAGE 182 LEDGER` (line 420), whereas the current `.py:40` banner correctly reads `STAGE 199`. The body math is otherwise consistent with the current script — every residual is `0` and the printed symbolic forms match the current code. This is informational/cosmetic only.

**Required change:** No source edit. Re-run the SymPy script to regenerate the saved output so its banner reads `STAGE 199`:
`python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage199_pairwise_orbit_transport_law_sympy_audit.py`
and confirm it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward; if the orchestrator's refresh already covers this, no action is needed.

**Verification command:** After refresh, `scripts/output/...stage199...sympy_audit.txt` line 3 reads `STAGE 199 — EXACT PAIRWISE ORBIT-TRANSPORT LAW…`, the ledger banner reads `STAGE 199 LEDGER`, and all residuals remain `0`. The Mathematica transcript is already fresh and correctly labeled `STAGE 199`; no Mathematica action needed.

## F2 — paper-side card-text lag (NO Codex action)

**Target (paper-side, do NOT edit):** `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_199.tex:11`

**Issue:** The card's `\stagefield{Verification}` says "Mathematica audit: none yet," but a passing independent `.wl` audit now exists. This is a prose currency lag, not a math discrepancy (script and paper agree). Codex does not edit paper/. Recorded for the orchestrator/user only; no directive action.
