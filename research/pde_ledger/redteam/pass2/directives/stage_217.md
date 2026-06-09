---
unit_id: 217
batch: VI.1
created_at: 2026-06-09T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 217

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing, RUN the affected script (`python3 <path>` for SymPy) and iterate until it exits 0 with all in-file checks passing.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage217_full_interior_five_coordinate_simplex_optimizer_and_finite_candidate_reduction_sympy_audit.txt`

**Issue:** The committed SymPy output `.txt` (mtime 2026-05-11, epoch 1778525386) predates the current SymPy script (mtime 2026-06-03, epoch 1780523951). The stale transcript carries an out-of-date banner — it prints `STAGE 200 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER` (output line 11) and header date `2026-05-11`, while the current `.py` banner (line 35) prints `STAGE 217 — ...`. The numeric content (Bézout 162, projected 750, all-zero residuals, EXIT_CODE 0) is still correct and consistent with the current script, so this is a freshness/label refresh only — no script-logic change.

**Required change:**
No edit to the `.py` source. Re-run the SymPy script and re-capture its stdout to overwrite the committed `.txt` so the banner line reads `STAGE 217 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER` and the header date/mtime are current. (The orchestrator's independent re-run regenerates this transcript; this directive records the cause.)

**Verification command:**
After re-run, the verifier confirms the `.txt` header date is current, line ~11 reads `STAGE 217 — UNIQUE FIVE-COORDINATE SIMPLEX INTERIOR OPTIMIZER`, `Lifted Bezout bound = 162`, `Projected one-chart Bezout bound = 750`, every residual line is `= 0`, and `EXIT_CODE: 0`.
