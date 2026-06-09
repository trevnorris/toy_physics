---
unit_id: 188
batch: V.3
created_at: 2026-06-09T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 188

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose document.

After editing, RUN the affected script (`python3 scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py`) and confirm it exits 0 with all in-file checks passing, regenerating the saved output.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.txt`

**Issue:** The committed SymPy output transcript is stale relative to the current script. The `.py` already prints the correct `STAGE 188` banner (`scripts/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.py:35` and `:198`), but the committed `.txt` still carries the pre-renumber banner `STAGE 171 — BRANCH-OBSERVABLE COMPLETION...` (`scripts/output/...sympy_audit.txt:3`) and `STAGE 171 LEDGER` (`:272`). The `.txt` mtime (Jun 1 11:23) is older than the `.py` mtime (Jun 3 15:59). The body math is otherwise current (every residual already prints `0`); this is a label-only refresh. No source edit is needed — the script source is already correct.

**Required change:**
No edit to the `.py` source. Simply re-run the SymPy script so the saved output is regenerated from the current source, replacing the stale `STAGE 171` banners with `STAGE 188`. If your run harness writes stdout to the committed `.txt` path, ensure the refreshed transcript lands at `scripts/output/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.txt`.

**Verification command:**
After Codex re-runs, the verifier confirms `scripts/output/moving_throat_pde_stage188_branch_observables_completion_sympy_audit.txt` line 3 reads `STAGE 188 — BRANCH-OBSERVABLE COMPLETION AND THE EXACT FIRST-ORDER OBSERVABLE COMPILER`, the ledger banner reads `STAGE 188 LEDGER`, the script exits 0, and the `.txt` mtime post-dates the `.py` mtime.
