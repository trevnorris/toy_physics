---
unit_id: 204
batch: VI.1
created_at: 2026-06-09T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 204

Apply the finding below. After applying, append an `## Applied: F1` block under it with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts and their saved outputs.

After editing, RUN the affected script (`python3 <path>`) and confirm it exits 0 with all in-file checks passing.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.txt`

**Issue:** The committed SymPy saved output predates the current `.py` (output Date 2026-05-11 vs `.py` mtime 2026-06-03) and was captured from a pre-renumber revision: its banner reads `STAGE 187 — …` (lines 11, 229) whereas the current script prints `STAGE 204` (`.py:35,225`). The captured math still agrees with the current script (no numeric disagreement), so this is informational — the only defect is staleness/wrong stage label in the transcript. The Mathematica output is already fresh and correctly labelled STAGE 204; no `.wl` action needed.

**Required change:**
No source edit. Re-run the SymPy script and overwrite the saved output:
1. `python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.py`
2. Capture stdout to `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage204_explicit_log_ray_sweep_and_scalar_root_predictor_sympy_audit.txt` using the project's standard output-capture path (do not hand-edit the banner; let the current script regenerate it).

**Verification command:**
The verifier will run `redteam exec-sympy 204` and confirm: exit code 0; the regenerated `.txt` banner reads `STAGE 204 — EXPLICIT LOG-RAY SWEEP AND SCALAR ROOT PREDICTOR` (not `STAGE 187`); Date header ≥ `.py` mtime; every check still prints `= 0`.
