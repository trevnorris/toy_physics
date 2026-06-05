---
unit_id: 072
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 072

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named. Do NOT touch paper.tex, notes/, or any prose documents.

After editing (here: re-running the scripts), the scripts must exit 0 with all in-file checks passing; the orchestrator independently re-runs afterward.

## F1 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt`

**Issue:**
Both committed transcripts predate the current scripts and were produced by the pre-renumber (055) revision. The current scripts emit `STAGE 072` banners (sympy py:24, mathematica wl:26; sympy ledger header py:110), but the committed sympy transcript reads `STAGE 55` (line 3) and `STAGE 55 THEOREM LEDGER` (line 31), and the committed mathematica transcript reads `STAGE 055` (line 3). mtimes confirm staleness: both `.txt` (`2026-05-22`) are ~4 days older than their scripts (`2026-05-26`). The symbolic/numeric content of the transcripts otherwise matches what the current scripts produce; this is a banner/freshness drift only.

**Required change:**
No source edit. Re-run both scripts and overwrite the committed transcripts:
- `python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.py` → save stdout to `scripts/output/moving_throat_pde_stage072_explicit_branch_thresholds_sympy_audit.txt`
- `math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.wl` → save stdout to `mathematica/output/moving_throat_pde_stage072_explicit_branch_thresholds_mathematica_audit.txt`

After refresh, the sympy transcript line 3 must read `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES` and the ledger header `STAGE 072 THEOREM LEDGER`; the mathematica transcript line 3 must read `STAGE 072 — EXPLICIT BRANCH THRESHOLD SURFACES`. All residual lines remain `0`/PASS and both scripts exit 0.

**Verification command:**
The verifier re-runs `redteam exec-sympy 072` and `redteam exec-mathematica 072` and confirms the refreshed transcripts carry the `STAGE 072` banner and all checks pass.
