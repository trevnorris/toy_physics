---
unit_id: 069
batch: III.3
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 069

Apply the finding below. After applying, append an `## Applied: F1` block with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose document. Only the committed `.txt` transcripts are refreshed here (via re-running the existing scripts).

## F1 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt`

**Issue:** Both committed transcripts predate the scripts (output mtime `May 22 20:05`; both scripts mtime `Jun 3 15:59`) and carry a stale banner self-label `STAGE 052 — FINAL REDUCED SUPPORT/SOURCE VERDICT` (line 3 of each). The current scripts already emit the correct banner `STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT` (`scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py:65`, `mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl:33`). The transcript body otherwise matches the current scripts; only the stage-number self-label is stale.

**Required change:**
No source edit. Re-run both scripts and overwrite the committed transcripts so the banner reads `STAGE 069`:
- `python3 scripts/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.py` → capture stdout to `scripts/output/moving_throat_pde_stage069_final_reduced_verdict_sympy_audit.txt`
- `math -script mathematica/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.wl` → capture stdout to `mathematica/output/moving_throat_pde_stage069_final_reduced_verdict_mathematica_audit.txt`

Both must exit 0 with all checks passing (`= 0` / `PASS:` lines unchanged in substance).

**Verification command:**
The verifier runs `redteam exec-sympy 069` and `redteam exec-mathematica 069` and confirms line 3 of each refreshed `.txt` reads `STAGE 069 — FINAL REDUCED SUPPORT/SOURCE VERDICT`, all assertion lines still pass, and both scripts exit 0.
