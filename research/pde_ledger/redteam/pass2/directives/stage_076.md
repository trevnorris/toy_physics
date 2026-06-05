---
unit_id: 076
batch: III.4
created_at: 2026-06-05T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 076

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — stale_output

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt`

**Issue:** Both committed transcripts predate their scripts (`.txt` mtimes `02:15`/`02:19` vs `.py`/`.wl` mtime `02:38`, all on 2026-05-27) and carry the pre-renumber banner: the sympy transcript line 3 reads `STAGE 59 — EXACT n=5 WALL-DEPTH LOCK` and the mathematica transcript line 3 reads `STAGE 059 — EXACT n=5 WALL-DEPTH LOCK`, whereas the current source banners already read `STAGE 076` (py:30, wl:26). All other transcript content (closed forms, `= 0` residuals, the `n=3` guard residual `3*K*rho**2/4`, the `PASS` lines, and the final `Stage 076 Mathematica audit passed.` line) matches what the current scripts produce, so this is a freshness/label drift only, not a content disagreement.

**Required change:** No source edit. Re-run both scripts and overwrite the committed transcripts so the captured output reflects the current source:
- `python3 /var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.py` → overwrite `scripts/output/moving_throat_pde_stage076_n5_wall_depth_lock_sympy_audit.txt`.
- `math -script /var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.wl` → overwrite `mathematica/output/moving_throat_pde_stage076_n5_wall_depth_lock_mathematica_audit.txt`.

**Verification command:**
After Codex applies, the verifier re-runs `redteam exec-sympy 076` and `redteam exec-mathematica 076` and confirms: sympy transcript line 3 reads `STAGE 076 — EXACT n=5 WALL-DEPTH LOCK`; mathematica transcript line 3 reads `STAGE 076 — EXACT n=5 WALL-DEPTH LOCK`; both scripts exit 0; all residual/PASS/closed-form lines unchanged from the prior transcripts.
