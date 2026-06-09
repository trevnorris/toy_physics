---
unit_id: 207
batch: VI.1
created_at: 2026-06-09T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 207

Apply the finding below. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead.

Do NOT introduce new features, refactors, or stylistic changes. Do NOT touch paper.tex, notes/, or any prose document. The red-team only modifies scripts and their captured outputs.

After editing, RUN the affected script (`python3 <path>`) and confirm it exits 0 with all in-file checks passing. The orchestrator independently re-runs afterward.

## F1 — stale_output

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage207_primitive_ray_hessian_envelopes_and_certified_ray_table_sympy_audit.txt`

**Issue:** The committed SymPy output transcript (mtime 2026-05-11 12:49:16) is older than the SymPy script (mtime 2026-06-03 15:59:11) and still prints the stale numbering banner `STAGE 190 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE` (output line 11) and `STAGE 190 SYMPY AUDIT PASSED` (line 242). The current script (`scripts/moving_throat_pde_stage207_..._sympy_audit.py` line 35 / line 244) prints `STAGE 207 …`. The rest of the transcript content (sections I–VII, all residuals = 0, `EXIT_CODE: 0`) already matches the current script; only the stage label is stale. No script-math change is required — this is a refresh of the captured output.

**Required change:**
Re-run the SymPy script and re-capture its stdout to the output `.txt` so the banner reads `STAGE 207 …` instead of `STAGE 190 …`. Do not edit the `.py` math. Do not hand-edit the `.txt` banner in isolation — regenerate the whole transcript from the current script so it is a faithful capture.

**Verification command:**
After Codex re-captures, the verifier runs `redteam exec-sympy 207` and confirms the fresh `scripts/output/..._sympy_audit.txt` shows `STAGE 207 — PRIMITIVE-RAY HESSIAN ENVELOPES AND CERTIFIED RAY TABLE` and `STAGE 207 SYMPY AUDIT PASSED`, with mtime newer than the `.py` and `EXIT_CODE: 0`.
