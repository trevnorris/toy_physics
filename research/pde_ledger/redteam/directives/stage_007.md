---
unit_id: 007
batch: I.1
created_at: 2026-05-25T00:00:00Z
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 007

The single finding is `stale_output`. **No script edits are required.** The verifier's normal `redteam exec-sympy 007` and `redteam exec-mathematica 007` refresh will land fresh captures and resolve the finding automatically. Codex should take no action on this directive.

## F1 — stale_output

**Target:** (no script edit)
- `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt`

**Issue:**
Both saved `.txt` captures predate the current scripts (capture mtimes 2026-05-21; script mtimes 2026-05-25), and their content reflects the pre-H(w) script — the sympy capture has no `xi_eff` print lines, and the mathematica capture has no M2b/M4b/M4c/M7b/M8b/M8c/M10b/M10c/M11b/M11c entries.

**Required change:**
None for Codex. The verifier executes the scripts as part of its normal sweep, which overwrites the `.txt` files. The math in both scripts already passes (substantively verified by the auditor against analytic Gaussian-moment closed forms). Do not edit either script.

**Verification command:**
After the verifier runs `redteam exec-sympy 007` and `redteam exec-mathematica 007`, confirm:
- `scripts/output/moving_throat_pde_stage007_projection_reduction_comparison_sympy_audit.txt` contains a `xi_eff^(proj)` print line (sections 1, 3, 4 all should mention it), and ends in `STATUS: PASS`.
- `mathematica/output/moving_throat_pde_stage007_projection_reduction_comparison_mathematica_audit.txt` contains `PASS: M2b`, `PASS: M4b`, `PASS: M4c`, `PASS: M7b`, `PASS: M8b`, `PASS: M8c`, `PASS: M10b`, `PASS: M10c`, `PASS: M11b`, `PASS: M11c` lines, and ends in `STATUS: PASS`.

If both refreshed captures show PASS for every assertion (including the new H/ξ ones), F1 is resolved and the stage is clean.
