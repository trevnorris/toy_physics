---
unit_id: 167
batch: V.1
created_at: 2026-05-28T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 167

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — hardcoded_result (mislabeled banner)

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py:29`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl:26`

**Issue:** Both scripts print a banner that says `STAGE 150 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION`, but this is the stage-167 audit. The wrong stage number propagates into both saved transcripts and mislabels the captured evidence. The math is unaffected; this is a pure label correction.

**Required change:**
1. In the SymPy script, line 29, change:
   - before: `banner("STAGE 150 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION")`
   - after:  `banner("STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION")`
2. In the Mathematica script, line 26, change:
   - before: `banner["STAGE 150 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION"];`
   - after:  `banner["STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION"];`

Do not change any other line. Do not touch the inline stage-number comments (they are cosmetic renumbering churn and out of scope).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 167` and `redteam exec-mathematica 167` and confirm: the saved SymPy transcript line 11 and Mathematica transcript line 11 both read `STAGE 167 — BUNDLE TRANSPORT AND TANGENT-COMPENSATION`, both scripts exit 0, and all checks still PASS (6 sympy `expect_zero`, 15 mathematica `expectZero`/PASS lines).

## Applied: F1
files_changed:
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage167_bundle_transport_tangent_compensation_sympy_audit.py`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage167_bundle_transport_tangent_compensation_mathematica_audit.wl`
summary: Corrected the mislabeled banner from STAGE 150 to STAGE 167 in both the SymPy and Mathematica audit scripts.
deviation: none
