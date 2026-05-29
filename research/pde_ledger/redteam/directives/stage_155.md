---
unit_id: 155
batch: IV.6
created_at: 2026-05-27T00:00:00Z
findings_count: 2
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 155

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:92`

**Issue:** The SymPy script computes and prints g_fp, S_fp, R_fp, Pi_fp (lines 89–92) but the only assertion in the script (lines 100–101) tests the transport-law identity, which is internally self-consistent — it would pass even if the converged moments were silently off. The paper notes (§1, §2) elevate g_fp ≈ 0.693352419668063, S_fp ≈ 0.6216013167514007, R_fp ≈ 0.2827139049082381, and Pi_fp ≈ 1.4885734438300713 to load-bearing deliverables, and the Mathematica twin (lines 101–104) already anchors them with `expectApprox`. The SymPy side is asymmetric.

**Required change:**
In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py`, insert exactly the following four lines between the existing line 92 (`print("Pi_fp =", Pi_fp)`) and the existing blank line preceding line 94 (`dg = sp.Float(g_fp - g_star, 30)`). Keep the indentation at module level (no leading spaces):

```python
assert abs(g_fp - 0.693352419668063) < 1e-12, ("g_fp mismatch", g_fp)
assert abs(S_fp - 0.6216013167514007) < 1e-12, ("S_fp mismatch", S_fp)
assert abs(R_fp - 0.2827139049082381) < 1e-12, ("R_fp mismatch", R_fp)
assert abs(Pi_fp - 1.4885734438300713) < 1e-12, ("Pi_fp mismatch", Pi_fp)
```

The four numeric constants are taken verbatim from the existing Mathematica `expectApprox` calls (lines 101–104 of the `.wl` file) and from notes §1/§2; do not change them.

**Claim manifest:** n/a (additive assertions on an existing script).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 155` and confirm the four new `assert` lines appear in the script AND the script exits 0.

## F2 — paper_misalignment

**Subtype:** target_mismatch (banner label)

**Paper side:**
- `/var/projects/toy_physics/research/pde_ledger/paper/stages/stage_155.tex:1` quote: "`\section[Stage 155]{Stage 155: Family-1 Co-Evolving Fixed Point at Frozen Canonical Traction}`"

**Script side:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py:80` quote: `banner("STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT")`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl:26` quote: `banner["STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"];`

**Issue:**
Both audit scripts print a banner that says "STAGE 138" instead of "STAGE 155". This appears to be a leftover from copy-paste from an earlier stage's audit script. The paper card declares this unit to be stage 155. The math is unaffected, but the transcript labels are wrong, complicating downstream tracker correlation.

This sub-finding is purely script-side (the paper card is correct; the scripts are wrong about their own identity). It is therefore safe for Codex to apply mechanically without user resolution.

**Required change:**
- In `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage155_frozen_traction_fixedpoint_sympy_audit.py` line 80, change the string literal `"STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"` to `"STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"`.
- In `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage155_frozen_traction_fixedpoint_mathematica_audit.wl` line 26, change the string literal `"STAGE 138 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"` to `"STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT"`.

Do NOT touch the SymPy docstring at line 11 ("transport law from Stage 154"). The Stage 154 vs. notes §3 "Stage 239" cross-reference question is a paper/notes-side textual matter that the user will resolve separately; the script-side docstring is acceptable as-is for this directive.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 155` and `redteam exec-mathematica 155`. Fresh transcripts should begin with `STAGE 155 — FROZEN-TRACTION CO-EVOLVING FIXED POINT` instead of `STAGE 138 — …`.
