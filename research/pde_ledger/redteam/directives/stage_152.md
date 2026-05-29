---
unit_id: 152
batch: IV.6
created_at: 2026-05-27T00:00:00-06:00
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 152

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts (except when a follow-up directive explicitly authorizes a paper-side edit after user resolution).

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py:115-120`

**Issue:**
The SymPy script computes the four load-bearing scale deliverables from the notes — `delta Pi_act`, `delta T_m,act`, `lambda_eff^(Pi)`, `lambda_eff^(T)` — but only `print`s them. Its sole `expect_close` is a self-consistent covariance moment identity (already at lines 118-120). The Mathematica counterpart (`mathematica/moving_throat_pde_stage152_family1_actual_correction_mathematica_audit.wl:126-129`) asserts all four with anchored `expectApprox` calls against the exact notes-quoted constants. The SymPy script should mirror those anchored assertions so a future regression in any of the canonical inputs (`Pi_star`, `Sigma_m_star`, `gprime_star`, `AT`, `BT`, the lambda baselines) flips the SymPy exit code, not just the Mathematica exit code.

**Required change:**

Locate the line in `scripts/moving_throat_pde_stage152_family1_actual_correction_sympy_audit.py` that currently reads:

```python
print("lambda_eff^(Pi)   =", lam_Pi)
print("lambda_eff^(T)    =", lam_T)
```

(currently lines 114-115). Immediately after these two `print` lines, and before the existing `# Internal consistency checks` comment block at line 117, insert exactly the following four lines:

```python
expect_close("delta Pi_act scale", float(deltaPi), 0.907084414842908, tol=1e-6)
expect_close("delta Tm_act scale", float(deltaT), 0.271653979462338, tol=1e-6)
expect_close("lambda_eff^(Pi) scale", float(lam_Pi), 0.380487632771110, tol=1e-6)
expect_close("lambda_eff^(T) scale", float(lam_T), 0.378939241176339, tol=1e-6)
```

Notes:
- Use `float(...)` because `deltaPi`, `deltaT`, `lam_Pi`, `lam_T` are `mpmath.mpf` values; `expect_close` is defined with `abs(val - target)` and a float target, so cast before the subtraction to avoid mpmath-vs-float comparison issues.
- Do NOT modify the existing `expect_close("delta g from sigma1-linearized covariance check", ...)` at lines 118-120; keep it intact.
- Do NOT change tolerances, do NOT renumber comments, do NOT touch any other line.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 152` and confirm that:
1. The script still exits 0.
2. The transcript now contains four new lines of the form `delta Pi_act scale = 0.907... (target 0.907084414842908, diff <1e-9)` and similarly for `delta Tm_act scale`, `lambda_eff^(Pi) scale`, `lambda_eff^(T) scale`.
3. The existing `delta g from sigma1-linearized covariance check` line is still present.
