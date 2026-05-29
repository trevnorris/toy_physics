---
unit_id: 034
batch: II.1
created_at: 2026-05-21T00:00:00Z
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-21T23:32:50Z
findings_applied: 1
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 034

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py:84-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl:91-98`

**Issue:**
The script defines `gBreq_sq_over_varpi2` (resp. `gBReqSqOverVarpi2`) as the output of `sp.solve(...)` / `Solve[...]` applied to the linear equation `gB + alpha_mix == alpha(x)`. The unique solution returned by the CAS is — by construction — `alpha_x - alpha_mix`. The two subsequent `expect_zero` / `expectZero` calls (sympy lines 84-91, mathematica lines 91-98) then check residuals of the form `gBreq - (alpha_x - alpha_mix)` and `alpha_mix + gBreq - alpha_x`, both of which are algebraically guaranteed to be zero by the construction of `gBreq`. Neither assertion can fail unless the CAS's linear solver is broken; they do not exercise any physical claim about the support-loading formula. Replace them with a non-tautological cross-form check (solve the loading equation in the original `lambda` variable independently, then verify the substitution `lambda -> A - x` recovers the x-form solution).

**Required change:**

**Step 1 — SymPy file.** Replace lines 84-91 of `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`.

Before (lines 84-91):
```python
# Check that solving the microscopic loading equation really reproduces alpha(x).
expect_zero(
    "solved g_B,req^2/varpi^2 - (alpha(x) - alpha_mix)",
    gBreq_sq_over_varpi2 - (alpha_x - alpha_mix),
)
expect_zero(
    "alpha_mix + g_B,req^2/varpi^2 - alpha(x)",
    alpha_mix + gBreq_sq_over_varpi2 - alpha_x,
)
```

After:
```python
# Independent solve of the loading equation in the original lambda variable,
# checked against the x-form solution under lambda = A - x. This is non-trivial:
# it exercises the closed-form alpha_lam together with the substitution.
gBreq_lambda = sp.simplify(
    sp.solve(sp.Eq(gB_sq_over_varpi2 + alpha_mix, alpha_lam), gB_sq_over_varpi2)[0]
)
expect_zero(
    "lambda-form vs x-form support loading",
    gBreq_lambda.subs(lam, A - x) - gBreq_sq_over_varpi2,
)
```

**Step 2 — Mathematica file.** Replace lines 91-98 of `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl`.

Before (lines 91-98):
```mathematica
expectZero[
  "solved g_B,req^2/varpi^2 - (alpha(x) - alpha_mix)",
  gBReqSqOverVarpi2 - (alphaX - alphaMix)
];
expectZero[
  "alpha_mix + g_B,req^2/varpi^2 - alpha(x)",
  alphaMix + gBReqSqOverVarpi2 - alphaX
];
```

After:
```mathematica
gBSolutionLambda = Solve[alphaMix + gBsqOverVarpi2 == alphaLambda, gBsqOverVarpi2, Reals];
If[Length[gBSolutionLambda] != 1, fail["support-loading lambda-form solve count", Length[gBSolutionLambda]]];
gBReqLambda = FullSimplify[gBsqOverVarpi2 /. First[gBSolutionLambda], Assumptions -> $Assumptions];
expectZero[
  "lambda-form vs x-form support loading",
  (gBReqLambda /. lambda -> A - x) - gBReqSqOverVarpi2
];
```

Notes for the Mathematica edit:
- `alphaLambda` was defined earlier (line 34) and is still in scope.
- The current `$Assumptions` (set at lines 79-82) covers `x`, `kappa0Sq`, `kappa1Sq`, etc. but does not include `lambda`. This is OK: the `Solve` treats `lambda` as a free parameter, and after the substitution `lambda -> A - x` no `lambda` symbol remains, so `FullSimplify` under the current assumptions can reduce the residual to zero.
- Do NOT modify `$Assumptions` to re-add `lambda` constraints — the substitution eliminates `lambda` from the residual before simplification.

Do not touch any other line of either file. Do not rename existing variables. Do not change the docstring or the banner strings.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 034` and `redteam exec-mathematica 034` and confirm:
- The new line `lambda-form vs x-form support loading = 0` appears in the SymPy output.
- The new line `PASS: lambda-form vs x-form support loading` appears in the Mathematica output.
- The lines `solved g_B,req^2/varpi^2 - (alpha(x) - alpha_mix) = 0` and `alpha_mix + g_B,req^2/varpi^2 - alpha(x) = 0` (and the corresponding `PASS:` lines in the Mathematica output) are no longer present.
- Both scripts exit 0.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`
  - `mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl`
- summary: Replaced tautological support-loading residual checks with independent lambda-form solves checked against the x-form solution.
- deviation: none
