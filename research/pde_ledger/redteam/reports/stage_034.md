---
unit_id: 034
batch: II.1
auditor_model: claude-opus-4-7[1m]
audit_date: 2026-05-21T00:00:00Z
verdict: findings
stop_cold: null
findings_count: 1
scripts_checked:
  sympy: present
  mathematica: present
  engines_agree: true
  outputs_fresh: true
---

# Audit unit 034 red-team report

## Files reviewed

- sympy: `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py`
- mathematica: `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl`
- sympy output: `/var/projects/toy_physics/research/pde_ledger/scripts/output/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.txt`
- mathematica output: `/var/projects/toy_physics/research/pde_ledger/mathematica/output/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.txt`

## What the script claims to verify

The unit checks the closed-form "softening-depth" normal form for the selected eigenvalue branch of a rank-1 two-mode secular problem. Specifically: (1) the closed-form expressions `alpha(x) = x(x+DeltaK)/(kappa0_sq*(x+DeltaK)+kappa1_sq*x)`, `s_-(x)` and `N_-(x)` (all functions of the softening depth `x`) are algebraically equal to the original-variable forms `alpha(lambda)`, `s_-(lambda)`, `N_-(lambda)` under the substitution `lambda -> A - x`; (2) `d alpha/dx` matches a hand-written manifestly-positive rational form, and the reciprocity `s_-(x) * d alpha/dx == 1` holds exactly; (3) solving the support-loading partition `gB^2/varpi^2 + alpha_mix == alpha(x)` for the loading reproduces `alpha(x) - alpha_mix`. Both engines exercise identical claims and report PASS for every assertion.

## Assertion inventory

| # | Script | Line | Form | Anchored to claim? |
|---|---|---|---|---|
| A1 | sympy | 61 | `simplify(alpha_lam.subs(lam, A-x) - alpha_x) == 0` | yes |
| A2 | sympy | 62 | `simplify(s_lam.subs(lam, A-x) - s_x) == 0` | yes |
| A3 | sympy | 63 | `simplify(N_lam.subs(lam, A-x) - N_x) == 0` | yes |
| A4 | sympy | 70 | `simplify(dalpha_dx - dalpha_target) == 0` | yes |
| A5 | sympy | 71 | `simplify(s_x * dalpha_dx - 1) == 0` | yes |
| A6 | sympy | 84-87 | `gBreq_sq_over_varpi2 - (alpha_x - alpha_mix) == 0` | no (tautological — by construction) |
| A7 | sympy | 88-91 | `alpha_mix + gBreq_sq_over_varpi2 - alpha_x == 0` | no (tautological — by construction) |
| A8 | mathematica | 65 | `expectZero[(alphaLambda /. lambda -> A-x) - alphaX]` | yes |
| A9 | mathematica | 66 | `expectZero[(sLambda /. lambda -> A-x) - sX]` | yes |
| A10 | mathematica | 67 | `expectZero[(nLambda /. lambda -> A-x) - nX]` | yes |
| A11 | mathematica | 75 | `expectZero[dAlphaDx - dAlphaTarget]` | yes |
| A12 | mathematica | 76 | `expectZero[sX*dAlphaDx - 1]` | yes |
| A13 | mathematica | 91-94 | `expectZero[gBReqSqOverVarpi2 - (alphaX - alphaMix)]` | no (tautological — by construction) |
| A14 | mathematica | 95-98 | `expectZero[alphaMix + gBReqSqOverVarpi2 - alphaX]` | no (tautological — by construction) |

## Findings

### F1 — tautological_check

**Severity:** low
**Files:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage034_softening_depth_normal_form_sympy_audit.py:74-91`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage034_softening_depth_normal_form_mathematica_audit.wl:78-98`

**What's wrong:**
In SymPy lines 76-78 the script constructs

```
alpha_total = gB_sq_over_varpi2 + alpha_mix
gBreq_sq_over_varpi2 = sp.solve(sp.Eq(alpha_total, alpha_x), gB_sq_over_varpi2)[0]
```

i.e. it linearly solves `gB + alpha_mix == alpha_x` for the single unknown `gB`. The unique solution returned by SymPy is, by definition, `alpha_x - alpha_mix`. Lines 84-91 then assert

```
gBreq_sq_over_varpi2 - (alpha_x - alpha_mix) == 0
alpha_mix + gBreq_sq_over_varpi2 - alpha_x == 0
```

Both residuals are algebraically guaranteed to be zero because `gBreq_sq_over_varpi2` was just constructed as the solution to that exact linear equation. Neither assertion can fail unless SymPy's solver is broken; they are not exercising any physical claim about the support loading. The Mathematica file mirrors the same pattern at lines 85-98 with `Solve[alphaMix + gBsqOverVarpi2 == alphaX, ...]` followed by the identical pair of residual checks.

**Why this matters:**
The docstring item 4 advertises the script as verifying "the exact required support-loading formula `g_B^2/varpi^2` in the softening-depth variable", but the two assertions claimed to do this only verify that the CAS's linear solver returned a correct rearrangement. They contribute zero adversarial signal: even if `alpha_x` or `alpha_mix` were defined incorrectly, the checks would still pass. A non-trivial check would, for example, evaluate the support-loading formula at a concrete numeric point or compare the lambda-form and x-form of the loading equation independently.

**Required change:**

Two-part edit, applied symmetrically to both engines.

1. **Remove the tautological residual checks.** In the SymPy file delete lines 84-91 (the two `expect_zero` calls). In the Mathematica file delete lines 91-98 (the two `expectZero` calls). Keep the `print`/`Print` of `gBreq_sq_over_varpi2` / `gBReqSqOverVarpi2` as informational output — these are not assertions.

2. **Replace with a substantive cross-form verification of the loading equation.** Add, immediately after the deleted block, a check that solving the loading equation in the original `lambda` variable yields the same result (modulo `lambda = A - x`). Concretely, in SymPy:

   ```python
   gBreq_lambda = sp.simplify(sp.solve(sp.Eq(gB_sq_over_varpi2 + alpha_mix, alpha_lam), gB_sq_over_varpi2)[0])
   expect_zero(
       "lambda-form vs x-form support loading",
       gBreq_lambda.subs(lam, A - x) - gBreq_sq_over_varpi2,
   )
   ```

   In Mathematica:

   ```mathematica
   gBSolutionLambda = Solve[alphaMix + gBsqOverVarpi2 == alphaLambda, gBsqOverVarpi2, Reals];
   If[Length[gBSolutionLambda] != 1, fail["support-loading lambda-form solve count", Length[gBSolutionLambda]]];
   gBReqLambda = FullSimplify[gBsqOverVarpi2 /. First[gBSolutionLambda], Assumptions -> $Assumptions];
   expectZero[
       "lambda-form vs x-form support loading",
       (gBReqLambda /. lambda -> A - x) - gBReqSqOverVarpi2
   ];
   ```

   This check is non-tautological because `gBreq_lambda` is constructed by independently solving a different equation (in `lambda`, using `alpha_lam`'s closed form), and the residual genuinely depends on both the substitution `lambda -> A - x` and the algebraic equivalence `alpha_lam.subs(lam, A-x) == alpha_x` being correctly carried through the solve.

**Verification:**
After Codex applies the edits, the SymPy output should show a new line
```
lambda-form vs x-form support loading = 0
```
and the lines `solved g_B,req^2/varpi^2 - (alpha(x) - alpha_mix) = 0` and `alpha_mix + g_B,req^2/varpi^2 - alpha(x) = 0` must be gone. The Mathematica output mirrors with a new `PASS: lambda-form vs x-form support loading` line and the prior two PASS lines for the tautological residuals must be gone. Both scripts must still exit 0.

## Independent-derivation check (Mathematica)

The Mathematica file uses parallel variable names (`alphaLambda`/`alpha_lam`, `s1`/`S1`, `sLambda`/`s_lam`, `alphaX`/`alpha_x`, etc.) and applies the same algebraic steps in the same order. However, what the unit actually verifies is the algebraic identity between the closed-form expression in `x` (the "claim") and the substituted lambda-form. Both engines independently invoke their own CAS simplifiers (`sp.simplify(sp.expand(...))` vs `FullSimplify[Together[Expand[...]]]`) to reduce the residuals to zero. Since the closed-form expressions ARE the claim being verified, both scripts hardcoding them is appropriate, and the engine-independent simplification of the residual constitutes a real cross-check. I do not flag `mathematica_transliteration`. (Borderline call; the surrounding choreography is parallel, but the load-bearing simplifier step is engine-native in each.)

## Engine cross-check

Both engines produce identical final closed forms (e.g. `alpha(x) = x*(DeltaK+x)/(kappa0_sq*(DeltaK+x)+kappa1_sq*x)`, `s_-(x) = (kappa0_sq*(DeltaK+x)+kappa1_sq*x)^2 / (kappa0_sq*(DeltaK+x)^2 + kappa1_sq*x^2)`, `N_-(x) = beta0*(...)^4 / (kappa0_sq*(A-x)*(...)^2)`). Every residual asserted in both files reduces to exactly 0 in both engines. No engine disagreement. Outputs (May 11) are newer than scripts (Apr 3 / May 11) — fresh.

## Verdict justification

Five of the seven assertions (A1-A5 in sympy, A8-A12 in mathematica) substantively verify the unit's algebraic claims about the softening-depth normal form: the lambda-to-x equivalence of `alpha`, `s_-`, `N_-`, the derivative identity, and the reciprocity `s_x * dalpha/dx == 1`. These survived attack — I tried to find domain assumption mismatches (sympy declares `x` nonnegative and `lam` positive, but never substitutes at `x=A` or `lam=A` where divergences would occur; Mathematica is stricter with `0 <= x < A` and `0 < lambda < A`) and an algebraic identity flaw (substituted `lambda = A - x` by hand in `alpha_lam = 1/(ks0/(A-lam) + ks1/(A+DeltaK-lam))` and got `1/(ks0/x + ks1/(DeltaK+x))` which simplifies to `x*(DeltaK+x)/(ks0*(DeltaK+x)+ks1*x)` — matches `alpha_x` exactly).

The support-loading "verification" (A6-A7 in sympy, A13-A14 in mathematica) is tautological: it asks the CAS to solve a one-unknown linear equation and then asserts the rearrangement of that same equation is zero. This cannot fail. Verdict is `findings` with one low-severity `tautological_check`. Not stop-cold; the fix is local and does not affect downstream units (the closed forms `alpha_x`, `s_x`, `N_x` carried forward are unchanged).

## Self-test notes

- **Variable independence**: Proposed numeric residual check substitutes all symbols simultaneously; no `sp.diff` against a non-occurring variable.
- **Symmetry/parity**: N/A — no unbounded-domain integrals here.
- **Trivial-case pre-check**: At `ks0=2, ks1=3, A=5, DeltaK=1, x=1, Chi=1, OmegaU=1, Delta0=2`: `alpha_x = 1*2/(2*2+3*1) = 2/7`, `alpha_mix = 1/(1*2) = 1/2`, so `gBreq = 2/7 - 1/2 = 4/14 - 7/14 = -3/14`. The residual `gBreq - (alpha_x - alpha_mix)` is then `-3/14 - (-3/14) = 0` — passes trivially, confirming the construction. The lambda-form cross-check: `alpha_lam` at `lambda = A - x = 4` is `(5-4)*(5+1-4)/(2*(5+1-4)+3*(5-4)) = 1*2/(2*2+3) = 2/7` — matches `alpha_x`, so `gB_lambda(lam=4) = gBreq(x=1) = -3/14`. Non-trivial because the lambda-form is constructed in independent symbolic variables and the residual genuinely depends on the substitution being correct.
- **Path specifications**: N/A — no missing-script findings.
