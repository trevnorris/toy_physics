---
unit_id: 028
batch: II.1
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 028

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage028_loaded_profile_selection_mathematica_audit.wl:47-48` swapped the two `expectZero` assertions for bare `Print` statements that just display `FullSimplify[kappa0^2 - kappa1^2, ...]` and `FullSimplify[2*kappa0*kappa1, ...]`. Confirmed in the diff and in the current file at lines 47-48.

**Assessment:**
The edit matches the directive exactly. The saved Mathematica output now shows `kappa0^2 - kappa1^2 = 56/(9*Pi^2)` and `2 kappa0 kappa1 = (-16*Sqrt[2])/(3*Pi^2)` as printed values with no `PASS:` lines for them, matching the sympy treatment on lines 87-88 of the .py file. Script still exits 0.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
`scripts/moving_throat_pde_stage028_loaded_profile_selection_sympy_audit.py:204-207` now substitutes the hand-written `alpha_crit_expected` (line 197 = `K0*K1/(K1*kappa0**2 + K0*kappa1**2)`) into `det_eff` for both the below-threshold sweep and the at-threshold assertion. The label printed is now `det(alpha_crit_expected) = 0` and the assertion key is `det(alpha_crit_expected)`.

**Assessment:**
The substitution is no longer tautological: `alpha_crit_expected` is a quoted closed form, not a `sp.solve` return value. Because `det_eff = K0 K1 - alpha (K1 kappa0^2 + K0 kappa1^2)` is linear in alpha, substituting the closed form produces `K0 K1 - K0 K1 = 0` only if the closed form is correct; any wrong closed form would leave a non-zero residual. The saved sympy output (lines 91-92) shows `det(alpha_crit_expected) = 0` printed twice (once from the bare `print`, once from `expect_zero`'s own diagnostic line) — correct behavior for this script's helper. The `det_below` factored expression on lines 93-98 of the output is unchanged structurally (`eps*(Keta + 6*TOmega)*(Keta*L^2 + 6*L^2*TOmega + Pi^2*Tw)/L^2`), matching what the directive predicted.

### F3 — tautological_check

**Classification:** resolved

**What changed:**
Both engines replaced the `-rhs - 2 alpha (-kappa0 kappa1)/(...)` sign-rewrite tautology with explicit positivity checks on `kappa0` and `kappa1`:
- `scripts/.../sympy_audit.py:151-153`: `expect_zero("kappa0 sign check (kappa0 > 0)", sp.simplify(kappa0 - sp.Abs(kappa0)))` and `expect_zero("kappa1 sign check (kappa1 < 0)", sp.simplify(kappa1 + sp.Abs(kappa1)))`.
- `mathematica/.../mathematica_audit.wl:95-98`: `Print["-tan(2 theta) = ", ...]` plus `expectZero["kappa0 sign check (kappa0 > 0)", kappa0 - Abs[kappa0]]` and `expectZero["kappa1 sign check (kappa1 < 0)", kappa1 + Abs[kappa1]]`.

**Assessment:**
The new checks are non-tautological. `kappa0 - Abs[kappa0]` collapses to 0 only when `kappa0 >= 0`; with `kappa0 = 2 Sqrt[2]/Pi` (concrete positive numeric literal) the absolute value resolves and the residual is 0. Symmetrically `kappa1 + Abs[kappa1]` is 0 only for `kappa1 <= 0`. If the sign convention on either kappa had been flipped (e.g., a typo making `kappa1 = +4/(3 Pi)`), the corresponding check would NOT collapse to zero. Both engines now show `PASS: kappa0 sign check (kappa0 > 0)` and `PASS: kappa1 sign check (kappa1 < 0)` in their saved outputs (sympy lines 78-79; Mathematica lines 28-31). The `-tan(2 theta)` value is preserved as a `print`/`Print` for inspection.

### F4 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Two genuinely independent Mathematica-native derivations were inserted into `mathematica/.../mathematica_audit.wl`:
- Lines 66-74: `eigvalsDirect = Eigenvalues[kEff]` with `expectZero["Eigenvalues[kEff] sum vs trace", FullSimplify[Total[eigvalsDirect] - (lambdaMinus + lambdaPlus), ...]]` and the corresponding product check.
- Lines 115-119: `alphaCritSolved = alpha /. First[Solve[detEff == 0, alpha]] /. ConditionalExpression[value_, _] :> value` with `expectZero["alphaCrit solved vs ratio closed form", FullSimplify[alphaCritSolved - K0*K1/(K1*kappa0^2 + K0*kappa1^2), ...]]`.

**Assessment:**
Both insertions match the directive. They constitute substantive independent derivations: `Eigenvalues[]` invokes Mathematica's own characteristic-polynomial routine (independent of the hand-written `(trExpected ± Sqrt[disc])/2` ansatz), and `Solve[detEff == 0, alpha]` derives `alphaCrit` by Mathematica's linear-equation solver (independent of the quoted `K0*K1/(K1*kappa0^2 + K0*kappa1^2)`). The saved Mathematica output (lines 16-19 and 41-42) shows all three new PASS lines. The downstream `alphaCrit`, `alphaCritClosed`, and `detEff /. alpha -> alphaCrit` lines are unchanged, as the directive required.

The deviation Codex flagged on F4 ("Stripped Mathematica's `ConditionalExpression` wrapper from the `Solve` result so the inserted `expectZero` check reduces to literal zero") is acceptable: `detEff` is linear in `alpha`, so `Solve` returns a single solution `alpha = K0*K1/(K1*kappa0^2 + K0*kappa1^2)`; any `ConditionalExpression` wrapper would only carry a non-zero-denominator clause (`K1*kappa0^2 + K0*kappa1^2 != 0`), which is implied by the existing real-positive assumptions on the constituent symbols. Stripping the wrapper does not mask a hidden constraint; it just lets `FullSimplify` see a plain expression. The substantive check (does Mathematica's `Solve` produce the same closed form as the hand-written ledger expression?) is preserved.

## Exec log assessment

**SymPy:** exit log not captured (`stage_028_sympy.log` missing — only the diff patch is present in `exec_logs/`). Inferred from the saved output: the script ran to completion. Notable lines from the saved output:
- `trace - expected = 0`, `det - expected = 0`, `characteristic-factorization check = 0`
- `dE/dtheta - stationarity_expected/2 = 0`
- `kappa0 sign check (kappa0 > 0) = 0`, `kappa1 sign check (kappa1 < 0) = 0` (new F3 checks)
- `det(alpha_crit_expected) = 0` (new F2 check)
- `alpha_crit - expected = 0`, `weak-loading coefficient - k0*k1/DeltaK = 0`, `strong-loading limit - tan(2 theta_max) = 0`, `tan(theta_max) + sqrt(2)/3 = 0`

`expect_zero` raises `AssertionError` on a non-zero residual (line 47-48 of the script), so a clean output through the last line implies exit 0.

**Mathematica:** exit log not captured (`stage_028_mathematica.log` missing). Inferred from the saved output, which ends with `Stage 028 Mathematica audit passed.` and was produced by `Exit[0]`. All PASS lines:
- `PASS: trace - expected`, `PASS: det - expected`
- `PASS: Eigenvalues[kEff] sum vs trace`, `PASS: Eigenvalues[kEff] product vs determinant` (new F4 checks)
- `PASS: characteristic factorization`
- `PASS: dE/dtheta - stationarity_expected/2`
- `PASS: kappa0 sign check (kappa0 > 0)`, `PASS: kappa1 sign check (kappa1 < 0)` (new F3 checks)
- `PASS: weak-loading coefficient - kappa0 kappa1/deltaK`
- `PASS: strong-loading limit - tan(2 theta_max)`, `PASS: tan(theta_max) + Sqrt[2]/3`
- `PASS: alphaCrit solved vs ratio closed form` (new F4 check)
- `PASS: alpha_crit - finite-throat closed form`, `PASS: det(alpha_crit)`

`expectZero` calls `fail[...]` (Exit[1]) on a non-zero residual, so the closing `Stage 028 Mathematica audit passed.` implies exit 0.

**Output freshness:** confirmed. Script mtimes are `17:05`, output mtimes are `17:07` (newer). The saved outputs reflect the post-fix scripts.

## Material-change assessment

`material_change`: false.

No derived closed-form quantity changed value. The F2 fix replaced a substitution into `det_eff` with a substitution of the same algebraic value (the closed form Mathematica was already using); F1 demoted assertions to prints; F3 added two new sign assertions that cannot affect any downstream symbolic identity; F4 added two new independent derivations that cross-check the same `lambdaMinus`, `lambdaPlus`, and `alphaCrit` values that were already being used downstream. The `tan(2 theta)`, `weak-loading coefficient`, `strong-loading limit`, `tan(theta_max)`, `alpha_crit`, and `det(alpha_crit*(1-eps))` outputs in both engines are byte-identical to the pre-fix forms recorded in the original auditor's cross-check table.

## Side observations (non-blocking)

- The sympy saved output prints `det(alpha_crit_expected) = 0` on two consecutive lines (91-92). This is because line 206 of the script does `print("det(alpha_crit_expected) =", det_at_expected)` and line 207 calls `expect_zero("det(alpha_crit_expected)", det_at_expected)`, whose helper also prints `<name> = <simplified>`. Identical structure to the pre-fix code (which also printed `det(alpha_crit) = 0` twice); no regression.
- The `exec_logs/` directory only contains `stage_028_diff.patch`. No separate sympy/mathematica log files were written by the orchestrator, but the saved output `.txt` files post-date the script edits, so freshness is satisfied another way.
- Banner in `mathematica_audit.wl:26` still reads `"STAGE 011 — LOADED PROFILE SELECTION"`, which is a legacy label inherited from the renumber. Not part of any finding; the auditor saw this and didn't flag it.

## Verdict justification

All four findings have been addressed in a substantively non-tautological way. F1 demoted tautological numeric identities to bare prints. F2 swapped a `sp.solve`-roundtrip check for a substitution of the independently-quoted closed form. F3 replaced an algebraic identity rewrite with explicit sign checks on `kappa0` and `kappa1` that would catch a flipped sign convention. F4 added Mathematica-native `Eigenvalues[]` and `Solve[]` derivations whose agreement with the hand-written closed forms is now a real cross-engine check rather than a re-simplification of the same expression. Both engines exit cleanly (inferred from the saved outputs ending in the success-banner / completing all assertions without `AssertionError`) and the cross-engine quantities remain identical to the pre-fix forms, so no downstream unit needs re-auditing.
