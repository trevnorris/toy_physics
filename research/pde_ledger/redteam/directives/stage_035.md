---
unit_id: 035
batch: II.1
created_at: 2026-05-21T21:52:09Z
findings_count: 2
stop_cold: null
applied: true
applied_at: 2026-05-21T23:35:06Z
findings_applied: 2
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 035

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:111`

**Issue:** Line 111 contains `expectZero["support split residual", gBReqSqOverVarpi2 - (alphaReqTarget - alphaMix)];`. The first argument to this `expectZero` is, by the definition on line 109 (`gBReqSqOverVarpi2 = FullSimplify[alphaReqTarget - alphaMix, ...]`), exactly `FullSimplify[E] - E` where `E = alphaReqTarget - alphaMix`. After the internal `FullSimplify` in `expectZero`, this is identically zero by construction. The assertion cannot fail no matter what `alphaReqTarget` or `alphaMix` evaluate to, so it verifies nothing. The SymPy counterpart does not contain this assertion (lines 96-99 of the `.py` only print the expression). Remove the Mathematica-only tautological assertion. Do not introduce a substitute claim, because no `g_B,req^2/varpi^2` closed-form target is defined in either script.

**Required change:**

Delete line 111 of `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl` in its entirety. The surrounding context (the assignment on line 109 and the `Print` on line 110) is retained so the value is still computed and emitted in the output transcript.

Before:
```
alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
gBReqSqOverVarpi2 = FullSimplify[alphaReqTarget - alphaMix, Assumptions -> $Assumptions];
Print["g_B,req^2 / varpi^2 = ", fmt[gBReqSqOverVarpi2]];
expectZero["support split residual", gBReqSqOverVarpi2 - (alphaReqTarget - alphaMix)];
```

After:
```
alphaMix = FullSimplify[Chi^2/(OmegaU^2*Delta0), Assumptions -> $Assumptions];
gBReqSqOverVarpi2 = FullSimplify[alphaReqTarget - alphaMix, Assumptions -> $Assumptions];
Print["g_B,req^2 / varpi^2 = ", fmt[gBReqSqOverVarpi2]];
```

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 035` and confirm:
- The output transcript no longer contains the lines `support split residual = 0` and `PASS: support split residual` (these currently appear at lines 37-38 of the saved output `.txt`).
- The script still exits 0.
- The `Print["g_B,req^2 / varpi^2 = ", ...]` line still emits the value.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl`
- summary: Removed the Mathematica-only tautological support split residual assertion while retaining the computed and printed value.
- deviation: none

## F2 — mathematica_transliteration

**Target:** `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl:50-54,64-69,74-83,92-101,113-119`

**Issue:** The Mathematica script's closed-form targets (`fTarget`, `dFTarget`, `softConstTarget`, `alphaReqTarget`, `alphaCritTarget`, `fSeriesTarget`, `alphaSeriesTarget`) are hard-coded as the same literal polynomials/rational functions used in the SymPy script. Worse, the LHS quantities they are compared against are themselves derived via the same algebraic choreography (same intermediate `x = A*xi`, `deltaKSub = A*delta`, same `alphaX`, same `nX`). As a result, the two engines do not provide independent cross-checks of those literal coefficients: both scripts simply confirm that `D[f, xi]` (or `Limit[...]`, or `Series[...]`) computed by their respective engines equals the SAME literal target. If the literal had a wrong coefficient (e.g., `121*xi^3` should be some other value), both engines would compute their own `D[f, xi]` correctly and both would compare against the same wrong literal — neither would catch it.

The fix is structural: do not change the literal targets themselves (those ARE the claim under test), but ensure the LHS in each `expectZero` is computed by Mathematica from the physical premises (`nX`, `f`, `alphaX`) and NOT by `FullSimplify`-ing the same literal target. This makes the Mathematica script a genuine independent verification of the literal coefficients.

**Required change:**

Edit the Mathematica script as follows. Replace the named ranges only — do not touch other lines.

1. **Lines 64-71** (currently): replace
   ```
   dF = FullSimplify[D[fTarget, xi], Assumptions -> $Assumptions];
   dFTarget = FullSimplify[
     (9*delta + 11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)/
       (81*(1 - xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3),
     Assumptions -> $Assumptions
   ];
   Print["dF/dxi = ", fmt[dF]];
   expectZero["dF/dxi - manifestly positive form", dF - dFTarget];
   ```
   with
   ```
   dF = FullSimplify[D[f, xi], Assumptions -> $Assumptions];
   dFTarget = FullSimplify[
     (9*delta + 11*xi)^3*(81*delta^3 + 72*delta^2 + 189*delta^2*xi + 297*delta*xi^2 + 121*xi^3)/
       (81*(1 - xi)^2*(9*delta^2 + 18*delta*xi + 11*xi^2)^3),
     Assumptions -> $Assumptions
   ];
   Print["dF/dxi = ", fmt[dF]];
   expectZero["dF/dxi - manifestly positive form", dF - dFTarget];
   ```
   That is, change `D[fTarget, xi]` to `D[f, xi]` on line 64. This makes `dF` an independent computation from `f` (which was derived from `nX` on line 50), not a re-simplification of the same literal `fTarget`.

2. **Lines 74-77**: replace
   ```
   softConst = Block[
     {$Assumptions = Element[{delta}, Reals] && delta > 0},
     FullSimplify[Limit[(1 - xi)*fTarget, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
   ];
   ```
   with
   ```
   softConst = Block[
     {$Assumptions = Element[{delta}, Reals] && delta > 0},
     FullSimplify[Limit[(1 - xi)*f, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
   ];
   ```
   Change `fTarget` to `f` inside the `Limit`.

3. **Lines 85-89**: replace
   ```
   epsSoft = Symbol["epsSoft"];
   softScaledSeries = FullSimplify[
     Normal[Series[(epsSoft*fTarget /. xi -> 1 - epsSoft), {epsSoft, 0, 0}]],
     Assumptions -> delta > 0 && epsSoft > 0
   ];
   ```
   with
   ```
   epsSoft = Symbol["epsSoft"];
   softScaledSeries = FullSimplify[
     Normal[Series[(epsSoft*f /. xi -> 1 - epsSoft), {epsSoft, 0, 0}]],
     Assumptions -> delta > 0 && epsSoft > 0
   ];
   ```
   Change `fTarget` to `f` inside the `Series`.

4. **Line 92**: replace
   ```
   alphaReq = FullSimplify[alphaX, Assumptions -> $Assumptions];
   ```
   leave this line unchanged; `alphaReq` is already derived from `alphaX` (the physical definition), not from `alphaReqTarget`. No edit needed here.

5. **Lines 97-100**: replace
   ```
   alphaCrit = Block[
     {$Assumptions = Element[{A, delta}, Reals] && A > 0 && delta > 0},
     FullSimplify[Limit[alphaReqTarget, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
   ];
   ```
   with
   ```
   alphaCrit = Block[
     {$Assumptions = Element[{A, delta}, Reals] && A > 0 && delta > 0},
     FullSimplify[Limit[alphaReq, xi -> 1, Direction -> "FromBelow"], Assumptions -> $Assumptions]
   ];
   ```
   Change `alphaReqTarget` to `alphaReq` (the Mathematica-derived `alphaX`).

6. **Lines 113-114**: replace
   ```
   fSeries = FullSimplify[Normal[Series[fTarget, {xi, 0, 2}]], Assumptions -> delta > 0];
   ```
   with
   ```
   fSeries = FullSimplify[Normal[Series[f, {xi, 0, 2}]], Assumptions -> delta > 0];
   ```
   Change `fTarget` to `f`.

7. **Line 118**: replace
   ```
   alphaSeries = FullSimplify[Normal[Series[alphaReqTarget, {xi, 0, 2}]], Assumptions -> A > 0 && delta > 0];
   ```
   with
   ```
   alphaSeries = FullSimplify[Normal[Series[alphaReq, {xi, 0, 2}]], Assumptions -> A > 0 && delta > 0];
   ```
   Change `alphaReqTarget` to `alphaReq`.

The pattern is uniform: anywhere the LHS of an `expectZero` was computed by acting on a literal target, switch the input to the Mathematica-derived physical quantity (`f` for `fTarget`, `alphaReq` for `alphaReqTarget`). Do NOT change the literal targets themselves (those define the claim being tested). Do NOT change the SymPy script.

**Verification command:**

After Codex applies, the verifier will run `redteam exec-mathematica 035` and confirm:
- The script still exits 0.
- All `expectZero` checks still print residual `= 0` and `PASS`.
- The `Print["dF/dxi = ", fmt[dF]]` output line should be unchanged (mathematically equivalent), since `D[f, xi]` and `D[fTarget, xi]` produce the same result when `f - fTarget` is verified zero earlier.
- The structural property that now matters: if an auditor were to introduce an artificial error into one literal target (e.g., change `121*xi^3` to `122*xi^3` in `dFTarget`), the script would now FAIL with a nonzero residual, because `dF = D[f, xi]` is derived independently. (The verifier does not need to actually run this perturbation — the code structure guarantees it.)

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage035_dimensionless_normalization_locus_mathematica_audit.wl`
- summary: Switched the Mathematica LHS computations for derivatives, limits, and series from literal targets to the derived `f` and `alphaReq` quantities.
- deviation: none
