---
unit_id: 084
batch: III.4
created_at: 2026-05-22T00:00:00-06:00
findings_count: 3
stop_cold: null
applied: true
applied_at: 2026-05-25T00:31:44-06:00
findings_applied: 3
findings_blocked: 0
verification_status: pending
---

# Codex directive — unit 084

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:** `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:65-66`

**Issue:** Lines 65-66 assert `xiF1FromUpsilon - 1369*upsilonW == 0` and `xiF1FromTheta - 136900*thetaW == 0`. Because lines 58-59 declare `xiF1FromUpsilon = lambdaEll^2*upsilonW` and `xiF1FromTheta = 100*lambdaEll^2*thetaW` with `lambdaEll = 37` on line 55, these residuals are `(37^2 - 1369)*upsilonW = 0` and `(100*37^2 - 136900)*thetaW = 0` by integer arithmetic alone. The assertion cannot fail short of changing the literal `37` on line 55. Replace with a substantive cross-route consistency check.

**Required change:**
Replace exactly lines 65 and 66 (the two `expectZero[...]` calls for the `Xi_F1` identities) with a single cross-route consistency check:

Before (current lines 65-66):
```
expectZero["Xi_F1(Upsilon) - 1369 Upsilon_w", xiF1FromUpsilon - 1369*upsilonW];
expectZero["Xi_F1(Theta) - 136900 Theta_w", xiF1FromTheta - 136900*thetaW];
```

After:
```
expectZero["Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)", (xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta];
```

Do not change any other lines. Do not delete the print lines 63-64 above.

**Verification command:** After Codex applies, the verifier will re-run `redteam exec-mathematica 084` and confirm: (a) the new `expectZero` line appears at what is now line 65 of the script, (b) the captured output contains `Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta) = 0` followed by `PASS:`, (c) the script exits 0, (d) the number of `PASS:` lines decreases by 1 (from 5 to 4) — except that F2 and F3 below add more PASS lines, so the final count depends on which findings are applied.

## Applied: F1

- files_changed:
  - `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`
- summary: Replaced the two tautological Xi_F1 literal checks with the requested cross-route consistency check.
- deviation: none

## F2 — hardcoded_result

**Target:** `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:80-81`

**Issue:** The two `expectApprox` calls on lines 80-81 subtract two hardcoded numeric `zeta` literals (declared on lines 68-72) and compare against a third hardcoded numeric target. Both sides are pure literals; the check verifies IEEE float subtraction, not physics. The two unused `zeta` values `zetaMinusJ` and `zetaPlusJ` (lines 70-71) are only printed (lines 76-77) and never asserted on. Replace with inequality checks that exercise the ordering structure.

**Required change:**
Replace exactly lines 80 and 81 with four `expectZero` inequality checks:

Before (current lines 80-81):
```
expectApprox["natural-window ordering gap", zetaPlusChi - zetaMinusChi, ToExpression["0.00130621926024`20"], 10^-12];
expectApprox["hard-ceiling gap", zetaMaxF1 - zetaPlusChi, ToExpression["9.6717311`10*^-8"], 10^-12];
```

After:
```
expectZero["chi-window ordering positive (zeta_+^chi > zeta_-^chi)", If[TrueQ[zetaPlusChi > zetaMinusChi], 0, 1]];
expectZero["hard-ceiling gap positive (zeta_max^F1 > zeta_+^chi)", If[TrueQ[zetaMaxF1 > zetaPlusChi], 0, 1]];
expectZero["J-window ordering positive (zeta_+^J > zeta_-^J)", If[TrueQ[zetaPlusJ > zetaMinusJ], 0, 1]];
expectZero["fail-side J below chi (zeta_+^J <= zeta_+^chi)", If[TrueQ[zetaPlusJ <= zetaPlusChi], 0, 1]];
```

Do not change lines 68-78 (the value declarations and prints). The `TrueQ` wrapper ensures undecidable comparisons return `False`, so `If[..., 0, 1]` returns `1` and `expectZero` fails loudly with residual `1`.

**Verification command:** After Codex applies, the verifier will re-run `redteam exec-mathematica 084` and confirm: (a) the four new `expectZero` lines appear in place of the old `expectApprox` pair, (b) the captured output contains four new `... = 0` / `PASS:` line pairs for the inequality checks and no longer contains `natural-window ordering gap diff = ...` / `hard-ceiling gap diff = ...`, (c) the script exits 0.

## Applied: F2

- files_changed:
  - `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`
- summary: Replaced the two hardcoded numeric gap checks with the four requested ordering inequality checks.
- deviation: none

## F3 — insufficient_verification

**Target:** `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:72-73`

**Issue:** Lines 39-46 define `omegaPe`, `zetaPhys`, `zetaReq`, `qMap`, `rQuad` as symbolic objects, but only `zetaReq` and `qMap` are exercised (line 53). The physical Pe→zeta map `zetaPhys` is computed but never pinned to any numeric value, and `kappaF1 = 12321/5`, `etaF1 = 37` (lines 56-57) are declared but never substituted anywhere — they are dead code. Per the upstream derivation in `scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:32-33`, the carried-forward `zeta_max^(F1) = 2.46752922945601` is defined as `limit(A_F1 * Omega^2, Pe, oo)` with `A_F1 = (kappa_F1 + pi^2/4) / (kappa_F1 + y_F1^2)`, `y_F1 = root of y*tan(y) = 37 near 1.53`, and `kappa_F1 = 12321/5`. The status-only consolidation should pin `zetaMaxF1` to this `Pe -> oo` limit of `zetaPhys` to demonstrate that the numeric value is consistent with the physical formula already defined in this script.

**Required change:**
Insert five new lines as a single block immediately AFTER line 72 (the existing `zetaMaxF1 = ToExpression["2.46752922945601`20"];`) and BEFORE the blank line 73. The bindings require `kappaF1`, `etaF1`, and `zetaMaxF1` to already exist, which is why the insertion must follow line 72.

Insert verbatim (after line 72, before line 73):
```
yF1 = y /. FindRoot[y*Tan[y] - 37, {y, 153/100}, WorkingPrecision -> 40];
zetaPhysF1Limit = Limit[zetaPhys /. {kappa -> kappaF1, eta -> etaF1, y -> yF1}, pe -> Infinity];
zetaPhysF1Numeric = N[zetaPhysF1Limit, 20];
Print["zeta_phys(Pe->oo, kappa_F1, eta_F1, y_F1) = ", fmt[zetaPhysF1Numeric]];
expectApprox["zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)", zetaPhysF1Numeric, zetaMaxF1, 10^-10];
```

Do not modify any other lines. After insertion, the blank line that was line 73 becomes line 78, and the subsequent `Print` block (former lines 74-78) shifts down accordingly.

BLOCKED IF: the binding `y_F1 = root of y*tan(y) = 37 near y = 1.53`, `kappa_F1 = 12321/5`, `Pe -> oo` does not produce a numeric value matching `zetaMaxF1 = 2.46752922945601` to 10 digits when computed by Mathematica. If `Limit[..., pe -> Infinity]` returns an unevaluated form, try `Series[... , {pe, Infinity, 0}] // Normal` instead, or wrap the substitution in `FullSimplify` before taking the limit. The upstream SymPy derivation is in `scripts/moving_throat_pde_stage080_family1_zeta_thresholds_sympy_audit.py:28-33`; the math is sound. If the Mathematica numeric does not agree, append a `## Blocked: F3` block describing the disagreement (with both numeric values printed) rather than committing a wrong check.

**Verification command:** After Codex applies, the verifier will re-run `redteam exec-mathematica 084` and confirm: (a) the new `expectApprox` line for `zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)` appears, (b) the captured output contains the corresponding `... diff = <small number>` and `PASS:` lines, (c) the script exits 0. If `Limit[...]` is slow, the elapsed-seconds annotation in the output header will reflect that; this is acceptable as long as it terminates.

## Applied: F3

- files_changed:
  - `mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl`
- summary: Inserted the requested Family-1 Pe-to-infinity zeta_phys limit calculation and approximate match check against zeta_max^(F1).
- deviation: none
