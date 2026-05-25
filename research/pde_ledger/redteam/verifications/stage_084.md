---
unit_id: 084
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T00:40:00-06:00
verdict: verified
sympy_exit: n/a
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 084

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:65` — the two tautological assertions on the old lines 65-66 (`xiF1FromUpsilon - 1369*upsilonW` and `xiF1FromTheta - 136900*thetaW`) are removed and replaced by the single requested cross-route consistency check:
```
expectZero["Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)", (xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta];
```

**Assessment:**
Edit matches the directive verbatim. The new assertion is non-tautological: `xiF1FromUpsilon = lambdaEll^2*upsilonW` and `xiF1FromTheta = 100*lambdaEll^2*thetaW`, so `(xiF1FromUpsilon /. upsilonW -> 100*thetaW) - xiF1FromTheta = 100*lambdaEll^2*thetaW - 100*lambdaEll^2*thetaW = 0` only because the two routes share the same `lambdaEll` and the same `100` factor. Changing either to `100*lambdaEll^2` vs `99*lambdaEll^2`, or using a different exponent on one side, would now be caught — unlike the prior literal-difference form, which only encoded `37^2 = 1369`. Output line 15 confirms `Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta) = 0` followed by `PASS:` on line 16. Print lines 63-64 are preserved per directive. No collateral edits.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:84-87` — the two `expectApprox` calls on the old lines 80-81 (which compared hardcoded float subtractions against hardcoded float targets) are removed; in their place are the four requested `expectZero[..., If[TrueQ[...], 0, 1]]` inequality checks: chi-window ordering, hard-ceiling gap, J-window ordering, and J-vs-chi fail-side consistency.

**Assessment:**
Edit matches the directive verbatim, including the `TrueQ` wrapping. The new assertions exercise the genuine ordering claims `zeta_+^chi > zeta_-^chi`, `zeta_max^F1 > zeta_+^chi`, `zeta_+^J > zeta_-^J`, and `zeta_+^J <= zeta_+^chi`. These are non-tautological: if any of the five imported numeric `zeta` literals were perturbed in a way that flipped an ordering relation, the corresponding `If` branch would return `1` and `expectZero` would `Exit[1]`. The previously unused `zetaMinusJ`/`zetaPlusJ` (printed but never asserted) are now exercised by checks 3 and 4. Output lines 27-34 confirm all four checks pass with residual `0`. Value declarations and print lines 68-78 are preserved per directive.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
`mathematica/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.wl:72-76` — five new lines inserted directly after the old line 72 (`zetaMaxF1 = ...`):
```
yF1 = y /. FindRoot[y*Tan[y] - 37, {y, 153/100}, WorkingPrecision -> 40];
zetaPhysF1Limit = Limit[zetaPhys /. {kappa -> kappaF1, eta -> etaF1, y -> yF1}, pe -> Infinity];
zetaPhysF1Numeric = N[zetaPhysF1Limit, 20];
Print["zeta_phys(Pe->oo, kappa_F1, eta_F1, y_F1) = ", fmt[zetaPhysF1Numeric]];
expectApprox["zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)", zetaPhysF1Numeric, zetaMaxF1, 10^-10];
```

**Assessment:**
Edit matches the directive verbatim (including the `FindRoot` seed `y = 153/100`, `WorkingPrecision -> 40`, and tolerance `10^-10`). The previously-dead `zetaPhys`, `kappaF1`, `etaF1` are now exercised: the limit computes `omegaPe^2 -> 1` as `Pe -> Infinity`, giving `zetaPhys -> (kappa + Pi^2/4)/(kappa + y^2)` evaluated at `(kappa = 12321/5, y = root of y*tan(y) = 37)`. Output line 19 prints `zeta_phys(Pe->oo, ...) = 2.46752922945601223332958450157053542039...` and line 20 reports `diff = 2.23e-15` (well under the `10^-10` tolerance), followed by `PASS:` on line 21. The Mathematica numeric agrees with the upstream `zetaMaxF1 = 2.46752922945601` to roughly 14 significant digits — well past the BLOCKED-IF threshold of 10 digits — so no block was warranted. The `Limit::alimv` warning about assumptions involving the limit variable is benign here: the limit still evaluates symbolically to the closed-form value.

## Exec log assessment

**SymPy:** exit=n/a. No SymPy script exists for unit 084 (status-only carve-out per the manifest, as established by the original audit).

**Mathematica:** exit=0 (inferred — no captured `.log` file exists; orchestrator captured only the diff). The canonical output `mathematica/output/moving_throat_pde_stage084_full_reduced_pde_writeup_mathematica_audit.txt` ends with `Stage 084 Mathematica audit passed.` and contains 7 PASS lines and 0 FAIL lines:
- `PASS: inverse demand map`
- `PASS: Xi_F1(Upsilon|Upsilon_w->100 Theta_w) - Xi_F1(Theta)`
- `PASS: zeta_phys at Family-1 (Pe->oo limit) matches upstream zeta_max^(F1)` (diff = 2.23e-15)
- `PASS: chi-window ordering positive (zeta_+^chi > zeta_-^chi)`
- `PASS: hard-ceiling gap positive (zeta_max^F1 > zeta_+^chi)`
- `PASS: J-window ordering positive (zeta_+^J > zeta_-^J)`
- `PASS: fail-side J below chi (zeta_+^J <= zeta_+^chi)`

One warning surfaces in the run: `Limit::alimv: Warning: Assumptions that involve the limit variable are ignored.` This is benign — the limit still evaluates to the expected closed-form value and the F3 check passes with sub-femto-digit residual.

**Output freshness:** Script mtime 1779690721; output mtime 1779690790. Output is 69 seconds newer than the script — confirmed re-generated post-fix.

## Material-change assessment

`material_change`: false.

The edits change *the form of assertions* without altering any derived symbolic content or carried-forward numeric value. The printed `zeta_phys`, `zeta_req`, `Q`, `R_quad`, `kappa_F1`, `eta_F1`, and the five carried `zeta` literals are unchanged. The F3 insertion adds a new printed value `zeta_phys(Pe->oo, kappa_F1, eta_F1, y_F1) = 2.46752922945601...`, but that value is a *consistency check* against the already-imported `zetaMaxF1` from stage 081 — it does not introduce a new constant that downstream units would consume. No downstream unit depends on outputs unique to 084.

## Side observations (non-blocking)

- The banner on line 32 still reads `"STAGE 067 — FULL REDUCED MOVING-THROAT PDE WRITE-UP SKELETON"`, propagating the original cosmetic mislabel into the output (line 3 of the .txt). The auditor explicitly excluded this from findings per the doc-alignment exclusion; just noting it persists.
- The `Limit::alimv` warning in the captured output is harmless for this case but could theoretically obscure a real issue in future edits — Mathematica is ignoring `pe > 0` because `pe` is the limit variable. Since the limit still evaluates correctly to a 20-digit numeric matching the upstream value to 14+ digits, no action needed here.

## Verdict justification

All three findings resolved verbatim per the directive: F1 replaces tautological literal-difference checks with a substantive cross-route consistency check that exercises the shared `lambdaEll` between the two `Xi_F1` formulas; F2 swaps hardcoded float-subtraction `expectApprox` calls for four genuine ordering inequalities that also bring the previously-unused `zetaMinusJ`/`zetaPlusJ` into the assertion set; F3 introduces a new `expectApprox` pinning the script's own `zetaPhys` symbolic form (evaluated at `Pe -> Infinity`, `kappa = 12321/5`, `eta = 37`, `y = y_F1`) to the upstream-imported `zetaMaxF1`, agreeing to ~14 significant digits. The output is fresh, the script ends with the success banner, and 7/7 PASS lines appear with no failures. No regressions in the diff. Material change is false because no symbolic content or carried-forward numeric value was altered.

stage 084: verified
