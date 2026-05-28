---
unit_id: 143
batch: IV.5
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-27T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 143

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- SymPy (`scripts/moving_throat_pde_stage143_equal_normalized_singular_limit_sympy_audit.py`):
  - Lines 19-23: added `expect_equal(name, lhs, rhs)` helper.
  - Lines 25-29: added `expect_positive(name, expr)` helper.
  - Lines 50-51: `expect_positive("pi**2 - 2*pi > 0", ...)` and `expect_positive("pi**2/2 - 4 > 0", ...)`.
  - Lines 53-54: Taylor-series check `expect_equal("exp remainder leading term is Pi**3/6", exp_rem_series.coeff(Pi, 3), sp.Rational(1, 6))`.
  - Lines 61-62: `expect_equal("lim_{Pi->0+} g_Pi == 2/pi", g0, 2/pi)` and `expect_equal("lim_{Pi->oo} g_Pi == 1", ginf, sp.Integer(1))`.
  - Lines 77-79: `expect_equal` calls for `R_infty == (1-r)**2/(1+r**2)`, `S_infty == 1`, and `lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))`.
- Mathematica (`mathematica/moving_throat_pde_stage143_equal_normalized_singular_limit_mathematica_audit.wl`):
  - Lines 33-37: added `expectEqual` helper.
  - Lines 39-43: added `expectPositive` helper.
  - Lines 71-72: `expectPositive["Pi^2 - 2*Pi > 0", ...]` and `expectPositive["Pi^2/2 - 4 > 0", ...]`.
  - Lines 73-74: `expRemSeries = Normal[Series[...]]` followed by `expectEqual["exp remainder leading term is piM^3/6", ..., 1/6]`.
  - Lines 82-83: `expectEqual` calls for the two endpoint limits.
  - Lines 100-102: `expectEqual` calls for `R_infty`, `S_infty`, and `tHatRatio`.

**Assessment:**
All eight new assertions specified in the directive are present on both engines and fire in the exec logs. The endpoint-limit assertions are non-tautological (they gate the result of `sp.limit` / `Limit[]`, not a hand-set value). The positivity assertions on `pi**2 - 2*pi` and `pi**2/2 - 4` are non-tautological because `sp.simplify(pi**2 - 2*pi)` returns the symbolic value `pi*(-2+pi)` and the `.is_positive` check uses sympy's numerical knowledge of `pi`; same on the Mathematica side via `Simplify[val > 0]` (these are concrete numeric inequalities about π that a typo in the constant would break). The Taylor-coefficient check `coeff(Pi, 3) == 1/6` non-trivially gates the exp series. The `S_infty == 1` assertion on the Mathematica side would have been tautological under the old `sInf = 1;` definition; under F2's fix it now gates the limit of the full `sQ` formula. No collateral edits beyond what the directive specified.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
Mathematica lines 85-99 replace the previous algebraic-substitution block. Now:
- Line 85: `sQ` defined (same closed form, but now actually used in a limit).
- Line 86: `rQ = (gPi - r)^2/(1 + r^2)` defined as a dynamical object in `piM`.
- Line 87: `sigma0 = piM/(1 - rQ*sQ)`.
- Line 88: `that = Sqrt[(9/20)*sigma0]`.
- Lines 90-94: `Clear[piInf2, piInf3, piInf4, piInf5]` then `rInf`, `sInf`, `sigmaRatio`, `tHatRatio` all computed via `Limit[... /. piM -> piInfN, piInfN -> Infinity]` with positivity assumptions.

The hardcoded `sInf = 1;` is gone.

**Assessment:**
The replacement matches the directive's specification line-for-line in structure and intent. The exec log confirms the new path produces the expected closed-form symbolic outputs: `R_infty = 1 - (20*Pi*Sqrt[4107 - 100*Pi^2])/4107` (equal to SymPy's `(-sqrt(4107 - 100*pi**2) + 10*pi)**2/4107`; the auditor already showed these are algebraically identical), `S_infty = 1`, `lim Sigma0/Pi = 4107/(20*Pi*Sqrt[4107 - 100*Pi^2])`, `lim That/sqrt(Pi) = (111*Sqrt[3/Pi])/(20*(4107 - 100*Pi^2)^(1/4))`. Numerics still 30-digit-agree with SymPy. F1's three downstream assertions on these objects now gate genuine `Limit` outputs rather than hardcoded literals.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Mathematica lines 60-65 (immediately after the existing `expectZero["numerator - exact positive decomposition", num - decomp];` on line 58) add the independent positivity check via `Reduce`:
```
numPositiveCheck = Reduce[num > 0, piM, Reals] /. {(Element[piM, Reals] && piM > 0) -> True, (piM > 0) -> True};
Print["Reduce[num > 0, piM, Reals] = ", fmt[numPositiveCheck]];
If[TrueQ[Simplify[numPositiveCheck === True || numPositiveCheck === (piM > 0)]],
  pass["num > 0 for piM > 0 via Reduce"],
  fail["num > 0 for piM > 0 via Reduce", numPositiveCheck]];
```

**Assessment:**
The block matches the directive verbatim. Exec log line 12-13 confirms `Reduce[num > 0, piM, Reals] = True` followed by `PASS: num > 0 for piM > 0 via Reduce`. This is a structurally independent positivity certification: `Reduce` over the reals does not subtract a hand-rolled decomposition, so it provides genuine second-engine cross-corroboration of the inequality $\mathfrak g_\Pi<1$ via a different mechanism than the SymPy line-by-line algebraic decomposition. Non-tautological — a typo in `num` or `gPi` would cause `Reduce` to either fail or return a different constraint, which the `If[TrueQ[...]]` guard would reject.

## Exec log assessment

**SymPy:** exit=0 (inferred: log ends cleanly with the "STAGE 143 LEDGER" banner and 30-digit numerics; no AssertionError, no traceback). Notable lines:
- `numerator - exact positive decomposition = 0` (A1 identity check passes).
- `exp remainder leading term is Pi**3/6: lhs - rhs = 0`.
- `lim_{Pi->0+} g_Pi == 2/pi: lhs - rhs = 0` and `lim_{Pi->oo} g_Pi == 1: lhs - rhs = 0`.
- `R_infty == (1-r)**2/(1+r**2): lhs - rhs = 0`, `S_infty == 1: lhs - rhs = 0`, `lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty)): lhs - rhs = 0`.
- `R_infty = 0.145454452260420126101421595368`, `lim That/sqrt(Pi) = 0.725669130700713219781041125011`.

**Mathematica:** exit=0 (confirmed by final lines `Stage 143 Mathematica audit passed.` and `Exit[0]` in the script body). Notable lines:
- `PASS: numerator - exact positive decomposition`.
- `Reduce[num > 0, piM, Reals] = True` + `PASS: num > 0 for piM > 0 via Reduce`.
- `PASS: Pi^2 - 2*Pi > 0`, `PASS: Pi^2/2 - 4 > 0`, `PASS: exp remainder leading term is piM^3/6`.
- `PASS: lim_{piM->0+} g_Pi == 2/Pi`, `PASS: lim_{piM->oo} g_Pi == 1`.
- `PASS: R_infty == (1-r)^2/(1+r^2)`, `PASS: S_infty == 1`, `PASS: lim That/sqrt(Pi) == sqrt((9/20)/(1-R_infty))`.
- 30-digit numerics `0.1454544522604201261014215953679127185274370895574784984548` and `0.7256691307007132197810411250114145098613931689620003558428`.

**Output freshness:** confirmed. mtimes:
- `scripts/.../sympy_audit.py`     2026-05-27 19:50:00
- `scripts/output/...sympy_audit.txt`  2026-05-27 19:51:40 (newer)
- `mathematica/.../mathematica_audit.wl` 2026-05-27 19:50:40
- `mathematica/output/...mathematica_audit.txt` 2026-05-27 19:53:33 (newer)

Both exec logs were regenerated after the corresponding script edits.

## Material-change assessment

`material_change`: false.

No derived quantity changed: `R_infty`, `S_infty`, `lim Sigma0/Pi`, `lim That/sqrt(Pi)` all retain the same closed forms and 30-digit numerics as before the fix. F2 replaced the path by which the Mathematica engine reaches `S_infty = 1` (hardcoded → limit) but the value is the same. F3 added a new positivity check without changing any printed value. F1 added assertions that gate already-printed quantities; no new outputs.

## Side observations (non-blocking)

- Auditor noted (self-test "cosmetic" remark) that both scripts previously had a stage-126 banner. The orchestrator notes say "Banner already handled by Cluster A. Both engines now read STAGE 143." Current scripts confirm: SymPy banner reads "STAGE 143 — EQUAL-NORMALIZED BRANCH IS A SINGULAR LIMIT" (line 34) and Mathematica reads the same (line 45). Non-blocking — already fixed.
- The directive's front-matter still has `applied: false` and there are no `## Applied: F1/F2/F3` blocks appended to the directive. Despite the missing bookkeeping, all three findings are materially addressed in the code and exec logs. Non-blocking for verification verdict, but the orchestrator may want to flip the front-matter flag.

## Verdict justification

All three findings are materially resolved: F1 added 8+ non-tautological assertions on both engines that exercise every paper deliverable (strict inequality, endpoint limits, R_infty, S_infty, traction divergence ratio); F2 replaced the hardcoded `sInf = 1` with a `Limit[sQ ...]` and constructed the dynamical `rQ`, `sigma0`, `that` objects so the Mathematica side now traces the same derivation chain as SymPy; F3 added an independent `Reduce[num > 0, piM, Reals]` positivity check that satisfies the second-engine policy. Both exec logs are fresh, exit 0, and show every directive-specified PASS line. No derived quantity changed, so no downstream units are stale.
