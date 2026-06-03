---
unit_id: 245
batch: VIII.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-06-03T10:20:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 5
findings_total: 5
material_change: false
---

# Verification — unit 245

## Per-finding outcomes

### F1 — insufficient_verification (vacuous support-blindness self-test)

**Classification:** resolved

**What changed:**
`scripts/...stage245...sympy_audit.py:195-203` — the pre-existing zero-asserts on the seven absent-variable objects (lines 187-193) are kept, and a positive control is appended: `U_bad = a_V*(f_U+Lam)/Delta_UV`, `dU_bad_dLam = diff(U_bad, Lam)`, then `assert dU_bad_dLam != 0` and `assert simplify(dU_bad_dLam - a_V/Delta_UV) == 0`. Mirrored in `.wl:206-215` (M7) via `expectTrue[controlDerivative =!= 0]` + `expectZero[controlDerivative - aV/(aU aV - chiUV^2)]`.

**Assessment:**
Genuinely non-vacuous now. The control exercises the differentiation route on an expression that *does* contain `Lam` and pins the exact nonzero result `a_V/Delta_UV`, demonstrating the check has discriminating power; the seven zero-derivative checks are thereby meaningful rather than guaranteed-zero. The vacuous-only situation flagged in the report is gone. The F1 trap is NOT re-introduced in the `.wl`.

### F2 — insufficient_verification (R_target not derived from selected-branch identity)

**Classification:** resolved

**What changed:**
`.py:106-111,123` — introduces `Lambda_0` (positive), defines `R_target_from_id = Lambda_0*(1-eps_eta)/T2` and `R_target_ref_from_id = Lambda_0*(1-eps_eta_ref)/T2_ref`, forms `R_ratio_derived`, and asserts `simplify(R_ratio_derived - R_ratio_paper) == 0`. Mirrored in `.wl:136-149` (M4).

**Assessment:**
This is a genuine implication from the selected-branch identity `R_target*T^2 = Lambda_0(1-eps_eta)`. `R_ratio_derived` is built independently from the identity (NOT from `R_ratio`), so it is not the old `1-1=0` round-trip. The old self-referential `R_exact_check` and the weak A10/A11 are retained but are no longer load-bearing. Lambda_0 cancels in the ratio as expected. Correct.

### F3 — tautological_check (Session-I rebuild was an inverse round-trip)

**Classification:** resolved

**What changed:**
`.py:242-245` — the two `< 1e-12` round-trip asserts on `U_rebuilt`/`V_rebuilt` are removed and replaced with `abs(float(R_ratio_obs) - 0.87984149) < 5e-9` and `abs(float(R1_obs) - (-0.12762119)) < 5e-9`. The `eps_eta_obs` anchor (line 241) is unchanged. Mirrored in `.wl:230-232` (M8).

**Assessment:**
The two genuine session anchors (0.87984149, -0.12762119) are now ASSERTED. Critically, `R_ratio_obs` (line 217) and `R1_obs` (line 227) are computed directly from the observed `(U_obs, V_obs)` via the compiler formulas, NOT from reconstructed parameters — so these are honest cross-checks, not `x = g(g^-1(x))`. The `_rebuilt` quantities survive as prints only. Correct.

### F4 — insufficient_verification (drain nonnegativity not asserted)

**Classification:** resolved

**What changed:**
`.py:90-95` — concrete opposite-sign admissible point (`a_U=2.5, a_V=3.0, chi_UV=-0.76, f_U=0.33`, giving `Delta=6.9224>0`) with `assert float(drain_neg_chi) > 0`; value ≈ 0.003938. Mirrored in `.wl:125-132` (M3) using exact rationals.

**Assessment:**
Genuinely exercises the "nonnegative even when U and V have opposite signs" claim by sampling the `chi_UV < 0` branch (the Session-I sign pattern). Non-vacuous numeric probe on an admissible point. Correct.

### F5 — missing_verification_script (missing_mathematica)

**Classification:** resolved

**What changed:**
New file `mathematica/...stage245...mathematica_audit.wl` (238 lines) covering M1–M8, exits 0 with all PASS.

**Assessment:**
Genuinely INDEPENDENT, not a transliteration: native `Solve` for the stationary packet, `Det` for the Hessian, native `Series`/`Normal`/`Coefficient` for the first-order packet (vs the `.py`'s custom `coeff_linear`), `R_target` derived from the selected-branch identity (the route the `.py` originally skipped), positive-control support-blindness (M7), and structural drain nonnegativity (M3). Session-I numbers (M8) asserted directly from `uObs`/`vObs` with no parameter reconstruction — so the `.wl` does NOT reproduce the F3 defect, and M7 does NOT reproduce the F1 trap. Uses the project `expectZero`/`expectTrue`/`expectApprox` idioms and strips `ConditionalExpression`.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `D_UV at opposite-sign point = 0.003937884170811954` (F4 probe positive)
- `R_target/R_ref (from identity) = (eps_eta_ref*exp(...) - 1)*exp(-...)/(eps_eta_ref - 1)` matches the paper form (F2 derivation)
- `control d/dLam U_bad = a_V/(a_U*a_V - chi_UV**2)` (F1 control fires nonzero)
- `R_target/R_ref(obs) = 0.8798414919352429`, `R1_obs = -0.1276211900000000` (F3 anchors)
- `All symbolic checks passed.`

**Mathematica:** exit=0, 38 PASS lines, 0 FAIL. Notable lines:
- `M3 opposite-sign drain > 0 = True` → PASS (F4)
- `M4 R_target ratio from identity = 0` → PASS (F2 route)
- `M7 positive control detects support contamination = True` and `M7 positive control exact derivative = 0` → PASS (F1)
- `M8 R_target/R_ref Session-I` diff ≈ 1.94e-9 ≤ 5e-9 → PASS; `M8 R1 Session-I` diff = 0 → PASS (F3)
- `Stage 245 Mathematica audit passed.`
- Note: the M1 `Solve` output carries a `ConditionalExpression[...]` wrapper, but it is stripped by `cleanScalar`/`stripConditional` before the `expectZero` comparisons, which then read 0 — no spurious pass.

**Output freshness:** confirmed. SymPy output `.txt` mtime 1780502817 > script mtime 1780502400; Mathematica output `.txt` mtime 1780502817 > `.wl` mtime 1780502497. Both regenerated post-fix.

## Material-change assessment

`material_change`: false. No script *values* changed — every derived result (U_sol, V_sol, drain form, R_target ratio, first-order coefficients, Session-I numbers) is identical to pre-fix. The edits only widen the verification surface (added positive control, identity-based R_target derivation, drain-positivity probe, session-anchor asserts, and a second engine). No downstream unit depends on a changed number.

## Side observations (non-blocking)

- `.py` retains the now-redundant `R_ratio`, `R_exact_check`, and the weak A10/A11 tautological asserts. They are harmless (the directive explicitly allowed keeping them) and no longer load-bearing. Not a defect.
- `.py:102` declares `R_target_ref` as an unused symbol; cosmetic only.

## Verdict justification

All five findings are resolved with genuine, non-tautological checks confirmed by reading the code (not by trusting the passing logs): F1's positive control gives the support-blindness check real discriminating power; F2 derives R_target/R_ref from the selected-branch identity independently of the old round-trip; F3 asserts the two genuine Session-I anchors computed directly from the observed point rather than reconstructed parameters; F4 probes drain positivity on the opposite-sign branch; and F5 adds a genuinely independent Mathematica route (38 PASS, exit 0) that uses a different decomposition and reproduces neither the F1 nor the F3 defect. Both engines exit 0 with fresh outputs and no regressions in the diff. material_change is false (verification surface widened, no values touched).
