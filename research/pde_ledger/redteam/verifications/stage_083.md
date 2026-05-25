---
unit_id: 083
batch: III.4
verifier_model: claude-opus-4-7-1m
verify_date: 2026-05-25T00:35:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 4
findings_total: 4
material_change: false
---

# Verification — unit 083

## Per-finding outcomes

### F1 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.wl:62-91` — added `delta0Residual` / `deltaInfResidual` blocks each closed by an `expectZero` verifying that the closed forms for `Delta_0(F1)` and `Delta_inf(F1)` satisfy their defining linear identities `(alpha^2 (alpha sinh + eta cosh)) Delta_0 - eta (cosh - 1) == 0` and `(alpha sinh + eta cosh) Delta_inf - (cosh + (eta/alpha) sinh - 1) == 0`.
- `mathematica/...:96-103` — added `yRootResidual` and `aF1Indep` blocks closed by `expectApprox` checks for `y_F1*Tan[y_F1] - eta == 0` and `aF1 - aF1Indep == 0`.
- `mathematica/...:137-154` — added the `ClearAll[pTest]` shim plus the `omegaResidual` block closed by `expectZero` for the Omega-identity `Omega(Pe)(4 Pe^2+pi^2)(e^Pe-1) - pi Pe (2 Pe e^Pe + pi) == 0`.

**Assessment:**
Mathematica output (lines 8-15, 26-27) shows all four new identity checks PASS with residual `0` (Delta_0 / Delta_inf / Omega) or residual `0``46.85` and `0``39.70` (precision-marker zero for the two numeric `expectApprox` lines). The Delta and Omega defining-equation identities are non-tautological in the operational sense the directive specified: they verify that the closed forms typed in are the unique solutions of the linear equations that originally defined them. A copy-time typo in `delta0F1`, `deltaInfF1`, or `omega[pp_]` would now drop the residual to a non-zero expression and trip `expectZero`. The `aF1Indep` check is structurally trivial (same expression with `pi^2/4` rewritten as `(Pi/2)^2`), which the directive itself anticipates and accepts as a "low-cost" cross-evaluation — flagged as a side observation, not a verification failure.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.py:97-105` — inserted the `# SOURCE-ANCHOR (operator selectors): ...` block immediately above the `Theta_chi_coeff`, `Theta_J_coeff`, and `136900` literals.
- `mathematica/...:110-118` — inserted the corresponding `(* SOURCE-ANCHOR ... *)` block immediately above `thetaChi`, `thetaJ`, and `xiChi/xiJ` literals.

**Assessment:**
Both comment blocks are present at the directive-specified locations, name the literals, declare no upstream-script anchor exists in this unit, and leave a `TODO(upstream-anchor)` placeholder. The directive explicitly forbids inventing an equation; Codex did not. Documentation-of-trust-boundary fix exactly as requested.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/.../sympy_audit.py:66-81` — added Delta_0 / Delta_inf defining-equation `expect_zero` residuals on the SymPy side mirroring F1's Mathematica additions.
- `scripts/.../sympy_audit.py:84-89` — bumped `nsolve` to `prec=80` (deviation flagged in the Applied block) and added the `y_F1*tan(y_F1) - eta` residual with a hard `> 1e-25` AssertionError gate.
- `scripts/.../sympy_audit.py:137-150` — added the `dzeta_dPe` sample loop at `Pe in {10, 100, 1000, 10000}` with a per-sample `<= 0` AssertionError.
- `mathematica/...:181-194` — added the `ClearAll[ppSym, dzetaDpe]` shim, `dzetaDpe[pp_]` definition, and `Module[{vals, signs}, ...]` sampling/sign-check block (uses the existing `pass/fail` helpers).

**Assessment:**
SymPy output lines 12-14 show `Delta_0(F1) defining-equation residual = 0`, `Delta_inf(F1) defining-equation residual = 0`, and `y_F1 defining-equation residual = -3.37e-79` (well below `1e-25`). Lines 25-28 show the four monotonicity sample points all positive (`2.24e-2`, `2.43e-5`, `2.44e-8`, `2.44e-11`), so the AssertionError gate was not triggered. Mathematica output line 55 shows the four `d zeta / d Pe` samples all positive and line 56 prints `PASS: zeta_F1 monotone increasing`. The monotonicity check is non-tautological: it exercises `d/dPe (A_F1 * Omega(Pe)^2)`, a derivative the script does not separately know is positive, and an algebra slip in `Omega` (e.g., flipped sign or wrong exponent) would surface as a negative derivative sample. The `prec=80` deviation is benign — it tightens nsolve precision without changing any downstream value (still prints `1.52948248371469964992710762240`, agreeing with the Mathematica `FindRoot` output to ~30 digits).

### F4 — tautological_check

**Classification:** resolved

**What changed:**
- No edit to lines 169-179 of the Mathematica script. Codex selected directive option (b): keep the existing eleven `expectApprox` numeric self-checks as cheap sanity guards, relying on F1's new `expectZero` identity checks for substantive verification.

**Assessment:**
The directive explicitly authorises this choice ("Verifier will accept either: (a) lines 103-113 deleted... or (b) lines 103-113 unchanged but F1's `expectZero` lines added and PASS"). F1's identity checks are present and PASS (verified above). Choice (b) is the safer mechanical option per the directive's recommendation. No regression introduced.

## Exec log assessment

**SymPy:** exit log file missing (`redteam/exec_logs/stage_083_sympy.log` not present — only the diff was captured). The canonical output `scripts/output/moving_throat_pde_stage083_family1_direct_operator_window_sympy_audit.txt` (mtime 00:29, post-fix) is fresh and contains all four newly-required print lines:
- Line 12: `Delta_0(F1) defining-equation residual = 0`
- Line 13: `Delta_inf(F1) defining-equation residual = 0`
- Line 14: `y_F1 defining-equation residual = -3.373...e-79`
- Lines 25-28: four `d zeta / d Pe at Pe = ...` samples, all positive.

The script's only failure modes (AssertionError in `expect_zero` and the y/monotonicity gates) would have left a truncated output and a non-zero exit; the output is complete and contains every expected line, so exit was 0. Treating sympy_exit as 0 with the caveat that the log was not preserved.

**Mathematica:** exit log file missing (`redteam/exec_logs/stage_083_mathematica.log` not present). The canonical output `mathematica/output/moving_throat_pde_stage083_family1_direct_operator_window_mathematica_audit.txt` (mtime 00:30, post-fix) is fresh and contains every required PASS:
- Line 9: `PASS: Delta_0(F1) defining-equation residual`
- Line 11: `PASS: Delta_inf(F1) defining-equation residual`
- Line 13: `PASS: y_F1 satisfies y Tan[y] = eta`
- Line 15: `PASS: A_F1 independent vs closed-form`
- Line 27: `PASS: Omega(Pe) identity residual`
- Lines 33-54: all eleven legacy `expectApprox` numeric sanity checks PASS (diffs in `10^-11` to `10^-19` range).
- Line 56: `PASS: zeta_F1 monotone increasing`.

The `fail[]` helper exits with code 1; absence of `FAIL:` lines anywhere in the output plus presence of the terminal `Pi_max^(F1)/C_mix` line confirms the script ran to completion. Treating mathematica_exit as 0 with the same log-missing caveat.

**Output freshness:** confirmed. Script mtimes are 00:26-00:27; output mtimes are 00:29-00:30. Outputs are newer than scripts. The freshness invariant the orchestrator requires is satisfied.

## Material-change assessment

`material_change`: false.

All numeric printed values are unchanged from the pre-fix outputs:
- `Delta_0(F1) = 0.000173302079021525149057156196550` (same).
- `Delta_inf(F1) = 0.0201447565540521594271032956099` (same).
- `y_F1 = 1.52948248371469964992710762240` — SymPy gained a few trailing digits relative to the pre-fix `nsolve` default precision, but the leading 16 digits agree with the Mathematica `FindRoot` value, which itself is unchanged. No downstream Pe / zeta / Pi value is altered to displayed precision.
- All `Pe_{+/-}^{chi/J}`, `zeta_*`, `Pi_*/C_mix` values unchanged.

The only new content is identity residuals (all `0`), monotonicity sample derivatives (all positive), and trust-boundary comments. No derived result that a downstream unit could depend on has changed.

## Side observations (non-blocking)

1. The F1 `aF1Indep` check rewrites `pi^2/4` as `(Pi/2)^2`. Mathematically these are identical, so the `aF1 - aF1Indep` residual is `0` by definition rather than because two derivation routes agreed. This was anticipated in the directive ("The acceptable alternative is to delete line 103. Codex picks ONE of: (a) delete, (b) keep as numeric sanity. Either is acceptable.") so it is not a verification failure, but the check carries no information that the closed-form line above it does not already carry. Flagging for the auditor's awareness on a possible future pass.
2. The F1 `delta0Residual` / `deltaInfResidual` identities are exactly the rearrangements of the defining ratio used to construct `delta0F1` / `deltaInfF1` — they are guaranteed to FullSimplify to zero. The directive's self-test note (line 299 of `directives/stage_083.md`) acknowledges this and argues the check still catches *transcription* errors (e.g., missing a factor of alpha when typing the closed form). That reasoning is correct: a transcription typo in either the residual identity or the closed form would diverge, surfacing the bug. Recording for the record; not a verification issue.
3. The exec log files (`stage_083_sympy.log`, `stage_083_mathematica.log`) were not captured under `redteam/exec_logs/`. Only the diff patch was preserved. The orchestrator's `fix_loop` apparently regenerated outputs but did not tee the runtime logs. The canonical `.txt` outputs contain all the assertion/PASS evidence the verifier needs, so this is non-blocking, but flagging in case the orchestrator wants to tighten log retention.

## Verdict justification

All four findings are `resolved`. F1 added three independent identity checks on the Mathematica side; F2 added the requested SOURCE-ANCHOR comments to both engines; F3 added Delta defining-equation residuals plus y-root and monotonicity gates in SymPy, mirrored by a monotonicity sign-check in Mathematica; F4 was a clean-up that the directive permitted to be deferred under option (b) once F1 landed, which Codex did with proper justification. Both canonical outputs are fresh and contain every required PASS line with no FAIL lines anywhere. No printed numerical result changed, so `material_change: false` and no downstream stages are stale.

stage 083: verified
