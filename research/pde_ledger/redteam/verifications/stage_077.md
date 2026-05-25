---
unit_id: 077
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 077

## Per-finding outcomes

### F1 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:45-47` — three new lines inserted immediately after the `xi_*` definition/print block: `S_at_star = (1 + sp.tanh(xi_star)) / 2`, `rho_quartic_at_star = sp.simplify(1 - alpha_r * S_at_star**2)`, `expect_zero("1 - alpha_r * S(xi_*)**2", rho_quartic_at_star)`.
- `mathematica/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.wl:53-55` — three new lines inserted between `Print["xi_*(alpha_r=10) = ", ...]` and `banner["NUMERICAL FAMILY-1 EXTRACTION"]`: `sAtStar = (1 + Tanh[xiStar])/2`, `rhoQuarticAtStar = FullSimplify[1 - alphaR*sAtStar^2, Assumptions -> $Assumptions]`, `expectZero["1 - alphaR*S[xi_*]^2", rhoQuarticAtStar]`.

**Assessment:**
The edits match the directive verbatim (variable names, order, expression form). The new SymPy assertion goes through the existing `expect_zero` infrastructure (which calls `sp.simplify(sp.expand(...))` and raises `AssertionError` if non-zero), so a regression in the closed form for `xi_*` would now fail symbolically. The Mathematica assertion routes through `expectZero` which uses `FullSimplify[Together[Expand[...]]]` under `$Assumptions` (already `Element[{xi, alphaR}, Reals] && alphaR > 0`), and only `pass[name]` on `res === 0`. The SymPy output (line 10) shows `1 - alpha_r * S(xi_*)**2 = 0` with no AssertionError; the Mathematica output (lines 11-12) shows `1 - alphaR*S[xi_*]^2 = 0` followed by `PASS: 1 - alphaR*S[xi_*]^2`. Non-tautological: the assertion reduces a nested `1 - alpha_r * ((1+tanh(atanh(2/sqrt(alpha_r)-1)))/2)^2` expression through the engines' simplifiers, exercising real algebra rather than echoing a literal. No collateral edits in either file.

### F2 — tautological_check

**Classification:** resolved

**What changed:**
The diff shows the removal of `expectApprox["xi_* numeric check", xiCut, ToExpression["-0.38558106921542562403635498846713378847348301441599`50"], 10^-30];` from the `.wl` file. Lines 92-95 of the current file retain the four legitimate cross-engine `expectApprox` checks for `<rho>_chi`, `<rho^2>_chi`, `Theta_w^(chi)`, `Theta_w^(J)`. Line 85 (`Print["numeric xi_* = ", fmt[xiCut]]`) is preserved, keeping the value visible.

**Assessment:**
Exactly the line specified in the directive was deleted and nothing else. The Mathematica output (lines 14-32) shows no `xi_* numeric check diff` or `PASS: xi_* numeric check` strings, while the four remaining `expectApprox` PASS lines and the Jensen `expectTrue` PASS line are present. Combined with F1's new symbolic identity, the cut-point claim is now backed by a real check rather than a number-vs-self comparison.

### F3 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.py:82-111` — added `def expect_close(name, value, target, tol)` helper followed by four `expect_close` calls against the same 50-digit constants used by the Mathematica `expectApprox` checks: `<rho>_chi` (tol 1e-28), `<rho^2>_chi` (tol 1e-28), `Theta_w^(chi)` (tol 1e-26), `Theta_w^(J)` (tol 1e-27). Inserted between the two `print` lines for the Theta values (line 80) and the Jensen-ordering guard (line 112). Jensen ordering preserved exactly as before.

**Assessment:**
The block matches the directive verbatim (helper signature, constants, tolerances, order). The four constants are the cross-engine ground truth from the Mathematica targets on lines 92-95 of the `.wl`, not invented by SymPy, so the assertions are non-tautological. The SymPy output (lines 17-20) confirms the four diffs are all roughly 1e-50 to 1e-51, well below their respective tolerances, indicating the script ran cleanly (exit 0 implied — no AssertionError raised, no traceback in the output, and the Jensen check is downstream of the new asserts). The Jensen-ordering inequality on lines 112-113 is unchanged, satisfying the "keep intact" directive clause.

## Exec log assessment

**SymPy:** exit log file is absent (`stage_077_sympy.log` does not exist in `redteam/exec_logs/`). However, the canonical output file `scripts/output/moving_throat_pde_stage077_family1_theta_extraction_sympy_audit.txt` exists, is fresh (mtime 2026-05-22 23:26:30, newer than script mtime 23:24:45), and shows all expected new lines plus no truncation/traceback markers. Notable lines:
- `I_f - 1/3 = 0` (line 7) — symbolic identity preserved.
- `1 - alpha_r * S(xi_*)**2 = 0` (line 10) — F1 SymPy check passing.
- `<rho>_chi diff = 4.0091...e-51` through `Theta_w^(J) diff = 2.67...e-51` (lines 17-20) — F3 SymPy checks all well within tolerance.
The orchestrator-issued `redteam exec-sympy 077` run is the source of those outputs; absence of the explicit log file is a logging gap, not a script failure. Treating `sympy_exit: 0` based on the canonical .txt being regenerated successfully with all expected content.

**Mathematica:** exit log file is also absent (`stage_077_mathematica.log` does not exist). The canonical output `mathematica/output/moving_throat_pde_stage077_family1_theta_extraction_mathematica_audit.txt` is fresh (mtime 2026-05-22 23:26:43, newer than script mtime 23:24:45). Notable lines:
- `PASS: I_f - 1/3` (line 8)
- `PASS: 1 - alphaR*S[xi_*]^2` (line 12) — F1 Mathematica check passing.
- `PASS: <rho>_chi numeric check` ... `PASS: Theta_w^(J) numeric check` (lines 24, 26, 28, 30) — four legitimate cross-engine checks.
- `PASS: Theta_w^(chi) >= Theta_w^(J) > 0` (line 32) — Jensen ordering.
- `Stage 077 Mathematica audit passed.` (line 34) followed by `Exit[0]` in the script.
No `FAIL:` strings, no `xi_* numeric check` strings (F2 cleanly excised). Treating `mathematica_exit: 0`.

**Output freshness:** confirmed. Both `.txt` outputs (mtimes 23:26:30 and 23:26:43) are newer than their corresponding source files (mtimes 23:24:45 and 23:24:45). The orchestrator's fix_loop refreshed them after Codex's edits.

## Material-change assessment

`material_change`: false.

No printed numeric value or symbolic form changed. The numeric moments, Theta values, denominator, cut point, and `I_f` print identically before and after the diff (the diff only adds new assertion lines and removes one tautological assertion line; it does not alter any computation path or output value that downstream units consume). The new SymPy `expect_close` helper is local to the file and does not affect exported state. The new symbolic identity check in both engines reduces to `0` rather than substituting into any downstream expression. Downstream units that depend on `<rho>_chi`, `<rho^2>_chi`, `Theta_w^(chi)`, `Theta_w^(J)`, or `xi_*` see the same values they saw pre-fix.

## Side observations (non-blocking)

- The `redteam/exec_logs/` directory lacks both `stage_077_sympy.log` and `stage_077_mathematica.log`, even though the canonical output `.txt` files were refreshed and the diff patch was captured. If the orchestrator's logging path was supposed to write `*_sympy.log` and `*_mathematica.log` alongside `*_diff.patch`, that emission appears to have been skipped for unit 077. This is a tooling/logging gap, not a script defect — flagging so the orchestrator owner can decide whether to backfill.
- The SymPy `expect_close` helper is defined locally inside the script body (line 82) rather than next to `expect_zero` at the top (line 24). Stylistically inconsistent but functionally fine; the directive specified this placement, so it's a directive-prescribed choice, not a Codex deviation.

## Verdict justification

All three findings are resolved with edits matching the directive verbatim (no collateral changes, no deviations). The refreshed `.txt` outputs confirm: F1's new symbolic vanishing identity prints `0` in both engines and passes (PASS line on the Mathematica side), F2's tautological numeric self-check is excised from the Mathematica output, and F3's four new SymPy `expect_close` assertions all report sub-1e-50 diffs against the cross-engine 50-digit targets. No regressions visible in the diff. Both scripts have fresh output files newer than their source mtimes. The only caveat is missing `*_sympy.log`/`*_mathematica.log` files in `exec_logs/`, but the canonical refreshed `.txt` outputs (which the user authorizes as the verifier's source of truth per the hard rules) supply the required exit-status evidence (no `FAIL:`, no traceback, completion banners present).

stage 077: verified

