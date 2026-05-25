---
unit_id: 075
batch: III.4
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-25T00:00:00Z
verdict: verified
sympy_exit: n/a
mathematica_exit: n/a
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 075

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- SymPy: `scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:96-104` — the lone tautological `assert Upsilon_expr == 100 * Theta` is gone; replaced with two non-trivial asserts that build `Upsilon_*_from_Theta = sp.simplify(alpha_r**2 * Theta_<branch>)` using the actually-constructed `Theta_fail` / `Theta_suff` (which were themselves built from `Upsilon_<branch> / alpha_r**2` via the `Delta_0`/`Delta_inf` closed forms — lines 72-79). Both diff prints land on `0` in the output transcript (lines 31-32).
- Mathematica: `mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:95-96` — single tautological `expectZero["Upsilon_w(reference) - 100 Theta_w", alphaR^2*thetaW - 100*thetaW]` replaced by two `expectZero` calls against the constructed `upsilonFail/thetaFail` and `upsilonSuff/thetaSuff`. Output lines 23-26 show both report `0` and `PASS`.

**Assessment:**
The replacement matches the directive verbatim. The new SymPy asserts are non-tautological because `Upsilon_fail` and `Theta_fail` are independently constructed from `Pe_req / (Lambda_ell**2 * Deltainf)` and `Upsilon_fail / alpha_r**2` respectively — the residual `Upsilon_fail - alpha_r**2 * Theta_fail` would not collapse to `0` if either of those two construction steps had an algebraic bug (e.g. a wrong factor on the `alpha_r**2` rescaling, or a sign error in the `Delta`-substitution). Same logic in the Mathematica side. No collateral edits.

### F2 — hardcoded_result

**Classification:** resolved

**What changed:**
- SymPy: `scripts/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.py:40-61` — inserted free-symbol identity block exactly as the directive specified. New `alpha_sym, eta_sym = sp.symbols(...positive=True, real=True)` carry the closed forms as functions of *symbols*, and `sp.simplify` is asked to verify two algebraic identities; both `assert ... == 0` lines guard them. Output lines 5-6 show both identities print `0`.
- Mathematica: `mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:83-93` — `Module[{aSym, eSym, ...}, ...]` block with `ClearAll[aSym, eSym]` and `Assuming[aSym > 0 && eSym > 0, FullSimplify[...]]` for both `Delta_0` and `Delta_inf` identities. Output lines 19-22 show both report `0` and `PASS`. The nine `expectApprox` calls remain as informational per directive (lines 98-105).

**Assessment:**
The new identity checks use *symbolic* `alpha_sym`/`eta_sym` (SymPy) and `aSym`/`eSym` (Mathematica), not the numeric substitutions `alpha = sqrt(12321/5)` and `eta = 37`. Mathematica's `FullSimplify` independently re-derives the algebraic relation rather than evaluating the same numeric expression as SymPy. A wrong factor of 2 or sign flip in the stated closed forms would now fail at the identity step. The previously-flagged "calculator self-consistency" hole is closed. No deviation from directive.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- Mathematica: `mathematica/moving_throat_pde_stage075_family1_threshold_window_mathematica_audit.wl:78-82` — the required explanatory comment is inserted immediately above the `Module[{aSym, eSym, ...}, ...]` block from F2, noting that this identity check constitutes the independent-derivation leg required by the second-engine policy.

**Assessment:**
The comment is present verbatim from the directive. F3's structural correction was delegated to F2's free-symbol identity check, which is now in place (see F2 above). No additional edits required or made.

## Exec log assessment

**SymPy:** exit=n/a. The orchestrator-captured `redteam/exec_logs/stage_075_sympy.log` is absent; I did not re-run python per the prompt's rule. Instead I read the refreshed transcript at `scripts/output/moving_throat_pde_stage075_family1_threshold_window_sympy_audit.txt` (mtime May 22 23:16, fresher than script mtime May 22 23:13). Notable lines:
- Line 5: `Delta_0 algebraic identity (free alpha, eta) = 0` — F2 identity check passes.
- Line 6: `Delta_inf algebraic identity (free alpha, eta) = 0` — F2 identity check passes.
- Lines 31-32: `Upsilon_fail - alpha_r^2 * Theta_fail = 0` / `Upsilon_suff - alpha_r^2 * Theta_suff = 0` — F1 round-trip residuals print 0 (which corresponds to the two `assert ... == 0` lines not raising).

**Mathematica:** exit=n/a. The orchestrator-captured `redteam/exec_logs/stage_075_mathematica.log` is absent; I read the refreshed `.txt` instead (mtime May 22 23:16, fresher than script mtime May 22 23:13). Notable lines:
- Lines 19-22: `Delta_0 algebraic identity (free alpha, eta) = 0` / `PASS`, then `Delta_inf algebraic identity (free alpha, eta) = 0` / `PASS` — F2/F3 identity checks pass.
- Lines 23-26: `Upsilon_fail - alphaR^2 * Theta_fail = 0` / `PASS`, then `Upsilon_suff - alphaR^2 * Theta_suff = 0` / `PASS` — F1 round-trip identities pass.
- Lines 27-42: all nine `expectApprox` numeric checks report `diff = 0` and `PASS` (still informational, as the directive explicitly left them in place).
- Line 44: `Stage 075 Mathematica audit passed.`

**Output freshness:** confirmed. Both `.txt` outputs have mtime `May 22 23:16` (later than the corresponding scripts' `May 22 23:13`). Outputs are fresh and reflect the post-fix scripts.

## Material-change assessment

`material_change`: false.

No printed symbolic content or numeric value in either transcript changed. `Delta_0`, `Delta_inf`, `Upsilon_*`, `Xi_*`, `Theta_*` are identical to the pre-fix values (compare the audit report's table of pre-fix targets to the new output: `Delta_0 = 0.00017330207902152514906`, `Delta_inf = 0.020144756554052159427`, `Upsilon_fail/Pe_req = 0.036260561797293886969`, etc. — all unchanged). The edits only added new assertions and a comment; they did not alter any derivation route or printed quantity. No downstream stages need to be marked stale on numeric grounds.

## Side observations (non-blocking)

- The SymPy script's banner still reads `STAGE 58` and the Mathematica banner reads `STAGE 058`, while the file is unit 075 — this is cosmetic, predates the edits in this batch, and is outside the F1-F3 scope.
- The SymPy script defines `Theta = sp.symbols("Theta_w", positive=True, real=True)` on line 77 but never uses the bare `Theta` symbol after the directive's edit (the previous tautological assertion was the sole consumer). This is now an unused symbol; harmless and not in scope to remove.

## Verdict justification

All three findings are resolved with directive-faithful edits. The two new SymPy asserts and four new Mathematica `expectZero` calls now exercise (i) the closed-form algebraic identities for `Delta_0` and `Delta_inf` with free symbols `alpha`/`eta` (closing F2 and F3 simultaneously), and (ii) the `Upsilon_<branch> = alpha_r^2 * Theta_<branch>` round-trip on the actually-constructed thresholds rather than the trivial `100*Theta_w == 100*Theta_w` identity (closing F1). The refreshed output transcripts show every new check landing on `0` / `PASS`, and no previously-printed symbolic content or numeric value changed.

stage 075: verified
