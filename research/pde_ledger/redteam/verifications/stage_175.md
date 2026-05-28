---
unit_id: 175
batch: V.1
verifier_model: claude-opus-4-8[1m]
verify_date: 2026-05-28T16:30:00-06:00
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 175

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- sympy `scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:127-142`
- mathematica `mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:91-105`

The old `X - X` tautology line — `expect_zero("Sigma_N - (2 dln Lambda - dK)", Sigma_N_direct - (2 * dlog(expr_L) - kappa))` / `expectZero["Sigma_N - (2 dln Lambda - dK)", sigmaNDirect - (2*dlog[exprL] - kappa)]` — has been removed from both engines. The cached `expr_L`/`exprL` and `expr_ratio`/`exprRatio` variables are gone entirely (grep confirms no surviving references outside explanatory comments). `Sigma_N_direct` is now built from a freshly-named `expr_PoverDelta_phys = (P/Delta).subs(subs_hat).subs(subs_eps)` (sympy:132 / wl:95), and the surviving assertion is `Sigma_N - dln(Lambda^2/K)` (sympy:135 / wl:98), comparing `2*dlog(P/Delta) - kappa` against `dlog(Lambda^2/K)`.

This matches the orchestrator note: on run, the applier's cross-check variant (`2 dln(P/Delta) - 2 dln Lambda` plus the relabeled line) reduced to a simplify-commutes identity, so the orchestrator switched to the directive's explicitly-blessed MINIMAL resolution — keep only the `Sigma_N - dln(Lambda^2/K)` check and drop the two near-trivial lines in both engines.

**Assessment:**
Correct and complete. The X-X tautology is genuinely gone — there is no `expect_zero`/`expectZero` of the form `Sigma_N_direct - (2*dlog(<same arg>) - kappa)` anywhere; the only remaining mentions of those phrases are inside comments documenting the removal. The surviving `Sigma_N - dln(Lambda^2/K)` check is non-tautological: the LHS subtracts `kappa` explicitly while the RHS derives the matching `-kappa` from differentiating `log(Lambda^2/K) = 2 log(Lambda) - log(K)` (with `dlog(K) = kappa`); a wrong `kappa` coefficient or sign would break the identity, so the `kappa = delta_K` content is load-bearing. The substantive scale-invariance content additionally sits in the genuine homogeneity check `N0 - Lambda^2` (line 83/61) and is propagated into the common-shape branch `Sigma_N + dK` and the new F2 aggregate. This is the intended fix. Resolved.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- sympy `scripts/moving_throat_pde_stage175_wall_normalized_load_shape_sympy_audit.py:157-165`
- mathematica `mathematica/moving_throat_pde_stage175_wall_normalized_load_shape_mathematica_audit.wl:114-118`

A new weighted-aggregate no-go assertion `Xi_load (all shapes frozen) + dK` was added immediately after the common-shape-branch check in both engines, exactly as the directive specified: introduce non-negative weights `rho1, rho2`, form `Xi_load_frozen = (rho1 + rho2) * Sigma_N_common`, substitute `rho1 + rho2 -> 1`, and assert the residual `+ kappa == 0`. The sympy version uses `.subs(rho1 + rho2, 1)`; the `.wl` uses the primary `/. (rho1 + rho2) -> 1` replacement (not the `rho2 -> 1 - rho1` fallback).

**Assessment:**
Correct, matches the directive verbatim, and non-tautological in a meaningful sense. The check depends on two facts holding simultaneously: (a) `Sigma_N_common` actually equals `-kappa` (the per-port no-go from line 156/113), and (b) the `sum rho = 1` substitution firing on the `(rho1+rho2)*(-kappa)` form to yield `-kappa`. Both exec logs show the line landing `= 0` / `PASS`, confirming the replacement fired in both engines. This anchors the `sum_r rho_r^(N) = 1` aggregation step that the headline no-go theorem (`Xi_load = -delta_K`) relies on, which previously existed only as printed prose. Resolved.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- banner (step 1): `mathematica/...mathematica_audit.wl:31` now reads `banner["STAGE 175 — WALL-NORMALIZED LOAD/SHAPE FACTORIZATION"]`; the SymPy banner (`...sympy_audit.py:39`) likewise corrected. grep for `STAGE 158`/`Stage 158` returns no matches in either file.
- F1 fix (step 2): applied independently in the `.wl` differential block (see F1 above) so the `Sigma_N` check is no longer a transliterated tautology.
- step 3 (dlogSeries series-coefficient route): intentionally NOT applied — accepted as a policy mirror.

**Assessment:**
Resolved, policy-accepted for step 3. Steps 1 and 2 are correct and verified: the wrong "STAGE 158" banner is fixed in both engines, and the F1 tautology is removed from the `.wl` as well as the `.py`. Step 3 (replacing `dlog` with `dlogSeries[expr_] := Coefficient[Normal[Series[Log[expr], {eps, 0, 1}]], eps]`) was left as a policy mirror rather than forced — the `.wl` Sigma_N block keeps `dlog`, consistent with the project's MATHEMATICA_MIRROR_POLICY default and the unit-169 F2 disposition this batch. The directive itself made step 3 conditional ("If a full re-derivation is out of scope ... at minimum apply (c) the banner correction and the F1 fix independently in each engine"), so applying steps 1+2 and policy-accepting step 3 satisfies the directive's stated minimum. This is not a `partial` outcome: the directive's required mechanical change (banner + independent F1 fix in each engine) is fully in place, and the residual transliteration of the structurally-identical algebra is a known, accepted mirror policy rather than an unaddressed finding. Resolved.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `N0 - Lambda^2 = 0` (genuine homogeneity check, anchors Sigma_N content)
- `Sigma_N - dln(Lambda^2/K) = 0` (the surviving, non-tautological differential check post-F1)
- `Common-shape branch Sigma_N + dK = 0`
- `Xi_load (all shapes frozen) + dK = 0` (the new F2 weighted-aggregate check)
No traceback / no `AssertionError`; `expect_zero` raises on any non-zero residual, so a clean transcript with every line `= 0` implies exit 0.

**Mathematica:** exit=0. Notable lines:
- `PASS: N0 - Lambda^2`
- `Sigma_N - dln(Lambda^2/K) = 0` / `PASS: Sigma_N - dln(Lambda^2/K)`
- `PASS: Common-shape branch Sigma_N + dK`
- `Xi_load (all shapes frozen) + dK = 0` / `PASS: Xi_load (all shapes frozen) + dK`
- closing `Stage 175 Mathematica audit passed.` (script ends with `Exit[0]`; `expectZero` failure path calls `fail` -> `Exit[1]`, and no `FAIL` appears)
Both banners now print `STAGE 175 — WALL-NORMALIZED LOAD/SHAPE FACTORIZATION`.

**Output freshness:** confirmed. sympy script mtime 16:16:28 < sympy output mtime 16:16:57; mathematica script mtime 16:16:36 < mathematica output mtime 16:17:13. Both saved `.txt` outputs were re-generated after the post-fix script edits. (The MANIFEST mtimes are stale start-of-session snapshot values and are not used for the exit determination.)

## Material-change assessment

`material_change`: false.

No derived result that downstream units depend on was altered. F1 removed two near-trivial / tautological assertions and renamed an intermediate variable (`expr_ratio`/`exprL` -> `expr_PoverDelta_phys`) without changing any computed value — the surviving `Sigma_N - dln(Lambda^2/K)` identity already passed and still passes. F2 added a new internal self-consistency assertion (`Xi_load (all shapes frozen) + dK`) that only confirms the already-established per-port `Sigma_N = -kappa` aggregates correctly under `sum rho = 1`; it introduces no new constant or result. F3 is a cosmetic banner correction plus the F1 fix mirror. All edits are local to the unit-175 audit scripts and produce no new numeric/symbolic output consumed elsewhere. Downstream units > 175 do not need re-audit on account of unit 175.

## Side observations (non-blocking)

- The printed "Conclusions" block (sympy:172 / wl:126) still echoes `Sigma_N = d ln(Lambda^2/K) = 2 d ln Lambda - dK` as prose. This is fine — it is descriptive output, not an assertion, and the `= 2 d ln Lambda - dK` equality is the log-algebra identity, not a claim that an independent check exists. No action needed.
- The `dlog` helper applies `simplify`/`FullSimplify` to its argument before differentiating, which is why the engines remain structurally similar even after F1. This is the accepted mirror-policy posture for this batch (F3 step 3), not a defect.

## Verdict justification

All three findings are resolved. F1: the `X - X` tautology is genuinely removed from both engines (no surviving `expr_L`/`exprL` references; the remaining `Sigma_N - dln(Lambda^2/K)` check is load-bearing on `kappa = delta_K`), matching the orchestrator's blessed minimal resolution. F2: the weighted-aggregate `Xi_load (all shapes frozen) + dK` check exercising `sum rho = 1` is present in both engines and is non-tautological, landing `= 0` in both logs. F3: banner corrected to STAGE 175 in both engines and the F1 fix mirrored in the `.wl`; step 3 is a sanctioned policy mirror (not a gap). The git diff shows no collateral edits, both exec logs pass with all assertions `= 0`/`PASS`, and both outputs are fresh relative to the post-fix scripts. Verdict: verified.
