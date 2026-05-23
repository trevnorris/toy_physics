---
unit_id: 062
batch: III.3
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-22T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: false
---

# Verification — unit 062

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:67-96`: introduced `sigma, phi` channel symbols, built `S_parent = (1/2) Theta_sigma sigma^2 + Lambda_phi sigma phi + (1/2) KX phi^2`, eliminated `sigma_star = sp.solve(sp.diff(S_parent, sigma), sigma)[0]`, substituted back to get `S_eff_phi`, then extracted `gain_from_action = (KX - 2 * S_eff_phi.coeff(phi, 2))/KX` and asserted `expect_zero("G_micro from parent action vs closed form", gain_from_action - G_micro_closed)` against the hand-written closed-form `rho_star g_phi^2 Osp^2/(m cs_star_sq KX Nss)`. The coherence-factor `expect_zero` was removed and replaced by a single documentation `print` (line 88). For Xi_micro (lines 91-96), `kappa` is now declared as an independent symbol up front (line 66), `Xi_micro = sp.simplify(kappa * gain_from_action)` and `Xi_target = rho_star g_phi^2 Osp^2 L^2/(m cs_star_sq TX Nss)` are formed, then `kappa_solved = sp.solve(Xi_micro - Xi_target, kappa)` is computed and `assert kappa_solved == [KX * L**2 / TX]` enforces uniqueness of the dimensional definition.
- `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:57-93`: same parent-action construction via `Solve[D[sParent, sigma] == 0, sigma]` + `Coefficient[sEff, phi, 2]`, augmented with a second symbolic route via `SeriesCoefficient[Series[sEff, {phi, 0, 2}], 2]` (see F3). Coherence check downgraded to a `Print` (line 82). Xi: `kappaRules = Solve[xiMicro == xiTarget, kappa]` + `If[FullSimplify[kappaSolved == kX*ell^2/tX] === True, pass[...], fail[...]]`.

**Assessment:**
The three replacement assertions are non-tautological. `sigma_star` comes from `sp.solve` on the stationarity equation `diff(S_parent, sigma) = 0`; substituting back into `S_parent` and extracting the phi^2 coefficient is an actual algebraic computation, not a re-typed twin. A wrong sign on `Lambda_phi`, a wrong factor in `Theta_sigma`, or a wrong exponent on `Osp` would propagate into `S_eff_phi.coeff(phi, 2)` and break the residual against `G_micro_closed`. The kappa block likewise is substantive: `sp.solve(Xi_micro - Xi_target, kappa)` returns a unique solution only when `Xi_target` matches the structure of `kappa * gain_from_action`; if either were wrong the solver would return a different value (or empty list). The arithmetic check is `assert kappa_solved == [KX * L**2 / TX]`. Coherence drop is correct — `C2 := Osp^2/(Nss Npp)` was a definition, not a theorem.

Output confirms: `G_micro from parent action vs closed form = 0` (SymPy line 18 of output, Mathematica passes via `PASS:`); `sigma_star = -O_sp*g_phi*phi*rho_star/(N_ss*cs_star_sq*m)` and `S_eff_phi = K_X*phi^2/2 - O_sp^2*g_phi^2*phi^2*rho_star/(2*N_ss*cs_star_sq*m)` — the algebraic chain is visibly performed, not retyped. `kappa solved from Xi_micro = Xi_target: K_X*L**2/T_X` matches expectation.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
- `scripts/moving_throat_pde_stage062_parent_action_gain_sympy_audit.py:33-61`: introduced symbol `n_poly = sp.symbols("n_poly", positive=True, integer=True)`, rewrote `U_general = K * rho**n_poly / (n_poly - 1)` (preserving the stage's `U = K rho^5/4` normalization at n=5), defined `h_general`, `hprime_general`, `cs_sq_general = sp.diff(K * rho**n_poly, rho)/m` symbolically in `n_poly`, asserted `expect_zero("h'(rho) = m c_s^2/rho (general polytrope)", hprime_general - m*cs_sq_general/rho)` BEFORE specializing. Then `subs_n5 = {n_poly: 5}` is applied and the n=5 specialization is also asserted. Finally an inconsistency probe `cs_sq_wrong = sp.diff(K * rho**(n_poly + 1), rho)/m` followed by `assert sp.simplify(residual_wrong) != 0` confirms the assertion machinery would catch an exponent mismatch.
- `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:28-55`: analogous structure with `nPoly`, `uGeneral`, `csSqWrong`, `residualWrongN5`, and an `If[... === 0, fail[...], Print[...]]` inconsistency probe.

**Assessment:**
The general-`n_poly` form is now substantive: had the polytropic identity been stated incorrectly (e.g., `c_s^2` defined with `rho^(n+1)` instead of `rho^n`), the residual would not simplify to zero for general `n_poly`. The wrong-exponent probe is logged: SymPy output line 12 shows `Inconsistency probe (n+1 in c_s^2 only): K*rho**3*(5 - 6*rho)` — explicitly nonzero, confirming the test machinery can detect the wrong index. The `(n - 1)` normalization preserves the stage's original n=5 form (`U = K rho^5/4`) as verified by output line 7 `U(rho) = K*rho**5/4`. Both general and specialized residuals print `= 0`.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
- `mathematica/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.wl:69-79`: in addition to `Coefficient[sEff, phi, 2]` (matching SymPy's `.coeff`), Mathematica also computes `gainFromSeries = (kX - 2*SeriesCoefficient[Series[sEff, {phi, 0, 2}], 2])/kX` via the series-coefficient route. The two routes are cross-checked: `expectZero["Mathematica two-route consistency", gainFromAction - gainFromSeries]` BEFORE comparison to the closed-form. SymPy does not use `SeriesCoefficient` or `Series`. This is Option A from the directive.

**Assessment:**
The .wl file now contains a primitive (`SeriesCoefficient`/`Series`) that does not appear in the .py file. The Mathematica gain is derived twice via independent symbolic chains and cross-verified before being compared to the closed form, where SymPy only does it once. This satisfies the second-engine policy — an algebraic error in one route (e.g., a miscoded `Coefficient`) would now be flagged by disagreement with the `Series` route. Output line 20-21 of the Mathematica `.txt`: `Mathematica two-route consistency = 0 / PASS:`. Confirmed.

## Exec log assessment

**SymPy:** exit=0 (inferred — no `.log` file was captured at `redteam/exec_logs/stage_062_sympy.log`, but `scripts/output/moving_throat_pde_stage062_parent_action_gain_sympy_audit.txt` exists with mtime `2026-05-22 19:40:39` (script mtime `19:38:25`) and ends with `All Stage 45 symbolic checks passed.`). Notable lines:
- `h'(rho) = m c_s^2 / rho (general polytrope) = 0` — general-`n_poly` identity holds.
- `Inconsistency probe (n+1 in c_s^2 only): K*rho**3*(5 - 6*rho)` — wrong-exponent probe nonzero.
- `G_micro from parent action vs closed form = 0` — non-tautological gain comparison passes.
- `kappa solved from Xi_micro = Xi_target: K_X*L**2/T_X` — unique kappa solution matches definition.

**Mathematica:** exit=0 (same inference — `mathematica/output/moving_throat_pde_stage062_parent_action_gain_mathematica_audit.txt` mtime `2026-05-22 19:40:44`, script mtime `19:39:26`, ends with `Stage 062 Mathematica audit passed.` and the `.wl` ends with `Exit[0]`). Notable lines:
- `PASS: h'(rho) = m c_s^2 / rho (general polytrope)`
- `inconsistency probe nonzero: capitalK*(5 - 6*rho)*rho^3`
- `PASS: Mathematica two-route consistency`
- `PASS: gMicro from parent action vs closed form`
- `PASS: kappa solution equals kX ell^2/tX`

**Output freshness:** confirmed — SymPy output mtime `2026-05-22 19:40:39` > script mtime `2026-05-22 19:38:25`; Mathematica output mtime `2026-05-22 19:40:44` > script mtime `2026-05-22 19:39:26`. Both regenerated post-fix.

Caveat: the orchestrator captured only `stage_062_diff.patch` in `redteam/exec_logs/`, not the `_sympy.log` or `_mathematica.log` files. The saved `.txt` outputs in the script output directories serve as the de-facto exec evidence; both are fresh and consistent.

## Material-change assessment

`material_change`: false.

The closed-form result `G_micro = rho_star g_phi^2 Osp^2/(m cs_star_sq KX Nss)` and `Xi_micro = kappa G_micro` with `kappa = KX L^2/TX` are unchanged. Codex restructured HOW the assertions are derived (parent-action elimination, polytropic-index parameterization, two-route Mathematica check) without altering the published formulas. Downstream units that depend on the closed-form `G_micro`/`Xi_micro` numerics see no change.

## Side observations (non-blocking)

- The SymPy `Xi_micro` block at line 91-96 retains an `assert ... != 0` style inside a `print` chain that already does heavy lifting; minor stylistic only.
- Banner string says "STAGE 45" in both engines while the unit is 062; cosmetic only, predates this directive.
- The Mathematica `csSqGeneral` uses `D[capitalK*rho^nPoly, rho]` rather than `D[capitalK*rho^nPoly/(nPoly-1)*(nPoly), rho]` style normalization, mirroring the SymPy choice consistently. This is the chosen convention `c_s^2 = (1/m) d(K rho^n)/drho`, preserved as the directive instructed.

## Verdict justification

All three findings (F1 tautological_check, F2 insufficient_verification, F3 mathematica_transliteration) are resolved with substantive, non-tautological replacements. The parent-action Gaussian elimination genuinely derives the gain rather than retyping its closed form. The polytropic-index parameterization gives the EOS identity a real failure mode (verified by the wrong-exponent probe being nonzero). The Mathematica file now uses `SeriesCoefficient`/`Series` — primitives absent from the .py file — and cross-checks two independent symbolic routes before comparing to the closed form. Both scripts produce fresh outputs ending in their respective "passed" lines, and no regressions are visible in the diff. Verdict: verified.
