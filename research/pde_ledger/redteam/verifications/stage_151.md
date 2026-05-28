---
unit_id: 151
batch: IV.6
verifier_model: claude-opus-4-7[1m]
verify_date: 2026-05-28T00:00:00Z
verdict: verified
sympy_exit: 0
mathematica_exit: 0
findings_resolved: 3
findings_total: 3
material_change: true
---

# Verification — unit 151

## Per-finding outcomes

### F1 — tautological_check

**Classification:** resolved

**What changed:**
SymPy script (`scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:23-130`) rewritten to use `mpmath` numerical integration at 40 dps with concrete `Pi_star=1.50882951349316`, `r1=1.7`, `r2=-0.9`, and concrete `gprime`, `AT`, `BT`. Defines `Sigma_star`, `c_kernel`, `K_kernel`, `R_residual`, `mean(f)`, computes `Rbar, cbar, Kbar, cRbar, KRbar, CovcR, CovKR` by quadrature, builds `delta_Sigma = -Sigma_star*(R - Rbar)` (hand form), then asserts `delta_g_int = -CovcR` and `delta_S_int = -CovKR` via independent integration. Mathematica script (`mathematica/...wl:26-74`) rewritten to derive `deltaSigma` independently via `Series[Exp[-Phi]/Z, {epsilon, 0, 1}]`, then verifies it agrees with the hand form, then derives covariance shifts via `Integrate`.

**Assessment:**
Non-tautological. SymPy diffs in log are ~1e-42 (well below 1e-30 tol). Mathematica returns symbolic 0 after `FullSimplify`. Deviation from directive: SymPy uses `mpmath` (numeric) instead of `sympy.integrate` (symbolic) and helper renamed `expect_zero`→`expect_close`. Orchestrator note acknowledges this; functionally equivalent — the integrals would not numerically equal the covariance combinations unless the underlying identity holds.

### F2 — insufficient_verification

**Classification:** resolved

**What changed:**
Both scripts now assert `deltaPi - Cov(c,R)/gprime = 0` and `deltaT + AT*Cov(c,R) + BT*Cov(K,R) = 0` (SymPy lines 116-122; Mathematica lines 63-64). Both produce `= 0` lines in transcripts.

**Assessment:**
Bias and traction retunings now exercised by assertion as required.

### F3 — mathematica_transliteration

**Classification:** resolved

**What changed:**
Mathematica derives `deltaSigma` via `Series` expansion of `Exp[-Phi[x]]/Z` to first order in `epsilon` (wl:31-37), then verifies agreement with the SymPy-style hand form via `expectZero["deltaSigma + SigmaStar*(R - <R>)", ...]` (wl:54). SymPy keeps the hand-form path. Two genuinely different code paths reaching the same identities.

**Assessment:**
Independent derivation confirmed. The cross-check `deltaSigma + SigmaStar*(R - rBar) = 0` is a non-trivial bridge between the two paths and it PASSes in the Mathematica log.

## Exec log assessment

**SymPy:** exit=0. Notable lines:
- `<1>_* = 1 (canonical normalization) = 1.0   (target 1.0, diff 0.0)`
- `delta_g_int = -Cov(c,R) = 0.05967492458671911...   (target 0.05967..., diff 2.15e-42)`
- `<delta_Sigma>_*  (centering) = 1.38e-42   (target 0.0, diff 1.38e-42)`
- `deltaT = -AT*Cov(c,R) - BT*Cov(K,R) = -0.2487...   (target -0.2487..., diff 5.74e-42)`

All six numeric checks pass at tol=1e-30 with residuals ~1e-42.

**Mathematica:** exit=0. Notable lines:
- `<deltaSigma>_*  (centering, from Series) = 0   PASS`
- `deltaSigma + SigmaStar*(R - <R>) = 0   PASS`
- `deltaGInt + Cov(c,R) = 0   PASS`
- `deltaT + aT*Cov(c,R) + bT*Cov(K,R) = 0   PASS`

All six symbolic `expectZero` checks pass.

**Output freshness:** Confirmed. SymPy: script mtime 1779989270 < txt 1779989407. Mathematica: script 1779984192 < txt 1779989491.

## Material-change assessment

`material_change`: true.

The script now produces concrete numerical values for `deltaPi = -0.8352526754408...` and `deltaT = -0.2487255954387...` keyed to the chosen `Pi_star=1.50882951349316` (canonical Family-1 mouth bias, per docstring referencing stage 156) and the example residual coefficients `r1=1.7, r2=-0.9`. These concrete numbers were not previously asserted. Downstream stages that consume `deltaPi`/`deltaT` (any stage > 151, particularly stage 152 traction retuning) should be flagged stale.

## Side observations (non-blocking)

- SymPy script uses `mpmath` numeric quadrature rather than the symbolic SymPy integration prescribed in the directive. Orchestrator note explicitly approves this deviation. The choice is reasonable because the symbolic `sp.integrate(c*Sigma_star, ...)` would yield closed-form expressions involving `Si`/`Ci`, and high-precision numerical agreement is at least as strong as symbolic equivalence for non-tautology.
- The `expect_close` default tol parameter is `1e-15`, but all call sites override with `tol=1e-30`; observed residuals are ~1e-42, comfortably within either bound.
- Pi_star value `1.50882951349316` is hardcoded and referenced as "canonical Family-1 mouth bias (notes / stage 156)" — this introduces a forward dependency on stage 156. If stage 156's canonical Pi_star ever changes, this script's numerics will need a refresh.

## Verdict justification

All three findings resolved. Both scripts now exercise the paper's deliverables non-tautologically: SymPy via high-precision mpmath quadrature, Mathematica via independent `Series` expansion of the full mouth profile. The centering identity, both moment-shift identities, the bias retuning, and the traction retuning are all asserted (not just printed). The two engines reach the same boxed identities through genuinely different code paths.
