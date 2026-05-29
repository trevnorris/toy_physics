---
unit_id: 151
batch: IV.6
created_at: 2026-05-28T19:25:17-06:00
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-28T21:07:31-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 151

Apply each non-paper_misalignment finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

After editing, RUN the affected scripts (`python3 <path>` for SymPy, `math -script <path>` for Mathematica) and iterate until they exit 0 with all in-file checks passing. Getting the scripts to run cleanly is your job; the orchestrator independently re-runs afterward.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — insufficient_verification

**Target:** `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py` (whole file)

**Issue:** The SymPy script verifies exact symbolic identities (the paper's deliverables #2 and #3, which must hold for all `Π_* > 0` and all residual coefficients) only as a single-point floating-point spot-check: it imports `mpmath` (line 24), fixes numeric anchors `Pi_star = mp.mpf("1.50882951349316")`, `r1 = 1.7`, `r2 = -0.9`, `gprime`, `AT`, `BT` (lines 50-55), integrates with `mp.quad`, and accepts via `expect_close` with tolerance (nonzero diffs ~1e-42 in the saved transcript). A single numeric anchor cannot distinguish the true identity from one that merely coincides at that point. Additionally, the SymPy script hand-codes the linearized correction `delta_Sigma(x) = -Sigma_star(x)*(R_residual(x) - Rbar)` (lines 99-100) instead of deriving it from the exponential expansion `e^{-Φ_*} = e^{-Π_* x - R_*(x)}` (paper deliverable #1), which is therefore not exercised on the SymPy side at all. The Mathematica engine already does the full symbolic derivation; SymPy must match it.

**Required change:**
Rewrite the SymPy script body so the verification is exact-symbolic, not numeric. Concretely:

1. Replace the `import mpmath as mp` line (line 24) with `import sympy as sp`.

2. Remove the `expect_close` helper (lines 34-38). Add an exact-zero helper, e.g.:
   ```python
   def expect_zero(name, expr):
       res = sp.simplify(expr)
       print(f"{name} = {res}")
       assert res == 0, f"{name} nonzero: {res}"
   ```
   (If `sp.simplify` does not collapse a residual to a bare `0`, fall back to `assert sp.simplify(expr).equals(0)`; do not use any numeric tolerance.)

3. Replace the numeric anchors (lines 50-55) and `mp.pi/2` with symbolic declarations:
   ```python
   x = sp.symbols("x", real=True)
   eps = sp.symbols("epsilon", real=True)
   Pi_star = sp.symbols("Pi_star", positive=True)
   r1, r2 = sp.symbols("r1 r2", real=True)
   gprime = sp.symbols("gprime", real=True, nonzero=True)
   AT, BT = sp.symbols("A_T B_T", real=True)
   k = sp.pi / 2
   ```

4. Derive `Sigma_star` and `delta_Sigma` from first principles (mirroring the Mathematica `Series` derivation, not the hand form):
   ```python
   Phi = Pi_star * x + eps * (r1 * x + r2 * x**2)
   Z = sp.integrate(sp.exp(-Phi), (x, 0, 1))
   Sigma_full = sp.exp(-Phi) / Z
   series = sp.series(Sigma_full, eps, 0, 2).removeO()
   Sigma_star = series.coeff(eps, 0)
   delta_Sigma = series.coeff(eps, 1)
   ```
   Define the kernels and residual as symbolic expressions in `x`:
   ```python
   c_kernel = sp.cos(k * x)
   K_kernel = sp.cosh(k * (1 - x)) / sp.cosh(k)
   R_residual = r1 * x + r2 * x**2
   ```
   and a `mean` helper: `mean = lambda f: sp.integrate(Sigma_star * f, (x, 0, 1))`.

5. Replace the moment/covariance block (lines 80-122) with symbolic computations and exact-zero assertions, covering all of M1-M7 below. Keep informative `print` lines for the simplified moments/covariances if desired, but the load-bearing checks must be `expect_zero`.

6. Keep the closing `print("Theorem:")` narrative block (lines 127-129) as-is.

After the rewrite there must be no `mpmath`, no `mp.mpf`, no `mp.quad`, no `expect_close`, and no fixed numeric values for `Pi_star`, `r1`, `r2`, `gprime`, `AT`, `BT` anywhere in the file. Do not change the filename, the docstring topic, or the banner text — but update the docstring's "numerically integrating … to high precision" sentence (lines 11-14) to describe the symbolic derivation instead (e.g., "symbolically integrating over symbolic Pi_star, r1, r2 and asserting the exact covariance identities reduce to zero").

**Claim manifest** (the new script must independently verify each, exactly):

- **M1** (deliverable #1, currently missing): the series-derived first-order correction equals the centered hand form —
  `delta_Sigma + Sigma_star*(R_residual - Rbar) == 0`, where `Rbar = mean(R_residual)`.
- **M2** (normalization): `mean(1) - 1 == 0`, i.e. `sp.integrate(Sigma_star, (x,0,1)) - 1 == 0`.
- **M3** (centering / mass preservation): `sp.integrate(delta_Sigma, (x,0,1)) == 0`.
- **M4** (deliverable #2, δg): with `CovcR = mean(c_kernel*R_residual) - mean(c_kernel)*Rbar`,
  `sp.integrate(c_kernel*delta_Sigma, (x,0,1)) + CovcR == 0`.
- **M5** (deliverable #2, δS): with `CovKR = mean(K_kernel*R_residual) - mean(K_kernel)*Rbar`,
  `sp.integrate(K_kernel*delta_Sigma, (x,0,1)) + CovKR == 0`.
- **M6** (deliverable #3, δΠ): with `delta_g_int = sp.integrate(c_kernel*delta_Sigma, (x,0,1))` and `deltaPi = -delta_g_int/gprime`,
  `deltaPi - CovcR/gprime == 0`.
- **M7** (deliverable #3, δT_m): with `delta_S_int = sp.integrate(K_kernel*delta_Sigma, (x,0,1))` and `deltaT = AT*delta_g_int + BT*delta_S_int`,
  `deltaT + AT*CovcR + BT*CovKR == 0`.

Each Mi must be wrapped in an `expect_zero("<label>", <expr>)` call. No tolerances. The names should make the transcript self-documenting (e.g. `expect_zero("delta_Sigma + Sigma_star*(R - <R>)  [M1]", ...)`).

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 151` and confirm: (a) the refreshed `scripts/output/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.txt` shows exact `= 0` results (no `diff ...e-NN` tolerance lines) for M1-M7; (b) a check corresponding to M1 (`delta_Sigma + Sigma_star*(R - <R>)`) now appears, which had no SymPy counterpart before; (c) the script exits 0; (d) no `mpmath` import remains in the `.py`.

## Applied: F1

- files_changed:
  - `scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`
- summary: Rewrote the SymPy audit to derive the first-order source correction symbolically and assert exact zero residuals for M1-M7.
- deviation: none
