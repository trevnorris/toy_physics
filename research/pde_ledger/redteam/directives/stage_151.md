---
unit_id: 151
batch: IV.6
created_at: 2026-05-28T03:22:59Z
findings_count: 3
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex directive — unit 151

Apply each finding below in order. After applying, append an `## Applied: F<n>` block under that finding with: `files_changed`, `summary` (one sentence), and `deviation` (or "none").

If a finding's required change is ambiguous or unsafe to apply mechanically, append `## Blocked: F<n>` with a question instead — skip that finding, continue with the rest.

Do NOT introduce new features, refactors, or stylistic changes. Edit exactly the file:line ranges named.

Do NOT run python or mathematica. Only edit files.

Do NOT touch paper.tex, notes/, or any prose documents. The red-team only modifies scripts.

## F1 — tautological_check

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:23-60`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:26-56`

**Issue:**
The two existing `expect_zero` / `expectZero` calls in each script check `delta_g + CovcR == 0` and `delta_S + CovKR == 0`. Because `delta_g` is defined as `-(cRbar - cbar*Rbar)` and `CovcR` as `cRbar - cbar*Rbar` (and similarly for the `K` pair), these residuals are zero by literal substitution of the definitions — no physics is exercised. The paper requires the moment-shift identities to be derived from `δg = ∫₀¹ c(x) δΣ_*(x) dx` and `δS = ∫₀¹ K_q(x) δΣ_*(x) dx`, with `δΣ = -Σ_*(R - <R>_*)`. Replace the tautological setup with an integration-based derivation.

**Required change (SymPy, `scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`):**

Replace the body from line 23 to line 60 (everything after the `expect_zero` helper) with the following construction. Keep the `from __future__ import annotations`, `import sympy as sp`, `banner`, and `expect_zero` helper at the top untouched.

Before-state (current lines 23–60) is the existing tautological derivation; after-state must contain (in this order):

1. `banner("FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION")`.
2. Declare symbols:
   - `x = sp.symbols("x", real=True)`
   - `Pi_star = sp.symbols("Pi_star", positive=True)`
   - `r1, r2 = sp.symbols("r1 r2", real=True)`
   - `gprime, AT, BT = sp.symbols("gprime AT BT", real=True, nonzero=True)`
3. Define the canonical mouth source and the kernels:
   - `Sigma_star = Pi_star * sp.exp(-Pi_star*x) / (1 - sp.exp(-Pi_star))`
   - `c_kernel = sp.cos(sp.pi*x/2)`
   - `K_kernel = sp.cosh(sp.pi*(1 - x)/2) / sp.cosh(sp.pi/2)`
   - `R_residual = r1*x + r2*x**2`
4. Define the expectation operator:
   - `def mean(f): return sp.integrate(Sigma_star*f, (x, 0, 1))`
5. Compute the five moments:
   - `Rbar = sp.simplify(mean(R_residual))`
   - `cbar = sp.simplify(mean(c_kernel))`
   - `Kbar = sp.simplify(mean(K_kernel))`
   - `cRbar = sp.simplify(mean(c_kernel*R_residual))`
   - `KRbar = sp.simplify(mean(K_kernel*R_residual))`
   - `CovcR = sp.simplify(cRbar - cbar*Rbar)`
   - `CovKR = sp.simplify(KRbar - Kbar*Rbar)`
6. Build the linearized correction and verify centering:
   - `delta_Sigma = -Sigma_star*(R_residual - Rbar)`
   - `expect_zero("<delta Sigma>_*  (centering)", sp.integrate(delta_Sigma, (x, 0, 1)))`
7. Verify the two moment shifts via integration:
   - `delta_g_int = sp.integrate(c_kernel*delta_Sigma, (x, 0, 1))`
   - `delta_S_int = sp.integrate(K_kernel*delta_Sigma, (x, 0, 1))`
   - `expect_zero("delta_g_int + Cov(c,R)", delta_g_int + CovcR)`
   - `expect_zero("delta_S_int + Cov(K,R)", delta_S_int + CovKR)`
8. Verify the bias retuning and traction retuning:
   - `deltaPi = -delta_g_int/gprime`
   - `deltaT  = AT*delta_g_int + BT*delta_S_int`
   - `expect_zero("deltaPi - Cov(c,R)/gprime", deltaPi - CovcR/gprime)`
   - `expect_zero("deltaT + AT*Cov(c,R) + BT*Cov(K,R)", deltaT + AT*CovcR + BT*CovKR)`
9. Keep the existing closing theorem print and the existing `print` statements for `deltaPi` and `deltaT`.

**Required change (Mathematica, `mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl`):**

Replace the body from line 26 (the `banner["FIRST-ORDER ..."]` call) to line 51 (the closing theorem print) — but NOT the helper definitions at lines 1–24 nor the final `Exit[0]` at line 56 — with an independent integration-based derivation. Important: this must NOT be a line-by-line transliteration of the SymPy. To keep it independent, derive the linearized correction by `Series` expansion instead of direct substitution:

After-state (in order):

1. `banner["FIRST-ORDER SELF-CONSISTENT SOURCE CORRECTION"];`
2. Clear and declare symbols:
   - `Clear[x, piStar, r1, r2, gPrime, aT, bT, epsilon];`
   - `$Assumptions = piStar > 0 && Element[{r1, r2, gPrime, aT, bT}, Reals] && gPrime != 0;`
3. Build the source from the full mouth potential and Series-expand to first order:
   - `Phi[x_] := piStar*x + epsilon*(r1*x + r2*x^2)`
   - `unnorm[x_] := Exp[-Phi[x]]`
   - `Z = Integrate[unnorm[x], {x, 0, 1}, Assumptions -> piStar > 0]`
   - `SigmaFull[x_] := unnorm[x]/Z`
   - `SigmaSeries = Normal[Series[SigmaFull[x], {epsilon, 0, 1}]]`
   - `SigmaStar = Coefficient[SigmaSeries, epsilon, 0]` (must equal `piStar Exp[-piStar x]/(1 - Exp[-piStar])`)
   - `deltaSigma = Coefficient[SigmaSeries, epsilon, 1]` (this is the linearized correction; by construction equals `-SigmaStar*(R - <R>)` if the centering identity holds)
4. Define the kernels and residual:
   - `cKernel[x_] := Cos[Pi*x/2]`
   - `kKernel[x_] := Cosh[Pi*(1 - x)/2]/Cosh[Pi/2]`
   - `RResidual[x_] := r1*x + r2*x^2`
5. Define the expectation operator independently:
   - `mean[f_] := Integrate[SigmaStar*f, {x, 0, 1}, Assumptions -> piStar > 0]`
6. Compute moments:
   - `rBar = FullSimplify[mean[RResidual[x]]]`
   - `cBar = FullSimplify[mean[cKernel[x]]]`
   - `kBar = FullSimplify[mean[kKernel[x]]]`
   - `cRBar = FullSimplify[mean[cKernel[x]*RResidual[x]]]`
   - `kRBar = FullSimplify[mean[kKernel[x]*RResidual[x]]]`
   - `covCR = FullSimplify[cRBar - cBar*rBar]`
   - `covKR = FullSimplify[kRBar - kBar*rBar]`
7. Verify centering and that the Series-expanded correction agrees with `-SigmaStar*(RResidual[x] - rBar)`:
   - `expectZero["<deltaSigma>_*  (centering, from Series)", Integrate[deltaSigma, {x, 0, 1}, Assumptions -> piStar > 0]]`
   - `expectZero["deltaSigma + SigmaStar*(R - <R>)", deltaSigma + SigmaStar*(RResidual[x] - rBar)]`
8. Integrate to get the moment shifts and assert:
   - `deltaGInt = Integrate[cKernel[x]*deltaSigma, {x, 0, 1}, Assumptions -> piStar > 0]`
   - `deltaSInt = Integrate[kKernel[x]*deltaSigma, {x, 0, 1}, Assumptions -> piStar > 0]`
   - `expectZero["deltaGInt + Cov(c,R)", deltaGInt + covCR]`
   - `expectZero["deltaSInt + Cov(K,R)", deltaSInt + covKR]`
9. Verify the bias and traction retunings:
   - `deltaPi = -deltaGInt/gPrime`
   - `deltaT  = aT*deltaGInt + bT*deltaSInt`
   - `expectZero["deltaPi - Cov(c,R)/gPrime", deltaPi - covCR/gPrime]`
   - `expectZero["deltaT + aT*Cov(c,R) + bT*Cov(K,R)", deltaT + aT*covCR + bT*covKR]`
10. Keep the closing theorem print and `Exit[0]`.

The key independence: SymPy expands the linearized correction by hand (`-Sigma_star*(R - Rbar)`); Mathematica derives it from `Series[Exp[-Phi[x]]/Z, {epsilon, 0, 1}]` and checks consistency with the hand form. Both engines must reach the same numeric covariance shifts via different code paths.

**Verification command:**
After Codex applies, the verifier will run `redteam exec-sympy 151` and `redteam exec-mathematica 151` and confirm both exit 0; the new transcripts must show all the listed `expect_zero` / `expectZero` lines printing `= 0` and the engines must agree on the covariance forms.

## F2 — insufficient_verification

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py:42-48`
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:40-46`

**Issue:**
The bias retuning `δΠ = -δg/g_*'` and the traction retuning `δT̂ = A_T δg + B_T δS` (paper notes Section 3 boxes) are only `print`ed by both scripts, never asserted. If they were wrong (sign flip, missing factor), the script would still exit 0.

**Required change:**
The two new assertions

- SymPy: `expect_zero("deltaPi - Cov(c,R)/gprime", deltaPi - CovcR/gprime)` and `expect_zero("deltaT + AT*Cov(c,R) + BT*Cov(K,R)", deltaT + AT*CovcR + BT*CovKR)`
- Mathematica: `expectZero["deltaPi - Cov(c,R)/gPrime", deltaPi - covCR/gPrime]` and `expectZero["deltaT + aT*Cov(c,R) + bT*Cov(K,R)", deltaT + aT*covCR + bT*covKR]`

are already specified as step 8/9 of F1. Apply F2 by ensuring those four assertions are present after the F1 rewrite. No separate edit required if F1 is applied in full.

**Verification command:**
Same as F1 — the four new assertions must appear in the script source and produce `= 0` lines in the saved transcripts.

## F3 — mathematica_transliteration

**Target:**
- `/var/projects/toy_physics/research/pde_ledger/mathematica/moving_throat_pde_stage151_first_order_selected_correction_mathematica_audit.wl:26-46`

**Issue:**
The Mathematica script's algebra mirrors the SymPy script line-by-line (same five expectation-value symbols, same two covariance definitions, same sign convention on `delta_g`/`delta_S`, same assertion list). Codex must remove this mirroring under F1 by giving the Mathematica script a genuinely different derivation path.

**Required change:**
Use the `Series`-based linearization in step 3 of F1's Mathematica rewrite (i.e., derive `deltaSigma` from `Series[Exp[-Phi[x]]/Z, {epsilon, 0, 1}]` rather than from a hand-written `-SigmaStar*(R - rBar)`), and verify the two forms agree via the new check `expectZero["deltaSigma + SigmaStar*(R - <R>)", deltaSigma + SigmaStar*(RResidual[x] - rBar)]`. The SymPy script must keep the hand-form derivation (`delta_Sigma = -Sigma_star*(R_residual - Rbar)`) without a Series step. This gives two independent code paths reaching the same identities. No separate edit beyond F1.

**Verification command:**
Reviewer reads the two scripts side by side. The SymPy script's derivation must define `delta_Sigma` directly from `-Sigma_star*(R - Rbar)`; the Mathematica script must derive `deltaSigma` via `Series` of `Exp[-Phi[x]]/Z` and confirm the agreement. Both transcripts must show all covariance, bias, and traction `= 0` lines.
