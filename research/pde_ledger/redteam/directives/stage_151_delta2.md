---
unit_id: 151
batch: IV.6
iteration: 3
created_at: 2026-05-28T22:30:00-06:00
parent_directive: stage_151_delta1.md
supersedes: stage_151_delta1.md
findings_count: 1
stop_cold: null
applied: true
applied_at: 2026-05-28T23:09:55-06:00
findings_applied: 1
findings_blocked: 0
verification_status: pending
needs_user_resolution: false
---

# Codex delta-directive — unit 151 (iteration 3, supersedes delta1)

**Why this supersedes delta1:** delta1's fully-symbolic SymPy approach is intractable — SymPy cannot evaluate `∫₀¹ e^{-Pi_star·x}·cos(πx/2)·xⁿ dx` (or the `cosh` analogue) with a SYMBOLIC `Pi_star`; it hangs (>30 min, killed twice). This was empirically confirmed three ways (definite integrate, indefinite antiderivative + bound-eval, trig→exp rewrite). With a CONCRETE rational `Pi_star` the same integrals are exact in ~5-10 s. The Mathematica engine already does the full all-`Pi_star` symbolic proof and is correct/independent (prior audit). So the agreed resolution (orchestrator + Codex consult, session 019e7216-165e-7d90-8900-9cda8b87b513) is:

> **Mathematica** = full symbolic verifier (all `Pi_star`). **SymPy** = an EXACT, symbolic-in-`(r1,r2,A_T,B_T,gprime)`, **multi-point** cross-check at concrete rational `Pi_star` values. This is exact (no floats, no tolerance) and multi-point — materially stronger than the original single mpmath spot-check.

## F1-fix — rewrite the SymPy script as an exact multi-point cross-check

Rewrite `scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py` so it loops over several concrete rational `Pi_star` samples and asserts M1–M7 reduce to **exact 0** (as expressions in the still-symbolic `r1,r2,A_T,B_T,gprime`) at each sample. Do NOT touch the Mathematica script (it is correct and is the full-symbolic authority).

**Required structure:**

1. Keep `r1, r2` real, `gprime` real & nonzero, `A_T, B_T` real, `eps` real, `x` real — all SYMBOLIC. `k = sp.pi/2`.
2. Sample set (exact rationals, all > 0): `PI_SAMPLES = [sp.Rational(1,2), sp.Integer(1), sp.Rational(3,2), sp.Integer(2), sp.Rational(5,3)]`.
3. Exact-zero helper (NO `is_zero`, NO tolerance):
   ```python
   def expect_zero(name, expr):
       res = sp.simplify(sp.cancel(sp.together(expr)))
       print(f"{name} = {res}")
       assert res == 0, f"{name} nonzero: {res}"
   ```
   (If `simplify` ever leaves a non-bare-0 that is genuinely zero, fall back to `sp.factor`/`sp.expand`; do not use numeric tolerance and do not rely on `.is_zero` which can return `None`.)
4. For EACH `Pi_star` in `PI_SAMPLES`, build everything FRESH and INDEPENDENTLY (do not reuse a symbolic-`Pi_star` object; substitute the concrete value so the x-integrals are tractable):
   - `w0 = sp.exp(-Pi_star*x)`; `pert = r1*x + r2*x**2` (this is `R_residual`).
   - Derive `Sigma_star` and `delta_Sigma` by the eps-first method (tractable): `num = sp.series(w0*sp.exp(-eps*pert), eps, 0, 2).removeO()`; `Z0 = sp.integrate(num.coeff(eps,0),(x,0,1))`; `Z1 = sp.integrate(num.coeff(eps,1),(x,0,1))`; `ser = sp.series(num/(Z0+eps*Z1), eps, 0, 2).removeO()`; `Sigma_star = ser.coeff(eps,0)`; `delta_Sigma = ser.coeff(eps,1)`. (With concrete `Pi_star` every integral here is elementary and fast.)
   - Kernels: `c_kernel = sp.cos(k*x)`; `K_kernel = sp.cosh(k*(1-x))/sp.cosh(k)`.
   - `mean = lambda f: sp.integrate(Sigma_star*f, (x,0,1))`; `Rbar=mean(pert)`, `cbar=mean(c_kernel)`, `Kbar=mean(K_kernel)`; `CovcR = mean(c_kernel*pert) - cbar*Rbar`; `CovKR = mean(K_kernel*pert) - Kbar*Rbar`.
   - Assert (label each with the `Pi_star` value so the transcript is self-documenting), all `expect_zero`:
     - **M1**: `delta_Sigma + Sigma_star*(pert - Rbar)`
     - **M2**: `sp.integrate(Sigma_star,(x,0,1)) - 1`
     - **M3**: `sp.integrate(delta_Sigma,(x,0,1))`
     - `dg = sp.integrate(c_kernel*delta_Sigma,(x,0,1))`; `dS = sp.integrate(K_kernel*delta_Sigma,(x,0,1))`
     - **M4**: `dg + CovcR`
     - **M5**: `dS + CovKR`
     - **M6** (δΠ): with `deltaPi = -dg/gprime`, assert `sp.simplify(deltaPi - CovcR/gprime)` is 0 (gprime is a nonzero symbol; this is symbolic division by a symbol, fine). Equivalently you may assert the cleared form `gprime*deltaPi - CovcR == ... ` — either is acceptable as long as it is exact.
     - **M7** (δT_m): with `deltaT = AT*dg + BT*dS`, assert `deltaT + AT*CovcR + BT*CovKR`.
5. Update the module docstring: state plainly that the SymPy engine is an **exact multi-point cross-check** (identities verified exactly at the sampled `Pi_star` values, symbolic in `r1,r2,A_T,B_T,gprime`), and that the **full all-`Pi_star` symbolic proof is carried by the Mathematica engine**. Do NOT describe the SymPy side as "fully symbolic." Remove any leftover `mpmath`/`mp.mpf`/`mp.quad`/`expect_close`/numeric-anchor references.
6. Keep a closing `print("Theorem: ...")` summary if you like.

**Tractability + timeout:** per your standing instructions, RUN the script under the 10-minute cap (`timeout 600 python3 <path>`). Each `Pi_star` sample's integrals are ~5-10 s, so 5 samples should finish in ~1-3 min. If it approaches the cap, you may reduce to 4 samples `{1/2, 1, 3/2, 2}` — but never drop below 4, and never revert to a single point or to floats. Iterate until it exits 0 with all `[Pi=...] M1..M7` lines printing `= 0`.

Append an `## Applied: F1-fix` block when done (note how many `Pi_star` samples and the final runtime).

## Applied: F1-fix

- files_changed:
  - `scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py`
- summary: Rewrote the SymPy audit as an exact five-sample rational `Pi_star` cross-check for M1-M7, symbolic in `r1,r2,A_T,B_T,gprime`; final capped runtime was 126.24 seconds.
- deviation: none
