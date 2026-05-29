---
unit_id: 151
batch: IV.6
iteration: 2
created_at: 2026-05-28T21:30:00-06:00
parent_directive: stage_151.md
findings_count: 1
stop_cold: null
applied: false
verification_status: pending
needs_user_resolution: false
---

# Codex delta-directive — unit 151 (iteration 2)

Iteration 1 rewrote the SymPy script to be symbolic (good intent, M1–M7 in place), BUT it **hangs**: line 46 `Z = sp.integrate(sp.exp(-Phi), (x, 0, 1))` integrates `exp` of a *quadratic-in-x with `eps` inside* (`Phi = Pi_star*x + eps*(r1*x + r2*x**2)`). SymPy turns that into `erf`/`erfi` and then chokes when it `sp.series`-expands and simplifies it. (`exec-sympy 151` ran > 30 min with no output and had to be killed.)

The Mathematica side avoids this by expanding in `eps` first; the SymPy side must do the same.

## F1-fix — expand in epsilon BEFORE integrating (make it tractable)

Replace the derivation block (the 6 lines from `Phi = Pi_star * x + ...` through `delta_Sigma = series.coeff(eps, 1)`, currently lines 45–50) with the block below. Everything else in the file (the kernels, `mean`, the covariance definitions, and all M1–M7 `expect_zero` checks at lines 52 onward) stays EXACTLY as is — it already works once `Sigma_star` and `delta_Sigma` are defined from tractable, elementary integrals.

BEFORE (lines 45–50):
```python
Phi = Pi_star * x + eps * (r1 * x + r2 * x**2)
Z = sp.integrate(sp.exp(-Phi), (x, 0, 1))
Sigma_full = sp.exp(-Phi) / Z
series = sp.series(Sigma_full, eps, 0, 2).removeO()
Sigma_star = series.coeff(eps, 0)
delta_Sigma = series.coeff(eps, 1)
```

AFTER:
```python
# Expand in epsilon BEFORE integrating so SymPy only ever integrates
# polynomial * exp(-Pi_star*x) (elementary, fast). Integrating exp of the full
# quadratic exponent directly yields erf/erfi and is intractable here.
pert = r1 * x + r2 * x**2                       # O(eps^1) exponent perturbation (= R_residual)
w0 = sp.exp(-Pi_star * x)                        # zeroth-order (eps=0) weight
# Numerator exp(-Phi) = w0 * exp(-eps*pert); expand the eps-dependence to O(eps):
num = sp.series(w0 * sp.exp(-eps * pert), eps, 0, 2).removeO()
# Normalization Z = int_0^1 num dx, taken order-by-order (each integral elementary):
Z0 = sp.integrate(num.coeff(eps, 0), (x, 0, 1))
Z1 = sp.integrate(num.coeff(eps, 1), (x, 0, 1))
# Sigma_full = num / (Z0 + eps*Z1); expand to O(eps) (rational in eps, no integration):
series = sp.series(num / (Z0 + eps * Z1), eps, 0, 2).removeO()
Sigma_star = sp.simplify(series.coeff(eps, 0))   # = w0/Z0
delta_Sigma = sp.simplify(series.coeff(eps, 1))  # genuinely derived first-order term
```

Why this is still a faithful (non-tautological) derivation of M1: `delta_Sigma` is obtained from the first-order `eps`-coefficient of the *actual normalized source* `num/(Z0+eps*Z1)`, NOT hand-assumed to be `-Sigma_star*(R-<R>)`. M1 then checks that this independently-derived `delta_Sigma` equals the centered hand form — it would fail on any sign/normalization error. (Algebraically `Z1/Z0 = -<R>_*`, so `delta_Sigma = -Sigma_star*(R_residual - Rbar)` and M1 reduces to 0.)

IMPORTANT — per your standing instructions (codex.md): after editing, **RUN** `python3 scripts/moving_throat_pde_stage151_first_order_selected_correction_sympy_audit.py` and confirm it completes quickly (seconds, not minutes) and exits 0 with M1–M7 all printing `= 0`/PASS. If it still hangs or any Mi is nonzero, diagnose and iterate until it exits 0. Do not touch the Mathematica script (it is already correct and independent). Do not reintroduce `mpmath`/`expect_close`/numeric anchors. Append an `## Applied: F1-fix` block when done.
