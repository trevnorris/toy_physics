# Moving-Throat PDE — Stage 58: Explicit Family-1 Threshold Window and the Last Remaining Wall-Amplitude Datum

## Purpose

After Stages 56–57, the first explicit Family-1 throat-support branch is no longer symbolic in its geometry/support coordinates. It now has

`Lambda_ell = 37,`
`eta = 37,`
`chi_s = 37/2,`
`kappa = 12321/5.`

So the exact Stage-41/42 threshold machinery can now be evaluated numerically on a concrete branch.

This stage does two things:

1. it computes the actual operator thresholds on that explicit branch;
2. it shows that the only remaining unknown is one wall-depth amplitude datum, not the full triplet `(chi_s, Lambda_ell, Upsilon_w)`.

---

## 1. Exact threshold functions on the explicit branch

Stage 55 gave

`Upsilon_fail = Pe_req / [ Lambda_ell^2 Delta_inf(kappa,eta) ],`

`Upsilon_suff = Pe_req / [ Lambda_ell^2 Delta_0(kappa,eta) ].`

On the Family-1/healing-locked branch,

`Lambda_ell = 37,`

`eta = 37,`

`kappa = 12321/5.`

Therefore the exact kernel scales are

`Delta_0(12321/5,37) ≈ 1.73302079021525e-4,`

`Delta_inf(12321/5,37) ≈ 2.01447565540522e-2.`

So the explicit branch thresholds become

`Upsilon_fail ≈ 0.0362605617972939 * Pe_req,`

`Upsilon_suff ≈ 4.21495341569977 * Pe_req.`

Equivalently, the Stage-41 fixed-point coupling thresholds are

`Xi_fail = Pe_req / Delta_inf ≈ 49.6407091004953 * Pe_req,`

`Xi_suff = Pe_req / Delta_0 ≈ 5770.27122609299 * Pe_req.`

So the first explicit throat-support branch has a very wide "indeterminate" operator window in `Xi`, but a concrete and easily stated wall-amplitude window in `Upsilon_w`.

---

## 2. Large-`alpha` interpretation of the reference branch

On this branch,

`alpha := sqrt(kappa) = 111/sqrt(5) ≈ 49.6407091,`

which is already deep in the large-`alpha` regime.

The exact formulas then behave as

`Delta_inf ~ 1/alpha,`

`Delta_0 ~ eta / [ alpha^2 (alpha + eta) ],`

up to exponentially small corrections.

Numerically this is why

`Xi_fail ≈ alpha,`

while `Xi_suff` is much larger.

So the explicit reference throat is a strongly stiff / strongly localized support branch.

---

## 3. The wall-depth amplitude reduction

The remaining explicit branch unknown is `Upsilon_w`.

But on the actual Family-1 wall branch the radial soft-wall profile already carries a fixed dimensionless depth parameter

`alpha_r = 10.`

If the physical radial wall amplitude is written as

`V0 = alpha_r mu_*`

relative to the local Thomas–Fermi enthalpy/chemical-potential scale `mu_*`, then

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2)`
`          = alpha_r^2 Theta_w,`

where the only remaining microscopic amplitude datum is

`Theta_w := 4 rho_w^2 mu_*^2 / (hbar^2 c_(s,w)^2).`

So on the balanced Family-1 reference branch,

`Upsilon_w = 100 Theta_w.`

This is the sharpest explicit reduction of the wall branch so far: after fixing the geometry/support coordinates, the branch no longer depends on a free wall amplitude and a free wall depth separately. It depends on one dimensionless microscopic wall-depth datum `Theta_w`.

---

## 4. Explicit threshold window for `Theta_w`

Since `Upsilon_w = 100 Theta_w`, the explicit branch theorem becomes

`Theta_w <= Theta_fail`  -> fail,

`Theta_w >= Theta_suff`  -> succeed,

with

`Theta_fail = Upsilon_fail / 100`
`           ≈ 3.62605617972939e-4 * Pe_req,`

`Theta_suff = Upsilon_suff / 100`
`           ≈ 4.21495341569977e-2 * Pe_req.`

So the moving-throat placement problem is now no longer a three-parameter branch hunt.

It is one explicit microscopic question:

> does the actual Family-1 wall depth datum `Theta_w` on the true branch lie below, within, or above this window?

---

## 5. What Stage 58 changes

Before this stage, the explicit branch still looked as if it required solving for all three Stage-55 controls.

After this stage:

- the geometry ratio is fixed,
- the support/healing ratio is fixed,
- the support operator scales are fixed,
- and the only remaining branch datum is the microscopic wall-depth amplitude `Theta_w`.

That means the explicit-branch phase is now essentially finished.

What remains is one last microscopic closure question on this branch:
derive or estimate `Theta_w` from the real wall/throat PDE and compare it directly to the explicit threshold window above.
