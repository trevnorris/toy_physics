# Moving-Throat PDE — Stage 63: Explicit Family-1 Conversion of the Stage-61 `Pe_req` Window into Quadrupole-Demand Thresholds `zeta_req`

## Purpose

Stage 61 reduced the explicit Family-1 branch to a wall-depth comparison window in the transport-bias variable `Pe_req`:

`Pe_req <= Pe_suff^(chi)`  -> guaranteed success on the natural shell-weighted branch,

`Pe_req >= Pe_fail^(chi)`  -> guaranteed failure,

with the conservative lower-envelope pair `Pe_suff^(J)`, `Pe_fail^(J)` defined analogously.

Stage 62 then replaced the abstract `Pe_req` by the exact Family-1 support-ratio map

`zeta_F1(Pe) = A_F1 Omega_Pe^2.`

So the next honest step is to push the explicit Family-1 verdict entirely into the quadrupole-demand variable `zeta_req`.

This stage does that.

The main result is that the Family-1 branch now carries exact quadrupole-demand windows

`zeta_req <= zeta_suff^(chi)(lambda_mu)`  -> guaranteed success,

`zeta_req >= zeta_fail^(chi)(lambda_mu)`  -> guaranteed failure,

with corresponding conservative lower-envelope functions `zeta_suff^(J)`, `zeta_fail^(J)`.

For `lambda_mu = 1`, the natural Family-1 branch already reaches

`zeta_suff^(chi)(1) ≈ 2.46622291347846,`

while the hard Family-1 ceiling is only slightly larger,

`zeta_max^(F1) ≈ 2.46752922945601.`

So on the explicit branch the wall-depth supply is indeed not the dominant unresolved issue. The remaining serious question is whether the final selected quadrupole branch demands a support ratio above or below this narrow `O(2.46)` window.

---

## 1. Exact Stage-61 transport thresholds

Stage 61 gave the explicit Family-1 transport-bias windows

`Pe_suff^(chi) = 96.5285247264386 lambda_mu^2,`

`Pe_fail^(chi) = 11220.5441626259 lambda_mu^2,`

for the natural shell-weighted datum, and

`Pe_suff^(J) = 22.0062226330754 lambda_mu^2,`

`Pe_fail^(J) = 2558.01892349205 lambda_mu^2,`

for the conservative lower envelope.

---

## 2. Exact conversion to `zeta_req` thresholds

Using the explicit Family-1 demand map from Stage 62,

`zeta_F1(Pe) = A_F1 Omega_Pe^2,`

define

`zeta_suff^(chi)(lambda_mu) := zeta_F1(Pe_suff^(chi)),`

`zeta_fail^(chi)(lambda_mu) := zeta_F1(Pe_fail^(chi)),`

and similarly

`zeta_suff^(J)(lambda_mu) := zeta_F1(Pe_suff^(J)),`

`zeta_fail^(J)(lambda_mu) := zeta_F1(Pe_fail^(J)).`

Then the explicit Family-1 branch theorem becomes

`zeta_req <= zeta_suff^(chi)(lambda_mu)`  -> guaranteed success,

`zeta_req >= zeta_fail^(chi)(lambda_mu)`  -> guaranteed failure,

with the conservative version obtained by replacing `(chi)` with `(J)`.

So the Stage-61 wall-depth result is now written entirely in the quadrupole-demand language.

---

## 3. Numerical values at `lambda_mu = 1`

For the natural shell-weighted datum,

`zeta_suff^(chi)(1) ≈ 2.46622291347846,`

`zeta_fail^(chi)(1) ≈ 2.46752913273870.`

For the conservative lower envelope,

`zeta_suff^(J)(1) ≈ 2.44257571477179,`

`zeta_fail^(J)(1) ≈ 2.46752736855058.`

Compare these to the exact Family-1 ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So on the natural explicit branch with `lambda_mu = 1`, the guaranteed-success threshold already lies less than `0.00131` below the hard ceiling, and the guaranteed-failure threshold is essentially saturated.

This is the sharpest numerical statement yet of the Stage-61 verdict.

---

## 4. Large-`lambda_mu` limit

Because all four Stage-61 transport thresholds scale as `lambda_mu^2` and `Omega_Pe` approaches `pi/2`, the corresponding quadrupole-demand thresholds satisfy

`lim_(lambda_mu -> +infinity) zeta_suff^(chi) = zeta_max^(F1),`

`lim_(lambda_mu -> +infinity) zeta_fail^(chi) = zeta_max^(F1),`

and likewise for the `(J)` branch.

So increasing the wall-depth normalization beyond `O(1)` does not open an unlimited quadrupole-demand window.
It only drives the branch toward the same hard Family-1 ceiling found in Stage 62.

This is another precise sense in which wall-depth supply is no longer the dominant open issue.

---

## 5. What Stage 63 changes

Before this step, the explicit Family-1 result still spoke in the transport-bias variable `Pe_req`.

After this step, the explicit branch verdict is fully phrased in the same support-ratio variable `zeta_req` that the quadrupole normalization branch actually demands.

So the remaining serious theorem question is now extremely narrow:

> does the final selected quadrupole branch require `zeta_req` below the explicit Family-1 support window `zeta_suff^(chi)(lambda_mu)` (or at least below the hard ceiling `zeta_max^(F1)`), or does it demand more than the explicit branch can ever supply?
