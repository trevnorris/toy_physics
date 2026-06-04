# Moving-Throat PDE — Stage 078: Explicit Family-1 Branch Comparison and Closing Verdict for This Subprogram

## Purpose

Stage 075 reduced the explicit Family-1 reference branch to the threshold window

`Theta_fail ≈ 3.62605617972939e-4 * Pe_req,`

`Theta_suff ≈ 4.21495341569977e-2 * Pe_req.`

Stage 077 then extracted the actual explicit branch datum

`Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2,`

with conservative lower envelope

`Theta_w^(J) ≈ 0.927552032539308 lambda_mu^2.`

This stage compares them directly and states the branch-level verdict.

---

## 1. Natural quadratic branch comparison

Use the natural support-weighted datum `Theta_w^(chi)`.

The exact explicit branch inequalities are

`Theta_w >= Theta_suff`  -> guaranteed success,

`Theta_w <= Theta_fail`  -> guaranteed failure.

Insert `Theta_w = Theta_w^(chi)`.

Then the explicit Family-1 branch is guaranteed to succeed whenever

`Pe_req <= Pe_suff^(chi) := Theta_w^(chi) / 4.21495341569977e-2`
`                     ≈ 96.5285247264386 lambda_mu^2,`

and it is guaranteed to fail only if

`Pe_req >= Pe_fail^(chi) := Theta_w^(chi) / 3.62605617972939e-4`
`                     ≈ 11220.5441626259 lambda_mu^2.`

So on the natural explicit wall-depth extraction, the reference Family-1 branch already clears the branch theorem unless the required quadrupole demand is anomalously large.

---

## 2. Conservative lower-envelope comparison

If one insists on the stricter lower-envelope estimator `Theta_w^(J)`, then guaranteed success still holds whenever

`Pe_req <= Pe_suff^(J) := Theta_w^(J) / 4.21495341569977e-2`
`                   ≈ 22.0062226330754 lambda_mu^2,`

while guaranteed failure would require

`Pe_req >= Pe_fail^(J) := Theta_w^(J) / 3.62605617972939e-4`
`                   ≈ 2558.01892349205 lambda_mu^2.`

So even the conservative floor says the explicit Family-1 branch is not obviously wall-depth starved unless the required quadrupole demand is already very large.

---

## 3. Branch-level verdict

The explicit Family-1 subprogram can now be closed cleanly.

The dominant result is:

> on the first concrete moving-throat branch, the wall-depth datum is no longer the natural bottleneck.

Within the natural `n=5` enthalpy lock and the explicit shell-weighted extraction, the branch supplies an `O(1)` wall-depth strength,

`Theta_w^(chi) ≈ 4.069 lambda_mu^2,`

while the Stage-075 success window only demands

`Theta_w >= 4.2149534e-2 Pe_req.`

So the branch succeeds automatically for modest required demand and fails only for very large demand.

That means the explicit-branch phase has reached its natural finish line:

- the geometry is fixed,
- the healing/support scale is fixed,
- the operator thresholds are fixed,
- the wall-depth datum has been extracted from a concrete parent-wall closure,
- and the remaining unresolved quantity is now the quadrupole-side demand `Pe_req`, not the wall-depth supply.

---

## 4. What remains open after Stage 078

This closes the Family-1 explicit-branch placement subprogram.

The remaining serious problem is now the one the higher-order stack had already isolated:

> determine the actual quadrupole-side requirement `Pe_req` from the completed moving-throat / outgoing-normalization branch and compare it to the explicit success/failure bands above.

So from this point on, the real bottleneck is no longer the wall-depth amplitude. It is the final quadrupole normalization bridge.
