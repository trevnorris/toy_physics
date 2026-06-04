# Moving-Throat PDE — Stage 077: Shell-Weighted Extraction of `Theta_w` on the Explicit Family-1 Wall

## Purpose

Stage 076 reduced the explicit Family-1 wall-depth datum to

`Theta_w = 25 lambda_mu^2 rho_w^2`

in normalized Family-1 wall variables.

So the only remaining question on this branch is what effective wall density `rho_w` should actually be used on the active shell.

This stage answers that by using the concrete Family-1 radial Thomas–Fermi wall profile already frozen in the coupled-wall appendix together with the canonical wall-support weight used throughout Stages 070–075.

---

## 1. Explicit midplane radial wall profile

On the balanced Family-1 reference branch,

`alpha_r = 10,`

`epsilon_r = 0.05,`

`p_r = 2.`

At the symmetry midplane the endcap term is exponentially negligible, so the coupled full-profile wall reduces to the radial slice

`rho_r(x) = [ 1 - alpha_r S((x-1)/epsilon_r)^2 ]_+^(1/4),`

with

`S(xi) = (1 + tanh xi)/2.`

In the local shell coordinate `xi = (x-1)/epsilon_r`, this becomes

`rho_r(xi) = [ 1 - alpha_r S(xi)^2 ]_+^(1/4).`

The support layer only overlaps the interior flank up to the exact cut point

`xi_* = artanh( 2/sqrt(alpha_r) - 1 ).`

For `alpha_r = 10`,

`xi_* ≈ -0.3855810692.`

So the active support weight lies mostly on the inner edge of the wall, not at the formal center of the soft switch.

---

## 2. Canonical support weight carried from the explicit branch

To stay on the same branch as Stages 070–075, keep the canonical support profile

`chi_phi(xi) = S'(xi) = (1/2) sech^2 xi.`

Its exact normalization moment is

`I_f = int dxi chi_phi(xi)^2 = 1/3.`

So the natural shell-weighted average of any wall quantity `Q(xi)` is

`<Q>_chi := [ int dxi chi_phi(xi)^2 Q(xi) ] / [ int dxi chi_phi(xi)^2 ].`

---

## 3. Why the correct effective wall datum uses `<rho^2>_chi`

Stage 076 showed that the local wall-depth datum scales as

`Theta(xi) = 25 lambda_mu^2 rho_r(xi)^2.`

So the natural effective branch datum is the support-weighted average of `Theta`, not `Theta` evaluated on a point and not `Theta` built from the averaged density after the square is taken.

Therefore

`Theta_w^(chi) = 25 lambda_mu^2 < rho_r^2 >_chi.`

This is the quantity that preserves the actual quadratic wall-depth weighting of the explicit branch.

---

## 4. Numerical extraction on the `alpha_r = 10` Family-1 branch

Using the exact canonical support weight and the explicit Family-1 radial profile gives

`<rho_r>_chi ≈ 0.192619005556493,`

`<rho_r^2>_chi ≈ 0.162745294003265.`

So the effective wall-depth datum is

`Theta_w^(chi) = 25 lambda_mu^2 <rho_r^2>_chi`
`             ≈ 4.06863235008162 lambda_mu^2.`

This is the first concrete branch value of `Theta_w` derived directly from an explicit wall-support profile rather than left symbolic.

---

## 5. Conservative lower-envelope estimator

For bookkeeping it is also useful to record the stricter Jensen-style lower-envelope obtained by averaging the density first and squaring afterward:

`Theta_w^(J) := 25 lambda_mu^2 <rho_r>_chi^2`
`            ≈ 0.927552032539308 lambda_mu^2.`

Because

`<rho_r^2>_chi >= <rho_r>_chi^2,`

this satisfies

`Theta_w^(chi) >= Theta_w^(J).`

So `Theta_w^(J)` can be used as a very conservative branch floor, while `Theta_w^(chi)` is the natural quadratic-energy matching value.

---

## 6. What Stage 077 changes

Before this step, the explicit Family-1 branch still had one unresolved microscopic wall-depth amplitude.

After this step, that datum is no longer abstract. On the canonical explicit branch,

`Theta_w^(chi) ≈ 4.06863235008162 lambda_mu^2`

with the conservative lower envelope

`Theta_w^(J) ≈ 0.927552032539308 lambda_mu^2.`

So the final branch-level question is no longer “what is `Theta_w`?”
It is now only:

> where do these explicit branch values sit relative to the exact Stage-075 threshold window?
