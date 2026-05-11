# Moving-Throat PDE — Stage 076: Exact `n=5` Wall-Depth Lock for the Family-1 Branch

## Purpose

Stage 58 reduced the explicit Family-1 reference branch to one remaining microscopic datum,

`Theta_w := 4 rho_w^2 mu_*^2 / (hbar^2 c_(s,w)^2),`

through

`Upsilon_w = alpha_r^2 Theta_w,`

with the frozen reference value `alpha_r = 10`.

The next honest step is to stop treating `Theta_w` as an opaque constant.

The present stage does that by combining two already-carried exact identities:

1. the frozen `n=5` GNLS thermodynamic relation,
2. the Stage-57 healing-width lock.

The result is that on the reference branch the wall-depth datum is not an arbitrary amplitude. It is directly proportional to the square of the effective wall density on the active shell.

---

## 1. Exact `n=5` enthalpy–sound-speed identity

For the frozen GNLS EOS

`P(rho) = K rho^5,`

we have

`c_s^2(rho) = (1/m) dP/drho = (5K/m) rho^4,`

and

`U(rho) = (K/4) rho^5,`

`h(rho) = dU/drho = (5K/4) rho^4.`

Therefore the local enthalpy is exactly

`h(rho) = m c_s(rho)^2 / 4.`

So on the `n=5` branch, enthalpy and sound speed are not independent.

---

## 2. Local enthalpy lock for the Family-1 wall amplitude

Stage 58 wrote the radial wall depth as

`V0 = alpha_r mu_*`

relative to a local Thomas–Fermi enthalpy / chemical-potential scale `mu_*`.

The natural exact closure on the `n=5` branch is to write

`mu_* = lambda_mu h_w,`

where

`h_w := h(rho_w) = m c_(s,w)^2 / 4,`

and `lambda_mu` keeps track of whether one chooses the wall enthalpy itself (`lambda_mu = 1`) or a nearby local chemical-potential normalization.

Then

`Theta_w = 4 rho_w^2 mu_*^2 / (hbar^2 c_(s,w)^2)`
`        = lambda_mu^2 m^2 rho_w^2 c_(s,w)^2 / (4 hbar^2).`

So the only remaining microscopic input is now the effective wall density on the active shell.

---

## 3. Exact healing-lock reduction

Stage 57 already fixed the local healing-width closure

`ell = hbar / (2 m c_(s,w)).`

Insert this into the previous formula:

`m c_(s,w) / hbar = 1 / (2 ell),`

hence

`Theta_w = lambda_mu^2 rho_w^2 / (16 ell^2).`

So on the explicit moving-throat branch,

> the wall-depth datum is exactly the square of the active-shell density, measured in healing-width units, up to the one local normalization factor `lambda_mu`.

---

## 4. Reference-branch form

On the explicit Family-1 reference branch,

`ell / a = 1/20.`

In the dimensionless Family-1 wall coordinates already used by the coupled profile,

`x = r/a,`

so `a = 1` in the normalized wall-shape description and therefore

`1 / ell^2 = 400.`

Hence the branch-level datum becomes

`Theta_w = 25 lambda_mu^2 rho_w^2.`

This is the cleanest algebraic reduction of `Theta_w` so far.

---

## 5. What Stage 59 changes

Before this step, the explicit Family-1 branch still had one unresolved microscopic wall datum with no clean parent formula.

After this step, that datum is no longer opaque:

`Theta_w = 25 lambda_mu^2 rho_w^2`

in the normalized Family-1 wall variables.

So the only real remaining task on this branch is now to choose the correct effective wall density `rho_w` on the active shell and compare the resulting `Theta_w` to the explicit Stage-58 threshold window.
