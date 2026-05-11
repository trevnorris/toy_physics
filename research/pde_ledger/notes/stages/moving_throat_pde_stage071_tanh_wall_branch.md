# Moving-Throat PDE — Stage 071: Canonical Tanh-Wall Branch and Natural Local Mouth Closure

## Purpose

Stage 53 reduced the first explicit moving-throat support branch to generic wall-shape moments `I_f`, `I_g`.

The next honest step is to choose the first canonical wall profile and the first natural local mouth closure.

This stage does that.

The canonical explicit wall is the smooth finite-thickness `tanh` wall,

`f(xi) = (1 + tanh xi)/2,`

so the support profile is

`chi_phi = f'(xi) = (1/2) sech^2 xi.`

For this profile the shell moments are exact:

`I_f = int dxi f'(xi)^2 = 1/3,`

`I_g = int dxi f''(xi)^2 = 4/15,`

hence

`I_g / I_f = 4/5.`

Then a natural local Robin mouth closure is

`K_m = T_X / ell,`

which simply says that the mouth spring is set by the same support scale as the axial wall tension.

With that closure, the whole branch collapses to three parent dimensionless variables.

---

## 1. Canonical `tanh` wall moments

Take

`f(xi) = (1 + tanh xi)/2.`

Then

`f'(xi) = (1/2) sech^2 xi,`

`f''(xi) = - sech^2 xi tanh xi.`

The exact shell moments are

`I_f = int_(−inf)^(+inf) dxi f'(xi)^2 = 1/3,`

`I_g = int_(−inf)^(+inf) dxi f''(xi)^2 = 4/15.`

Therefore

`I_g / I_f = 4/5.`

So the Stage-53 branch formulas become completely explicit.

---

## 2. Explicit branch coefficients

For the canonical `tanh` wall,

`T_X = pi a^2 ell hbar^2 / (3 m rho_w),`

`K_X = 4 pi a^2 ( 5 m^2 c_(s,w)^2 ell^2 + hbar^2 ) / (15 ell m rho_w),`

`J_1 = 1 / (3 H_w),`
with
`H_w = m c_(s,w)^2 / rho_w.`

The exact geometry/support parameter is

`kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L / ell)^2.`

And the wall figure of merit remains

`W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).`

So the shape dependence of `W_wall` still cancels completely, while `kappa` retains one explicit wall-profile number, `4/5`.

---

## 3. Natural local mouth closure

To close the mouth compliance on the same local scale, take

`K_m = T_X / ell.`

This is not a universal theorem of the parent PDE; it is the first natural local Robin closure consistent with the same wall-shell support scale.

Then the Stage-40/41 Robin variable becomes

`eta := K_m L / T_X = L / ell.`

So the explicit branch now has

`eta = L / ell,`

`kappa = 4 (m c_(s,w) L / hbar)^2 + (4/5) (L / ell)^2.`

---

## 4. Three branch control parameters

Define the parent dimensionless ratios

`chi_s := m c_(s,w) L / hbar,`

`Lambda_ell := L / ell,`

`Upsilon_w := 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2).`

Then the entire canonical wall branch reduces to

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,`

`eta = Lambda_ell,`

`W_wall = Upsilon_w Lambda_ell^2.`

This is the main result of the stage.

The first explicit moving-throat branch is no longer described by the seven symbolic inputs `(a, L, ell, T_X, J_1, kappa, eta)`. It is described by the three parent dimensionless controls `(chi_s, Lambda_ell, Upsilon_w)`.

---

## 5. What Stage 54 changes

Stage 53 showed that the explicit branch data were derivable.
Stage 54 shows that, on the first canonical `tanh` wall with the first natural local mouth closure, those data collapse much further than expected.

What remains now is not a symbolic branch ledger.
What remains is an explicit three-parameter branch-placement problem.

That is exactly the right form for the next step, because it means we can now compare the explicit branch directly to the exact Stage-49 / Stage-52 success window instead of still talking in abstract support/source language.
