# Moving-Throat PDE — Stage 038: Dimensionless Continuum Placement Map, Exact Product Relation, and the Three-Lane Factorization of the Selected Quadrupole Branch

## Purpose

Stage 20 derived the reduced branch data directly from one explicit finite-throat continuum kernel.
That was already a meaningful advance, but it still left the placement problem written in a somewhat redundant microscopic language.

The next step is therefore to compress the Stage-20 continuum operator into the smallest useful **dimensionless kernel ledger**.

The main result is that the continuum placement problem collapses to five exact kernel ratios:

`eps_eta = c_(eta U)^2 / ( K_U K_eta^(eff) ),`

`eps_W   = c_(UW)^2 sigma / ( K_U K_W^(eff) ),`

`rho     = c_(UW) c_(eta U) / ( K_U c_(eta W) ),`

`Z_W     = c_(eta W)^2 / ( K_eta^(eff) K_W^(eff) ),`

`delta_0 = pi^2 T_w / ( L^2 K_eta^(eff) ).`

In terms of these, the actual selected-branch placement coordinates are

`delta = delta_0 / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ],`

with

`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`

So the continuum PDE-side placement problem is now fully compressed to one geometric anisotropy ratio, one mixed-sector stability ratio, one interference ratio, one wall-to-mixed overlap ratio, and one radiative demand scale.

Even better, these variables satisfy an exact product law,

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`
`               = 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W).`

That means three apparently independent microscopic lanes only redistribute the defect **along** a fixed product curve, while the mixed-sector stability ratio `eps_W` and the radiative demand scale `Lambda` set the product itself.

So Stage 21 turns the Stage-18/19 admissibility problem into a very small continuum-kernel map.

---

## 1. Exact dimensionless kernel ratios

Keep the Stage-20 continuum effective stiffnesses

`K_eta^(eff) = K_eta + 6 T_Omega,`

`K_W^(eff)   = K_W + pi^2 T_W / (4 L^2).`

Now define the exact dimensionless kernel ratios

`eps_eta := c_(eta U)^2 / ( K_U K_eta^(eff) ),`

`eps_W   := c_(UW)^2 sigma / ( K_U K_W^(eff) ),`

`rho     := c_(UW) c_(eta U) / ( K_U c_(eta W) ),`

`Z_W     := c_(eta W)^2 / ( K_eta^(eff) K_W^(eff) ),`

`delta_0 := pi^2 T_w / ( L^2 K_eta^(eff) ),`

and the radiative demand scale

`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`

The exact stability window is now simply

`0 < eps_eta < 1,`

`0 < eps_W   < 1,`

together with the natural nonvanishing transfer branch

`1 + rho != 0.`

---

## 2. Exact continuum placement formulas

Substituting the Stage-20 continuum formulas into the Stage-18/19 branch variables gives

`delta = delta_0 / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ].`

So the selected branch is now placed by only five exact dimensionless kernel ratios.

There is also a useful corollary for the outgoing transfer factor:

`beta_0`
`= (mu_W / mu_eta) [ K_eta^(eff) / K_W^(eff) ]`
`  Z_W (1 + rho)^2 / (1 - eps_W)^2.`

So among the inertial masses only the ratio `mu_W / mu_eta` survives in the conservative-to-outgoing transfer factor, while the geometry coordinate `delta` and the mixed baseline `M_mix` are purely stiffness/coupling ratios.

---

## 3. Exact product relation

Multiplying the exact placement formulas yields the simplest structural identity of the stage:

`R_target M_mix`
`= 8 Lambda (1 - eps_W) / pi^2`
`= 54 G c_s^5 K_W^(eff) (1 - eps_W) / (5 a^5 c^5 mu_W).`

Equivalently,

`R_target M_mix = N_Q^(target) Delta_0 / Omega_U^2.`

This is important because the wall-U dressing `eps_eta`, the wall-to-mixed overlap `Z_W`, and the interference ratio `rho` cancel completely from the product.
Those quantities do still matter, but they only redistribute the physical defect along a fixed product curve.

So the continuum kernel separates into two distinct tasks:

1. the pair `(eps_W, Lambda)` sets the product scale,
2. the trio `(eps_eta, Z_W, rho)` decides where on that product curve the actual defect lands.

That is the cleanest factorization seen so far in the moving-throat branch-placement problem.

---

## 4. Exact one-way parameter tendencies

The dimensionless continuum map also makes the monotone tendencies completely explicit.

### 4.1 Wall-U dressing

`d delta / d eps_eta = delta_0 / (1 - eps_eta)^2 > 0,`

`d M_mix / d eps_eta = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta)^2 (1 - eps_W) ] > 0,`

`d R_target / d eps_eta = - Lambda (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ] < 0.`

So stronger wall-U softening pushes the defect toward larger anisotropy and larger mixed baseline while lowering the normalization demand ratio.

### 4.2 Internal mixed blocking

`d M_mix / d eps_W = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W)^2 ] > 0,`

`d R_target / d eps_W = - 2 Lambda (1 - eps_eta) (1 - eps_W) / [ Z_W (1 + rho)^2 ] < 0.`

So stronger `U/W` blocking simultaneously raises the mixed baseline and lowers the normalization demand ratio.
It is the only stability ratio that also enters the exact product relation.

### 4.3 Wall-to-mixed overlap strength

`d M_mix / d Z_W = 8 (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ] > 0,`

`d R_target / d Z_W = - Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W^2 (1 + rho)^2 ] < 0.`

So increasing the wall-to-mixed overlap pushes the physical point upward in `M_mix` and downward in `R_target`.

### 4.4 Interference ratio

`d M_mix / d rho = 16 Z_W (1 + rho) / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`d R_target / d rho = - 2 Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^3 ].`

So on the natural nonvanishing transfer branch `1 + rho > 0`, constructive interference increases `M_mix` and decreases `R_target`, while destructive interference does the opposite.

---

## 5. Three-lane factorization of the continuum placement problem

The continuum map now splits cleanly into three structural lanes.

### 5.1 Geometry lane

`delta = delta_0 / (1 - eps_eta)`

depends only on the bare wall anisotropy and the wall-U softening ratio.
So the bare geometry lane is controlled entirely by the wall sector plus the brane-like internal doublet.

### 5.2 Mixed-stability/product lane

`R_target M_mix = 8 Lambda (1 - eps_W) / pi^2`

depends only on the radiative demand scale and the internal mixed-sector stability ratio.
So the overall product scale is set in the mixed lane, not by the wall overlap or interference bookkeeping.

### 5.3 Redistribution lane

At fixed product, the trio `(eps_eta, Z_W, rho)` redistributes the point between `R_target` and `M_mix`.
So these parameters control branch placement **along** the product curve but not the curve itself.

This is the sharpest structural decomposition reached so far.

---

## 6. Best current theorem gate after Stage 21

The moving-throat selected quadrupole branch is now described by two exact layers.

### Layer 1 — universal branch geometry

From Stages 18–19, the stable selected branch is controlled by the universal functions

`R_target = F(xi,delta),`

`M_mix <= G(xi,delta).`

### Layer 2 — continuum placement map

From Stages 20–21, the continuum operator places the actual defect at

`delta = delta_0 / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + rho)^2 / [ pi^2 (1 - eps_eta) (1 - eps_W) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps_W)^2 / [ Z_W (1 + rho)^2 ].`

So the remaining theorem gap is now very narrow.

It is no longer:

- “derive more microscopic couplings somehow,”
- or “guess which branch the PDE lands on.”

It is:

> compute the dimensionless kernel ratios `(eps_eta, eps_W, rho, Z_W, delta_0, Lambda)` from the completed moving-throat PDE and check whether the resulting point lies inside the exact Stage-18/19 admissible region.
