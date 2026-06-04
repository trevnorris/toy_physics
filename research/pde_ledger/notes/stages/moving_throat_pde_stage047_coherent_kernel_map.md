# Moving-Throat PDE — Stage 047: Exact Coherent-Kernel Dimensionless Map and the Support-Enhancement Factor

## Purpose

Stage 046 left the reduced theorem gap in its sharpest form so far:

> the first coherent local D/N kernel lands on the physical tracking branch, but the constructive split-`U` deformation lowers `R_tr` below `1` and therefore worsens the normalization target relative to the old flat branch.

That was already useful, but the actual kernel-strength question was still spread over too many dimensional couplings.

The next honest move is therefore to compress the coherent local kernel to the smallest exact **dimensionless placement map** and identify how the support lane enters the branch.

This stage does that.

The main result is that, on the coherent tracking branch identified in Stage 045, the first coherent local D/N kernel collapses to a very small parameter set:

- one wall/internal dressing ratio `eps_eta`,
- one split-`U` axial ratio `delta_U`,
- one common interference ratio `chi_0`,
- one mixed blocking ratio `eps_W`,
- one wall-to-mixed overlap ratio `Z_W`,
- one support-to-mixed amplitude ratio `zeta`,
- and one radiative demand scale `Lambda`.

In those variables the coherent branch is fully described by

`R_tr = [ 1 + chi_0/(1+delta_U) ] / (1 + chi_0),`

`eps = eps_W^(split) = eps_W [ 1 - (2/11) delta_U/(1+delta_U) ],`

`delta = [ delta_0 + eps_eta delta_U/(1+delta_U) ] / (1 - eps_eta),`

`M_mix = 8 Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - eps) ],`

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - zeta eps) ],`

`R_target = Lambda (1 - eps_eta) (1 - eps)^2 / [ Z_W (1 + chi_0)^2 ].`

So the entire support lane enters through one exact factor

`S(zeta;eps) := M_tr / M_mix = 1 + zeta(1-eps)/(1-zeta eps),`

with

`M_tr = M_mix S(zeta;eps).`

That turns the support problem into a one-parameter enhancement problem rather than a diffuse multidimensional closure question.

---

## 1. Coherent local interaction data

On the coherent local D/N branch from Stage 045 the mixed and support lanes couple through the same local source density,

`L_int^(coh)`
`= - int_0^L ds [ lambda_W W + lambda_phi phi ] [ eta - gamma U ].`

So the continuum couplings obey

`c_(etaW) = lambda_W,`

`c_(UW)   = gamma lambda_W,`

`c_(etaphi) = lambda_phi,`

`c_(Uphi)   = gamma lambda_phi.`

The wall/internal coupling `c_(etaU)` remains independent.

The split-`U` stiffness is

`K_(U1) = K_U (1 + delta_U),`

with

`delta_U = pi^2 T_U / (L^2 K_U).`

The effective wall and D/N support stiffnesses are

`K_eta^(eff) = K_eta + 6 T_Omega,`

`K_W^(eff)   = K_W   + pi^2 T_W   / (4 L^2),`

`K_phi^(eff) = K_phi + pi^2 T_phi / (4 L^2).`

As before,

`sigma = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2).`

---

## 2. Exact coherent dimensionless ratios

The natural coherent-kernel dimensionless ratios are

`eps_eta := c_(etaU)^2 / ( K_U K_eta^(eff) ),`

`eps_W   := gamma^2 lambda_W^2 sigma / ( K_U K_W^(eff) ),`

`Z_W     := lambda_W^2 / ( K_eta^(eff) K_W^(eff) ),`

`delta_0 := pi^2 T_w / ( L^2 K_eta^(eff) ),`

`chi_0   := gamma c_(etaU) / K_U,`

`zeta    := lambda_phi^2 K_W^(eff) / ( lambda_W^2 K_phi^(eff) ),`

`Lambda  := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W).`

Two exact coherent identities then follow immediately:

`rho_0 = sigma_0 = chi_0,`

`eps_phi = zeta eps_W,`

`Z_phi   = zeta Z_W.`

So the support lane is not independent of the mixed lane at the level of dimensionless couplings. It differs only by the one ratio `zeta`.

---

## 3. Exact split data on the coherent branch

The coherent tracking factor is

`R_tr = [ 1 + chi_0/(1+delta_U) ] / (1 + chi_0),`

with the same exact range as Stage 045:

`1/(1+delta_U) < R_tr < 1`

on the constructive branch `chi_0 > 0`, `delta_U > 0`.

The split mixed and support blocking ratios are

`eps = eps_W^(split)`
`    = eps_W [ 1 - (2/11) delta_U/(1+delta_U) ],`

`eps_phi^(split) = zeta eps.`

The geometry anisotropy is still

`delta = [ delta_0 + eps_eta delta_U/(1+delta_U) ] / (1 - eps_eta).`

So once the coherent kernel is imposed, the branch direction, the blocking, and the anisotropy are all explicit functions of the same small dimensionless set.

---

## 4. Exact mixed/support baselines and support-enhancement factor

The mixed baseline becomes

`M_mix = 8 Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - eps) ].`

Using `Z_phi = zeta Z_W` and `eps_phi^(split) = zeta eps`, the support baseline is

`M_supp = 8 zeta Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - zeta eps) ].`

Therefore the exact total baseline is

`M_tr = M_mix + M_supp`
`     = M_mix S(zeta;eps),`

with the exact support-enhancement factor

`S(zeta;eps)`
`:= 1 + zeta(1-eps)/(1-zeta eps).`

This is the cleanest structural simplification of the stage.

The support lane does not change the radiative demand ratio directly. It only enhances the total baseline by the one factor `S`.

### Exact properties of `S`

`S(0;eps) = 1,`

`dS/dzeta = (1-eps)/(1-zeta eps)^2 > 0`

for `0 < eps < 1` and `0 <= zeta < 1/eps`.

So the coherent support lane is a strictly monotone baseline enhancer.

---

## 5. Exact target ratio and product law on the coherent branch

The normalization-demand ratio remains purely mixed/outgoing:

`R_target = Lambda (1 - eps_eta) (1 - eps)^2 / [ Z_W (1 + chi_0)^2 ].`

So `R_target` is exactly independent of `zeta`.

Multiplying by the total coherent baseline gives

`R_target M_tr`
`= 8 Lambda (1 - eps) / pi^2 * S(zeta;eps).`

Equivalently,

`R_target M_tr`
`= 8 Lambda (1 - eps) / pi^2`
`  * [ 1 + zeta(1-eps)/(1-zeta eps) ].`

So the mixed lane still sets the bare product scale, while the support lane simply multiplies that scale by the one enhancement factor `S`.

This is the exact coherent-kernel replacement of the Stage-038 product law.

---

## 6. Best current theorem statement after Stage 047

The first coherent local D/N kernel is now compressed to an exact dimensionless placement map.

### What is exact now

- the coherent branch depends on one common interference ratio `chi_0`,
- the support lane differs from the mixed lane by only one ratio `zeta`,
- the exact tracking factor is `R_tr(chi_0,delta_U)`,
- the exact mixed blocking is `eps = eps_W^(split)`,
- the exact total baseline is `M_tr = M_mix S(zeta;eps)`,
- the target ratio `R_target` is independent of `zeta`,
- and the support lane is a strictly monotone enhancer of the total baseline.

### What this means physically

The reduced theorem gap is no longer “how do we represent the support lane?”
It is now much sharper:

> for the first coherent local D/N kernel on the coherent tracking branch, the support lane enters only by increasing the baseline through one exact enhancement factor, while the retarded demand ratio is still fixed by the mixed/outgoing lane.

So the next honest question is no longer about closure choice.
It is:

> **how much coherent support enhancement is needed to overcome the exact tracking-branch deficit before the stable branch softens out?**
