# Moving-Throat PDE — Stage 62: Exact Family-1 Map from Quadrupole Demand `zeta_req` to the Required Transport Bias `Pe_req`

## Purpose

Stage 61 closed the explicit Family-1 wall-depth subprogram and showed that the natural branch is not wall-depth starved. The remaining serious bottleneck was therefore no longer the supply side `Theta_w`, but the demand side `Pe_req` coming from the quadrupole-normalization branch.

The next honest step is to eliminate the abstract `Pe_req` in favor of the explicit support-ratio demand `zeta_req` carried from Stages 35 and 42.

This stage does that for the concrete Family-1 branch.

The main result is that, on the explicit Family-1 geometry/support branch,

`zeta_F1(Pe) = A_F1 Omega_Pe^2`

with a fixed branch factor

`A_F1 = (kappa_F1 + pi^2/4)/(kappa_F1 + y_F1^2),`

where

`kappa_F1 = 12321/5,`

`y_F1 tan(y_F1) = 37,  0 < y_F1 < pi/2.`

So the remaining demand problem is now the unique inverse relation

`zeta_F1(Pe_req) = zeta_req.`

That already yields a hard Family-1 quadrupole-demand ceiling:

`zeta_req <= zeta_max^(F1) := A_F1 pi^2 / 4`

with

`zeta_max^(F1) ≈ 2.46752922945601.`

So even before the final moving-throat quadrupole normalization is solved, the explicit Family-1 branch can only support the selected quadrupole branch if its required support ratio stays below this concrete `O(1)` ceiling.

---

## 1. Fixed Family-1 support-compliance factor

Stages 56–58 fixed the explicit Family-1 branch to

`kappa_F1 = 12321/5,`

`eta_F1 = 37.`

Let `y_F1` denote the unique Robin root on the constructive branch,

`y_F1 tan(y_F1) = 37,`

`0 < y_F1 < pi/2.`

Then the exact Robin support factor carried from Stage 40 is

`A_F1 := (kappa_F1 + pi^2/4) / (kappa_F1 + y_F1^2).`

Numerically,

`y_F1 ≈ 1.52948248371470,`

`A_F1 ≈ 1.00005192880220.`

So the explicit Family-1 branch is already slightly above the symmetric-twin baseline even at zero transport bias.

---

## 2. Exact Family-1 support-ratio map

Stage 39 gave the exact constructive transport overlap factor

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2) (exp(Pe)-1) ],`

with continuous extension `Omega_0 = 1`.

Therefore the explicit Family-1 lowest-lane support ratio is

`zeta_F1(Pe) = A_F1 Omega_Pe^2.`

This is the exact demand-side bridge that Stage 61 was still missing.

Its endpoint values are

`zeta_F1(0) = A_F1 ≈ 1.00005192880220,`

`lim_(Pe -> +infinity) zeta_F1(Pe) = A_F1 pi^2/4`
`                                  ≈ 2.46752922945601.`

So the full Family-1 constructive branch spans only the compact interval

`1.0000519288... <= zeta_F1(Pe) <= 2.4675292294... .`

That is the first exact explicit ceiling on the quadrupole-demand side.

---

## 3. Exact inversion problem for `Pe_req`

Given a selected quadrupole-branch demand `zeta_req`, the required transport bias is the unique constructive root of

`A_F1 Omega_(Pe_req)^2 = zeta_req.`

Equivalently,

`Omega_(Pe_req)^2 = zeta_req / A_F1.`

Because the carried constructive branch has `Omega_Pe` strictly increasing from `1` to `pi/2`, this inverse exists iff

`A_F1 <= zeta_req <= A_F1 pi^2/4.`

So for the explicit Family-1 branch:

- if `zeta_req < A_F1`, the demand is already met at zero transport bias;
- if `A_F1 <= zeta_req <= zeta_max^(F1)`, there is a unique constructive `Pe_req`;
- if `zeta_req > zeta_max^(F1)`, the Family-1 branch fails irrespective of wall depth.

This is a theorem-level sharpening of Stage 61: the wall-depth supply may be generous, but the branch still has a hard support-ratio ceiling.

---

## 4. Small-demand expansion

Using the carried Stage-39 expansion

`Omega_Pe = 1 + ((4-pi)/(2pi)) Pe + O(Pe^2),`

the Family-1 support-ratio map begins as

`zeta_F1(Pe)`
`= A_F1 [ 1 + ((4-pi)/pi) Pe + O(Pe^2) ].`

So for demands only slightly above the zero-bias baseline,

`zeta_req = A_F1 (1 + eps_z),`

the required transport bias is approximately

`Pe_req ~= (pi/(4-pi)) eps_z.`

This shows that the initial transport cost is linear in the excess quadrupole demand above the zero-bias Family-1 baseline.

---

## 5. What Stage 62 changes

Before this step, the remaining bottleneck was written as the abstract Peclet requirement `Pe_req`.

After this step, the explicit Family-1 branch no longer needs `Pe_req` as an independent concept.
It has an exact demand map

`zeta_req <-> Pe_req`

with a hard ceiling

`zeta_max^(F1) ≈ 2.46752922945601.`

So the remaining serious question is now much narrower:

> does the final selected quadrupole branch demand a support ratio `zeta_req` that stays below this explicit Family-1 ceiling?
