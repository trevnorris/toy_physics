# Moving-Throat PDE — Stage 072: Explicit Branch Placement Map and Threshold Surfaces

## Purpose

Stages 70–71 reduced the first explicit moving-throat support branch to three parent dimensionless variables:

`chi_s = m c_(s,w) L / hbar,`

`Lambda_ell = L / ell,`

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2),`

with

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,`

`eta = Lambda_ell,`

`W_wall = Upsilon_w Lambda_ell^2.`

So the exact Stage-66 / Stage-69 support/source theorem can now be written directly on the first explicit branch.

This stage does that.

---

## 1. Exact explicit-branch fail/succeed surfaces

The universal matched-branch theorem is

`W_wall <= Pe_req / Delta_inf(kappa,eta)`  -> fail,

`W_wall >= Pe_req / Delta_0(kappa,eta)`    -> succeed.

Insert the explicit branch formulas:

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2,`

`eta = Lambda_ell,`

`W_wall = Upsilon_w Lambda_ell^2.`

Then the first explicit moving-throat branch satisfies

`Upsilon_w <= Upsilon_fail(chi_s,Lambda_ell)`  -> fail,

`Upsilon_w >= Upsilon_suff(chi_s,Lambda_ell)`  -> succeed,

where

`Upsilon_fail`
`:= Pe_req / [ Lambda_ell^2 Delta_inf( 4 chi_s^2 + (4/5) Lambda_ell^2, Lambda_ell ) ],`

`Upsilon_suff`
`:= Pe_req / [ Lambda_ell^2 Delta_0( 4 chi_s^2 + (4/5) Lambda_ell^2, Lambda_ell ) ].`

So the whole first explicit moving-throat support/source theorem has now collapsed to a comparison of one physical wall-loading amplitude `Upsilon_w` against two explicit threshold surfaces.

---

## 2. Exact physical wall-amplitude thresholds

Because

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2),`

the explicit wall-amplitude thresholds are

`V0_fail^2`
`= hbar^2 c_(s,w)^2 Upsilon_fail / (4 rho_w^2),`

`V0_suff^2`
`= hbar^2 c_(s,w)^2 Upsilon_suff / (4 rho_w^2).`

So the first explicit branch no longer speaks in terms of abstract gains at all.

It speaks directly in terms of the physical wall amplitude `V0` and the two explicit support-geometry ratios `chi_s`, `Lambda_ell`.

---

## 3. Two asymptotic branch regimes

The explicit branch already shows two physically distinct regimes.

### (a) Shell-gradient dominated branch

If

`(4/5) Lambda_ell^2 >> 4 chi_s^2,`

then

`kappa ~ (4/5) Lambda_ell^2,`

so `alpha = sqrt(kappa) ~ 2 Lambda_ell / sqrt(5)`.

In this regime the threshold surfaces reduce to

`Upsilon_fail ~ 2 Pe_req / (sqrt(5) Lambda_ell),`

`Upsilon_suff ~ (4/5)(1 + 2/sqrt(5)) Pe_req.`

So the fail threshold decreases with increasing shell aspect ratio, while the sufficiency threshold saturates to a finite constant multiple of `Pe_req`.

### (b) Compression dominated branch

If

`4 chi_s^2 >> (4/5) Lambda_ell^2,`

then

`kappa ~ 4 chi_s^2,`

so `alpha ~ 2 chi_s`.

In this regime

`Upsilon_fail ~ 2 Pe_req chi_s / Lambda_ell^2,`

`Upsilon_suff ~ 4 Pe_req chi_s^2 (Lambda_ell + 2 chi_s) / Lambda_ell^3.`

If in addition `chi_s >> Lambda_ell`, this becomes

`Upsilon_suff ~ 8 Pe_req chi_s^3 / Lambda_ell^3.`

So compression-dominated branches are much harder to push across the universal success threshold.

---

## 4. What Stage 072 changes

At the end of the reduced support/source phase, the theorem window still involved the symbolic branch data `(kappa, eta, W_wall)`.

After Stages 70–72, the first explicit moving-throat branch is no longer expressed in those symbols.

It is expressed directly in the parent variables

`(chi_s, Lambda_ell, Upsilon_w)`

or equivalently in the physical branch quantities

- `L/ell`,
- `m c_(s,w) L / hbar`,
- `rho_w`,
- `c_(s,w)`,
- `V0`.

So the next phase is now sharply defined:

> compute the actual moving-throat branch values of `chi_s`, `Lambda_ell`, and `Upsilon_w` on the real throat, and compare them directly to the explicit surfaces `Upsilon_fail`, `Upsilon_suff`.

That is the first genuinely explicit branch-placement problem beyond the finished reduced support/source program.
