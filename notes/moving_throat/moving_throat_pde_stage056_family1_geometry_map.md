# Moving-Throat PDE — Stage 56: Family-1 Reference-Branch Geometry Map

## Purpose

Stage 55 reduced the first explicit moving-throat support theorem to the branch variables

`chi_s = m c_(s,w) L / hbar,`

`Lambda_ell = L / ell,`

`Upsilon_w = 4 rho_w^2 V0^2 / (hbar^2 c_(s,w)^2).`

The next honest step is to stop treating the geometry ratio `Lambda_ell` as free and evaluate it on the first concrete throat branch already present in the frozen constructive stack.

The relevant frozen inputs are:

- the carried preferred throat aspect ratio `L/a ≈ 1.85`,
- the balanced thin-layer-consistent Family-1 radial wall branch
  `epsilon_r = 0.05`,
- and the explicit thin-wall identification `ell = epsilon_r a`.

This stage does that map exactly on the chosen reference branch.

---

## 1. Reference-branch wall width

On the radial Family-1 reference branch,

`epsilon_r = 0.05 = 1/20.`

With the radial soft-wall coordinate written as

`xi = (r-a)/ell,`

the natural thin-layer identification is

`ell = epsilon_r a.`

So on the chosen balanced reference branch,

`ell / a = 1/20.`

---

## 2. Carried geometric aspect ratio

The lower-order stack already carries the preferred throat aspect ratio

`L / a ≈ 1.85.`

For the explicit reference branch used here, write that carried value as

`Lambda_* := L/a = 37/20.`

This is a **reference-branch numerical freeze**, not a new theorem of the unsolved moving-throat PDE.

---

## 3. Exact reference-branch value of `Lambda_ell`

Combine

`L/a = 37/20,`

`ell/a = 1/20.`

Then

`Lambda_ell := L/ell = (L/a)/(ell/a) = (37/20)/(1/20) = 37.`

So the first explicit moving-throat branch fixes

`Lambda_ell = 37.`

This is the first truly explicit branch value beyond the symbolic Stage-55 map.

---

## 4. Immediate consequence for the Robin mouth variable

Stage 54 adopted the first natural local mouth closure

`K_m = T_X / ell.`

So

`eta := K_m L / T_X = L/ell = Lambda_ell.`

Therefore the same reference branch fixes

`eta = 37.`

So the first explicit throat-support branch is now pinned to one concrete large-`eta` Robin regime.

---

## 5. What Stage 56 changes

Before this step, the first explicit branch still depended on the symbolic support-geometry ratio `Lambda_ell`.

After this step, the balanced Family-1 / thin-wall / preferred-aspect-ratio branch fixes it to

`Lambda_ell = 37,`

and therefore also fixes the Robin support variable to

`eta = 37.`

That means the remaining actual branch uncertainty is now concentrated in the support/healing scale `chi_s` and the wall-loading amplitude `Upsilon_w`.
