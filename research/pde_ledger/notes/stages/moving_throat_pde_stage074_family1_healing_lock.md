# Moving-Throat PDE — Stage 074: Healing-Length Lock and the Actual Reference-Branch Support Scale

## Purpose

Stage 56 fixed the first explicit Family-1 support-geometry ratio to

`Lambda_ell = 37,`

and therefore

`eta = 37.`

The next honest step is to stop treating the support/healing variable

`chi_s = m c_(s,w) L / hbar`

as independent.

A useful exact carry-forward fact already exists in the GNLS static compliance sector: the conservative scalar response obeys a Yukawa/Helmholtz law with

`ell_h^2 = hbar^2 / (4 m^2 c_s^2).`

If the active wall width on the explicit throat-support branch is identified with that same local conservative healing width on the support layer, then the branch acquires an exact support lock.

This stage does that.

---

## 1. Exact GNLS healing-width identity

The exact conservative static GNLS compliance law gives

`(1 - ell_h^2 nabla^2) eta = s / (rho_0 m c_s^2),`

with

`ell_h^2 = hbar^2 / (4 m^2 c_s^2).`

Equivalently,

`ell_h = hbar / (2 m c_s).`

This is the static conservative healing/compliance width of the GNLS medium.

---

## 2. Controlled wall-healing closure on the explicit throat branch

The canonical wall branch already uses a thin active shell of width `ell`.

The natural next closure is to identify that active support width with the local GNLS healing width on the wall-support layer:

`ell = ell_h = hbar / (2 m c_(s,w)).`

This is a **controlled local closure**, not yet a theorem of the full moving-throat PDE, but it is the first honest way to tie the support scale to the parent GNLS medium instead of leaving it free.

---

## 3. Exact support-scale lock

Using

`chi_s = m c_(s,w) L / hbar`

and

`ell = hbar / (2 m c_(s,w)),`

one gets

`chi_s = L / (2 ell) = Lambda_ell / 2.`

Since Stage 56 fixed `Lambda_ell = 37`, the same branch now fixes

`chi_s = 37/2 = 18.5.`

So the first explicit throat-support branch no longer has an independent support scale.

---

## 4. Exact `kappa` on the reference branch

Stage 54 gave

`kappa = 4 chi_s^2 + (4/5) Lambda_ell^2.`

Insert `chi_s = Lambda_ell/2`:

`kappa = 4 (Lambda_ell^2 / 4) + (4/5) Lambda_ell^2`
`      = (1 + 4/5) Lambda_ell^2`
`      = (9/5) Lambda_ell^2.`

With `Lambda_ell = 37`,

`kappa = (9/5) * 37^2 = 12321/5 = 2464.2.`

So the explicit branch now has

`chi_s = 37/2,`

`eta = 37,`

`kappa = 12321/5.`

This is the first point where the reference throat-support branch becomes a concrete point in `(kappa, eta)` space rather than a symbolic family.

---

## 5. Useful derived scale

Write

`alpha := sqrt(kappa).`

Then on the same branch

`alpha = sqrt(12321/5) = 179/sqrt(5) ≈ 49.6407091.`

This is the support-decay scale entering the exact Stage-41/42 kernel formulas.

---

## 6. What Stage 57 changes

After Stages 56–57, the explicit Family-1 throat-support branch has fixed

`Lambda_ell = 37,`
`eta = 37,`
`chi_s = 37/2,`
`kappa = 12321/5.`

So the only remaining explicit-branch unknown from the Stage-55 triplet is the wall-loading amplitude `Upsilon_w`.

That is a much sharper endpoint than the earlier symbolic branch map.
