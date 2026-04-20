# Moving-Throat PDE — Stage 77: Isotropic Geometry-Decoupling Theorem

## Purpose

Stage 75 showed that the only remaining reduced obstruction to the `3/4 + 1/4` conservative quadrupole module is dynamic contamination from the geometry lane through the two numbers

`eps_2 = Omega_Q^2 K_(g,2) / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole`.

So the next honest check is not broad anymore:

> on the actual **isotropic** moving-throat branch, can the geometry lane contribute nonzero `omega^2` or `omega^4` moments to the grouped real `P2` conservative module at linear order?

This stage answers that directly.

The main result is:

- on an isotropic reference throat, the scalar/geometry `l=0` lane and the grouped real `l=2` lanes are exactly orthogonal in the quadratic wall theory,
- therefore the dynamic `l=0` geometry lane cannot contribute to the grouped real `P2` conservative module at linear order,
- and the only geometry contribution allowed on that branch is the already-identified **static** completion.

Equivalently,

`K_(g,2) = 0`,

`K_(g,4) = 0`,

so

`eps_2 = eps_4 = 0`

on the natural isotropic branch.

---

## 1. Quadratic wall action on the isotropic reference throat

From the distributed wall lift, the quadratic geometry action is

`S_eta^(2) = (1/2) int dt dw dOmega [`
`              mu_eta(w) (partial_t eta)^2`
`              - T_w(w) (partial_w eta)^2`
`              - T_Omega(w) eta (-Delta_(S^2)) eta`
`              - K_eta(w) eta^2 ]`.

All coefficients depend only on the axial coordinate `w`, not on the angles. So the angular operator is exactly `O(3)`-invariant.

Expand the wall field into scalar and quadrupole pieces,

`eta(Omega,w,t) = g(w,t) Y_00(Omega) + sum_A q_A(w,t) Y_(2A)(Omega) + ...`,

where `A in {20,21c,21s,22c,22s}`.

Because the operator is isotropic, every bilinear cross term between `Y_00` and any `Y_(2A)` is proportional to one of the angular integrals

`int Y_00 Y_(2A) dOmega`,

`int grad_(S^2) Y_00 . grad_(S^2) Y_(2A) dOmega`,

or equivalently

`int Y_00 (-Delta_(S^2)) Y_(2A) dOmega`.

All of them vanish exactly:

- `Y_00` is orthogonal to every `l=2` harmonic,
- `grad_(S^2) Y_00 = 0` because `Y_00` is constant,
- and `(-Delta_(S^2)) Y_(2A) = 6 Y_(2A)` reduces the third integral back to the first.

So the quadratic wall action is block diagonal in the `(l=0)` and `(l=2)` sectors.

---

## 2. Exact block structure of the isotropic reduced wall theory

After angular integration, the quadratic reduced action takes the schematic form

`S_red^(2) = (1/2) int dt dw [ g D_0 g + sum_A q_A D_2 q_A ]`,

with no bilinear mixing term of the form `g M_A q_A`.

So, on the isotropic branch,

`M_(0<->2) = 0`.

That means the scalar/geometry lane can renormalize only the scalar sector, while the grouped real `P2` channels evolve independently through their own isotropic quadrupole operator.

This is stronger than a small-coupling statement. It is an exact linear selection rule of the isotropic quadratic wall theory.

---

## 3. Consequence for the conservative grouped-`P2` quadrupole module

Now interpret the 3PN static geometry completion exactly the way the 3PN result says it should be interpreted:

- the grouped real `P2` lane carries the dynamic quadrupole pole structure,
- the leftover geometry completion is a scalar/pair-side static remainder.

Because the isotropic quadratic wall theory has no `l=0 <-> l=2` bilinear mixing, the scalar geometry lane cannot feed any dynamic even moments into the isotropic grouped-`P2` conservative quadrupole module.

So on the natural isotropic branch the effective conservative quadrupole module has the form

`K_Q^cons(omega) = K_(g,0) + K_pole /(1 - omega^2/Omega_Q^2)`

with

`K_(g,2) = 0`,

`K_(g,4) = 0`.

Therefore the contamination numbers of Stage 75 vanish:

`eps_2 = Omega_Q^2 K_(g,2) / K_pole = 0`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole = 0`.

And the exact obstruction formula collapses back to

`c_pole = 1/4`,

`c_geom = 3/4`.

So the `3/4 + 1/4` conservative quadrupole module is not obstructed on the natural isotropic branch.

---

## 4. What would be required to make `eps_2` or `eps_4` nonzero?

This theorem also makes the failure channels precise.

Nonzero geometry contamination at `O(omega^2)` or `O(omega^4)` requires at least one of the following:

1. **explicit angular anisotropy** in the quadratic wall operator, so that `l=0` and `l=2` cease to be orthogonal,
2. **a genuine second dynamic `l=2` geometry pole** independent of the grouped-`P2` pole already being used as the quadrupole carrier,
3. **nonlinear/higher-order backreaction** beyond the linear isotropic branch, which can induce contamination only beyond the exact linear theorem.

None of those are present in the natural minimal isotropic branch frozen by the present reduced hierarchy.

So the actual check requested after Stage 76 comes out cleanly:

> on the isotropic linear moving-throat branch, the geometry lane is dynamically inert through `O(omega^4)` with respect to the grouped real `P2` conservative quadrupole module.

---

## 5. Best current theorem statement after Stage 77

Inside the present reduced hierarchy,

- isotropy makes the quadratic wall operator block diagonal in `l`,
- the scalar/geometry lane is `l=0`,
- the conservative quadrupole carrier is the grouped real `l=2` bundle,
- so the geometry lane contributes only the already-known static completion and no dynamic `omega^2` or `omega^4` contamination.

Therefore

`eps_2 = eps_4 = 0`

on the natural isotropic branch,

and the Stage-74 `3/4 + 1/4` split is recovered exactly.
