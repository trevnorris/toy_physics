# Moving-Throat PDE — Stage 091: Deriving the `3/4 + 1/4` Conservative Quadrupole Module from the Grouped-`P2` + Geometry Split

## Purpose

Stage 088 extracted

`rho_alpha = 4/3`,

`zeta_req = 1/3`,

from the minimal isotropic conservative quadrupole module, but it did so by taking the contact-plus-pole representation as the natural reading of the explicit support/source branch.

The next honest step is sharper:

> derive that same `3/4 + 1/4` split directly from the conservative **grouped real `P2` + geometry** organization already frozen by the 3PN program.

This stage does that.

The main result is:

- if the isotropic grouped-`P2` conservative branch is carried by one effective pole,
- and if the 3PN geometry completion is genuinely static through `O(omega^4)`,

then the minimal isotropic 2.5PN branch identity forces

`K_pole = K0/4`,

`K_geom = 3 K0/4`,

and therefore

`Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

So the `3/4 + 1/4` conservative module is no longer just a plausible parametrization. Under the minimal grouped-`P2` + static-geometry realization, it is forced.

---

## 1. Frozen input from the conservative hierarchy

The 3PN conservative split already says that the higher conservative payload is organized as

`free compiler image + grouped real P2 middle block + unique geometry completion`.

In the compact summary notation,

`Delta L3^GR = Delta l1 v^8 + L_(P2)^mid + Delta l15^(g) U^4`.

The important structural point is that the geometry lane appears there as a **static** completion, while the grouped-`P2` middle block carries the nontrivial higher-order conservative quadrupole structure.

The 2.5PN program, independently, already fixed the minimal isotropic quadrupole branch identity

`K0 K4 = 4 K2^2`,

equivalently in normalized language,

`u4 = 4 u2^2`.

So the only remaining move is to combine those two facts in the smallest consistent way.

---

## 2. Minimal grouped-`P2` + geometry realization

Take the isotropic conservative quadrupole module in the form

`K_Q^cons(omega) = K_geom + K_pole /(1 - omega^2/Omega_Q^2)`.

Here:

- `K_geom` is the static geometry completion,
- `K_pole` is the isotropic grouped-`P2` pole residue,
- `Omega_Q` is the effective isotropic grouped-`P2` pole.

This is the smallest realization compatible with the 3PN conservative split if the grouped-`P2` side is the only dynamic quadrupole lane and geometry contributes only the static completion.

Expanding at low frequency gives

`K_Q^cons(omega) = K0 + K2 omega^2 + K4 omega^4 + O(omega^6)`

with exact coefficients

`K0 = K_geom + K_pole`,

`K2 = K_pole / Omega_Q^2`,

`K4 = K_pole / Omega_Q^4`.

---

## 3. The branch identity forces the `3/4 + 1/4` split

Insert those coefficients into the minimal isotropic branch identity

`K0 K4 = 4 K2^2`.

This gives

`(K_geom + K_pole) * (K_pole / Omega_Q^4)`
`= 4 * (K_pole / Omega_Q^2)^2`.

Assuming the branch is nontrivial (`K_pole != 0`), the common factor `K_pole / Omega_Q^4` cancels and one finds

`K_geom + K_pole = 4 K_pole`.

So

`K_geom = 3 K_pole`.

Equivalently,

`K_pole = K0 / 4`,

`K_geom = 3 K0 / 4`.

Therefore the normalized conservative quadrupole response is forced to be

`Yhat_Q^cons(omega)`
`= K_Q^cons(omega) / K0`
`= 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

This is exactly the minimal isotropic conservative quadrupole module previously isolated from the outgoing `l=2` moment matching.

So the earlier `3/4 + 1/4` structure is now recovered directly from the grouped-`P2` + static-geometry split.

---

## 4. Immediate corollary for the support/source loading ratio

Stage 088 already proved that on the explicit support/source branch the static contact fraction and the finite conservative pole map to the loading ratio as

`c0 = alpha_mix / alpha_req`,

`c1 = (alpha_req - alpha_mix)/alpha_req`,

with `c0 + c1 = 1`.

Using the newly derived grouped-`P2` + geometry split,

`c0 = 3/4`,

`c1 = 1/4`.

So

`alpha_mix / alpha_req = 3/4`,

which gives exactly

`rho_alpha = alpha_req / alpha_mix = 4/3`,

and

`zeta_req = (alpha_req - alpha_mix)/alpha_mix = 1/3`.

So the Stage-088 loading-ratio extraction is now a direct corollary of the grouped-`P2` + static-geometry realization.

---

## 5. What this actually proves

This stage does **not** prove that the full moving-throat PDE has already produced the minimal isotropic branch.

What it proves is narrower and more useful:

> if the actual conservative grouped-`P2` branch is one isotropic pole and the geometry lane contributes only the static completion through `O(omega^4)`, then the `3/4 + 1/4` conservative quadrupole module is forced algebraically.

So the remaining reduced theorem gap is now extremely sharp:

- either the real moving-throat branch obeys that minimal grouped-`P2` + static-geometry realization,
- or the missing PDE must generate extra dynamic geometry moments or a richer isotropic grouped-`P2` pole structure.

That is exactly the right question to carry into the next phase.
