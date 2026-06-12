# Moving-Throat PDE — Stage 092: Exact Obstruction Formula if the Geometry Lane Carries Dynamic Even Moments

## Purpose

Stage 091 showed that the grouped-`P2` + geometry split forces the `3/4 + 1/4` conservative quadrupole module **if**

- the grouped-`P2` side is one isotropic conservative pole,
- and the geometry lane is static through `O(omega^4)`.

The natural next question is:

> what exactly changes if the geometry lane is **not** purely static, but carries its own `omega^2` and `omega^4` moments?

This stage answers that exactly.

The main result is that the pole fraction is then no longer fixed to `1/4`.
Instead it becomes

`c_pole = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`,

where

`eps_2 = Omega_Q^2 K_(g,2) / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole`.

So the `3/4 + 1/4` split is recovered **iff** the geometry lane is static at the relevant orders.

This is the cleanest reduced obstruction formula for the next phase.

---

## 1. Generalized isotropic grouped-`P2` + geometry ansatz

Allow the geometry lane to carry its own even moments:

`K_g(omega) = K_(g,0) + K_(g,2) omega^2 + K_(g,4) omega^4 + O(omega^6)`.

Keep the grouped-`P2` side as one isotropic conservative pole:

`K_P2(omega) = K_pole /(1 - omega^2/Omega_Q^2)`.

Then the total conservative isotropic quadrupole module is

`K_Q^cons(omega) = K_g(omega) + K_P2(omega)`.

Its low-frequency coefficients are

`K0 = K_(g,0) + K_pole`,

`K2 = K_(g,2) + K_pole/Omega_Q^2`,

`K4 = K_(g,4) + K_pole/Omega_Q^4`.

---

## 2. Exact branch identity with dynamic geometry

Imposing the same minimal isotropic branch identity

`K0 K4 = 4 K2^2`

now gives the exact relation

`(K_(g,0) + K_pole)( K_(g,4) + K_pole/Omega_Q^4 )`
`= 4 ( K_(g,2) + K_pole/Omega_Q^2 )^2.`

So the geometry-contact term is no longer forced to equal `3 K_pole` unless the dynamic geometry moments vanish.

Solving for `K_(g,0)` gives

`K_(g,0)`
`= 4 ( K_(g,2) + K_pole/Omega_Q^2 )^2 / ( K_(g,4) + K_pole/Omega_Q^4 )`
`  - K_pole.`

That is the exact obstruction formula.

---

## 3. Pole fraction in dimensionless contamination variables

Define the dimensionless geometry-contamination parameters

`eps_2 = Omega_Q^2 K_(g,2) / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4) / K_pole`.

Then the total static normalization on the minimal branch is

`K0 = 4 K_pole (1 + eps_2)^2 / (1 + eps_4)`.

So the grouped-`P2` pole fraction becomes

`c_pole = K_pole / K0`
`       = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`.

And therefore the geometry contact fraction is

`c_geom = 1 - c_pole`.

In the strict static-geometry limit

`eps_2 = eps_4 = 0`,

this reduces exactly to

`c_pole = 1/4`,

`c_geom = 3/4`.

So the `3/4 + 1/4` split is not a generic identity of “grouped-`P2` plus geometry” in the abstract.
It is the exact consequence of the **static-geometry** realization.

Provenance note. Stage 091 forces \(K_{\rm geom}=3K_{\rm pole}\) only on the static-geometry slice
\(\epsilon_2=\epsilon_4=0\); this stage deliberately frees the split unless those contamination variables vanish.

---

## 4. Small-contamination expansion

For small geometry contamination,

`|eps_2| << 1`,

`|eps_4| << 1`,

the pole fraction expands as

`c_pole`
`= 1/4 [ 1 + eps_4 - 2 eps_2 + O(eps^2) ]`.

So:

- positive `eps_4` raises the pole fraction,
- positive `eps_2` lowers it twice as strongly at first order.

This is the cleanest reduced sensitivity formula for the next theorem gate.

---

## 5. What this changes

Stage 091 already showed that the minimal isotropic `3/4 + 1/4` module is forced if the geometry lane is static.

Stage 092 sharpens the remaining gap:

> the real moving-throat derivation now only needs to answer whether the geometry lane is static through `O(omega^4)` on the natural isotropic branch, or else compute the two contamination numbers `(eps_2, eps_4)`.

If both vanish, the Stage-091 result is exact.
If they do not, the contact/pole fractions are still fixed — but by the obstruction formula above rather than by the simple `3/4 + 1/4` split.
