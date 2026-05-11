# Moving-Throat PDE — Stage 099: The Reduced Finish Line After the Geometry-Lane Check

## Purpose

This note records the sharpest reduced-program verdict reached so far.

After the grouped-`P2` conservative closure, the geometry-lane check, the explicit Family-1 support/source analysis, and the passive/outgoing quadrupole reduction, the remaining theorem gap is no longer diffuse.

It is one scalar normalization datum.

---

## 1. What is now fixed inside the reduced hierarchy

On the actual natural isotropic branch we now have:

1. `eps_2 = eps_4 = 0`, so the geometry lane is dynamically inert through `O(omega^4)` with respect to the grouped real quadrupole carrier.
2. The conservative quadrupole module is exactly

   `Yhat_Q^cons(omega) = 3/4 + (1/4)/(1 - omega^2/Omega_Q^2)`.

3. Therefore

   `rho_alpha = 4/3`,

   `zeta_req = 1/3` in the unblocked limit.

4. On the explicit Family-1 branch the support/source side is automatic even in the admissible blocked regime.

So the reduced program is no longer carrying separate uncertainties in

- geometry contamination,
- contact/pole splitting,
- support/source sufficiency,
- or explicit Family-1 viability.

---

## 2. Single remaining reduced theorem gate

Write the actual isotropic passive/outgoing branch in canonical invariant form as

`Kbar_Q^cons(omega) = Kbar_0 [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]`.

Then define

`N_Q := Kbar_0 / Kbar_0^target,`

with

`Kbar_0^target = 64 G Omega_Q^5 / (45 c^5)`

or equivalently

`Kbar_0^target = 54 G c_s^5 / (5 a^5 c^5)`

after using `Omega_Q = 3 c_s/(2a)`.

The full reduced GR-like point-particle 2.5PN closure on the actual isotropic branch is now equivalent to

`N_Q = 1`.

Everything else has already been reduced away.

---

## 3. Practical meaning

The completed moving-throat PDE no longer needs to answer a large family of loosely related questions in order to close the reduced theorem.

It only needs to determine the actual radiative normalization of the passive/outgoing grouped-`P2` quadrupole branch.

In other words, the open PDE task is now exactly one of the following equivalent statements:

- compute `Kbar_0` on the actual passive/outgoing branch,
- compute `Gammabar_5` on that branch,
- or compute the scalar defect `N_Q - 1`.

The three questions are equivalent on the actual isotropic branch.

---

## 4. Best current theorem statement

Inside the present reduced hierarchy, the moving-throat PDE program has reached the following finish line:

> the actual isotropic grouped-`P2` one-pole branch is conservatively clean, the explicit Family-1 support/source side is already sufficient, and the only remaining reduced theorem gap is the single passive/outgoing quadrupole normalization defect `N_Q - 1`.

That is the narrowest and strongest honest carry-forward statement available before solving the full moving-throat PDE normalization problem.
