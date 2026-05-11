# Moving-Throat PDE — Stage 095: First Nonzero Geometry Contamination Appears Only at Second Order in Anisotropy/Mixing

## Purpose

Stage 77 proved that on the natural isotropic branch the geometry lane is dynamically inert through `O(omega^4)` with respect to the grouped real `P2` conservative quadrupole module:

`eps_2 = eps_4 = 0`.

The next useful question is then:

> if isotropy is weakly broken and the scalar/geometry lane mixes with the grouped-`P2` quadrupole carrier, how do the contamination numbers turn on?

This stage answers that in the smallest exact reduced model.

The main result is:

- the first nonzero geometry contamination is **quadratic** in the mixing parameter,
- so the isotropic result is stable,
- and any deviation from the Stage-74 `3/4 + 1/4` split begins only at `O(chi^2)`.

---

## 1. Minimal mixed scalar-geometry / grouped-`P2` model

Take one grouped-`P2` quadrupole carrier `q` and one scalar/geometry mode `g`.

Write the conservative reduced kernel as

`D_q(omega) = K_stat + K_pole /(1 - omega^2/Omega_Q^2)`,

`D_g(omega) = G_0 + G_2 omega^2 + G_4 omega^4 + O(omega^6)`,

and let weak anisotropy or weak operator non-commutation generate a bilinear mixing term

`chi M_0 q g`.

Then the quadratic action is

`L = (1/2) q D_q q + (1/2) g D_g g + chi M_0 q g + J q`.

Integrating out the scalar/geometry mode gives the exact effective quadrupole kernel

`D_eff(omega) = D_q(omega) - chi^2 M_0^2 / D_g(omega)`.

So the whole contamination problem is encoded in the low-frequency expansion of the Schur-complement term.

---

## 2. Exact low-frequency contamination coefficients

Expand

`1 / D_g(omega)`

through `O(omega^4)`:

`1 / D_g(omega)`
`= 1/G_0`
`  - (G_2/G_0^2) omega^2`
`  + (G_2^2/G_0^3 - G_4/G_0^2) omega^4`
`  + O(omega^6)`.

So

`- chi^2 M_0^2 / D_g(omega)`
`= - chi^2 M_0^2 / G_0`
`  + chi^2 M_0^2 G_2 / G_0^2 * omega^2`
`  + chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3 * omega^4`
`  + O(omega^6)`.

The first term is only a static renormalization and can be absorbed into the geometry contact slot.
The dynamic contaminations are therefore

`K_(g,2)^eff = chi^2 M_0^2 G_2 / G_0^2`,

`K_(g,4)^eff = chi^2 M_0^2 (G_0 G_4 - G_2^2) / G_0^3`.

So the obstruction numbers of Stage 75 are

`eps_2 = Omega_Q^2 K_(g,2)^eff / K_pole`,

`eps_4 = Omega_Q^4 K_(g,4)^eff / K_pole`,

and both satisfy

`eps_2 = O(chi^2)`,

`eps_4 = O(chi^2)`.

That is the exact first stability theorem beyond Stage 77.

---

## 3. Consequence for the contact/pole fractions

Insert the small contamination into the Stage-75 obstruction formula

`c_pole = (1 + eps_4) / [ 4 (1 + eps_2)^2 ]`.

Expanding for small `chi` gives

`c_pole = 1/4 [ 1 + eps_4 - 2 eps_2 + O(chi^4) ]`.

Because both `eps_2` and `eps_4` are already `O(chi^2)`, the deviation from the `1/4` pole fraction is itself `O(chi^2)`.

So the Stage-74 `3/4 + 1/4` module is not just exact on the isotropic branch — it is also **perturbatively stable** against weak anisotropy/mixing.

---

## 4. Best current interpretation

Stages 77–78 together now give the exact reduced answer to the geometry-lane check:

1. on the natural isotropic branch,
   
   `eps_2 = eps_4 = 0`;

2. the first nonzero geometry contamination requires an explicit `l=0 <-> l=2` mixing source,
   and even then it appears only at
   
   `O(chi^2)`.

So the current reduced hierarchy supports the clean conservative reading:

- the grouped real `P2` branch is the dynamic quadrupole carrier,
- the geometry lane is static at the relevant order on the natural branch,
- and the Stage-74 `3/4 + 1/4` split is the correct actual branch value unless a genuine symmetry-breaking or extra `l=2` geometry pole is later found.

---

## 5. Best current theorem statement after Stage 78

Inside the present reduced hierarchy,

`eps_2 = eps_4 = 0`

exactly on the isotropic branch, and for weak symmetry breaking

`eps_2, eps_4 = O(chi^2)`.

So the dynamic geometry obstruction is absent on the actual isotropic branch and only enters at second order once an explicit mixing mechanism is turned on.
