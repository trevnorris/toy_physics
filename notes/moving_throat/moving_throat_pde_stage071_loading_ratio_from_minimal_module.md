# Moving-Throat PDE — Stage 71: Loading-Ratio Extraction from the Minimal Isotropic Quadrupole Module

## Purpose

Stage 70 left the explicit Family-1 support/source theorem in its cleanest reduced form:

`rho_alpha = alpha_req / alpha_mix`

is the only quantity still needed from the actual passive/outgoing quadrupole branch.

The 2.5PN quadrupole audit had already isolated the smallest viable isotropic conservative precursor as

`Y_Q^cons(omega) = c0 + c1 / (1 - omega^2 / Omega_Q^2),`

with the unique positive-real coefficients

`c0 = 3/4,`
`c1 = 1/4,`
`Omega_Q = 3 c_s / (2 a).`

The next honest calculation is therefore to identify how that conservative precursor maps onto the support/source loading ratio.

This stage shows that under the natural **contact-plus-pole** reading of the explicit support/source branch,

- the mixed baseline contributes the static contact fraction,
- the extra support lane contributes the finite conservative pole,
- and the loading ratio is therefore fixed exactly by the precursor coefficients.

The main result is

`rho_alpha = 4/3,`
`zeta_req = 1/3,`
`Pi_tr = (4/3) C_mix.`

So the minimal isotropic quadrupole precursor lands in the exact **symmetric-lowest-twin** support regime and does not require non-twin asymmetry.

---

## 1. Natural contact-plus-pole identification

On the explicit support/source branch, the mixed baseline loading `alpha_mix` is the conservative directional loading already present before any extra support pole is added.
If the final passive/outgoing quadrupole branch is represented by one extra conservative support pole, then the natural normalized conservative precursor is

`Y_Q^cons(omega)`
`= alpha_mix / alpha_req`
`  + (alpha_req - alpha_mix)/alpha_req * 1/(1 - omega^2/Omega_Q^2).`

Equivalently, in the pure loading-ratio variable

`rho_alpha := alpha_req / alpha_mix,`

the same precursor is

`Y_Q^cons(omega)`
`= 1/rho_alpha`
`  + (rho_alpha - 1)/rho_alpha * 1/(1 - omega^2/Omega_Q^2).`

So the contact fraction and pole residue are

`c0 = 1/rho_alpha,`
`c1 = (rho_alpha - 1)/rho_alpha,`

with

`c0 + c1 = 1`

as required by the normalized static limit.

This gives the exact inverse formulas

`rho_alpha = 1/c0 = 1/(1-c1),`

`zeta_req = rho_alpha - 1 = c1/c0.`

So the support/source loading ratio is directly encoded in the static-contact / pole-residue split of the conservative quadrupole precursor.

---

## 2. Matching to the minimal isotropic quadrupole module

The 2.5PN quadrupole audit already fixed the smallest viable isotropic conservative precursor to

`c0 = 3/4,`
`c1 = 1/4.`

Inserting these into the exact inverse formulas above gives immediately

`rho_alpha = 1 / (3/4) = 4/3,`

`zeta_req = c1/c0 = (1/4)/(3/4) = 1/3.`

So the natural contact-plus-pole interpretation of the minimal isotropic quadrupole branch fixes the explicit support demand exactly:

`alpha_req / alpha_mix = 4/3,`
`(alpha_req - alpha_mix)/alpha_mix = 1/3.`

In product language, since Stage 68 proved

`Pi_tr / C_mix = alpha_req / alpha_mix,`

the same result is

`Pi_tr = (4/3) C_mix.`

---

## 3. Regime classification

Stage 35 split the support regimes into:

- `Pi_tr <= C_mix`              : mixed-only already enough,
- `C_mix < Pi_tr <= 2 C_mix`    : symmetric lowest twin enough,
- `Pi_tr > 2 C_mix`             : non-twin asymmetry required.

Because the minimal isotropic branch gives

`Pi_tr = (4/3) C_mix,`

it lands exactly in the middle regime:

`C_mix < Pi_tr < 2 C_mix.`

So the minimal isotropic passive/outgoing quadrupole branch:

- does require extra support beyond the mixed baseline,
- but requires only the symmetric lowest-twin lane,
- and does **not** require non-twin asymmetry.

Equivalently, in support-ratio language,

`zeta_req = 1/3 < 1,`

so the exact symmetric lowest twin already suffices.

---

## 4. What this changes

Before this step, the remaining reduced question was still “what value of `rho_alpha` does the actual passive/outgoing branch select?”

After this step, the minimal isotropic quadrupole precursor provides a very concrete answer:

> if the actual passive/outgoing branch is the natural contact-plus-pole realization of the minimal isotropic conservative module, then it selects
>
> `rho_alpha = 4/3`.

That is a much stronger statement than the earlier support-side ceiling alone.

It says the explicit Family-1 branch is not merely *capable* of surviving the outgoing quadrupole demand.
On the natural minimal isotropic branch, it only has to accommodate a support ratio of one third above the mixed baseline.

---

## 5. Best current theorem statement after Stage 71

Under the natural unblocked contact-plus-pole identification of the passive/outgoing isotropic quadrupole branch,

`Y_Q^cons(omega) = c0 + c1/(1 - omega^2/Omega_Q^2),`

with

`c0 = 3/4,`
`c1 = 1/4,`

the explicit Family-1 support theorem is driven by the exact loading ratio

`rho_alpha = 4/3,`

equivalently

`zeta_req = 1/3,`
`Pi_tr = (4/3) C_mix.`

So the reduced moving-throat PDE has advanced from a vague outgoing-branch loading question to a sharp statement:

> **if the actual passive/outgoing quadrupole branch realizes the minimal isotropic conservative precursor in the natural contact-plus-pole way, then the explicit support/source side is comfortably compatible with it.**
