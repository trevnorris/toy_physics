# 5PN continuation notes — Stages 272–275

These stages pick up from the Stage-271 verdict that the actual isotropic grouped-`P2` / geometry branch is already conservatively clean and that the explicit Family-1 support/source side is automatic there. The point of the continuation is therefore **not** to reopen the bulk GNLS or to keep probing the scalar/geometry lane. It is to collapse the remaining reduced theorem gap all the way down to the outgoing quadrupole normalization factor.

The main new result is that the whole reduced 2.5PN / 4PN endgame can now be written as one outgoing scalar defect

a) `chi_Q - 1`,

with a clean canonical branch `chi_Q = 1`, and with an exact three-parameter isotropic DtN deformation algebra that shows which deformations can actually move it.

---

## Stage 272 — actual isotropic outgoing reduction

Start from the actual isotropic grouped-`P2` one-pole conservative carrier

`Y_Q^cons(omega) = N_Q [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2) ]`.

Expanding at low frequency gives

`K0 = N_Q`,

`K2 = N_Q / (4 Omega_Q^2)`,

`K4 = N_Q / (4 Omega_Q^4)`.

So once the actual isotropic one-pole branch is accepted, every conservative low-frequency coefficient is carried by one scalar defect `N_Q`.

Adding the retarded outgoing normalization factor gives

`Y_Q^ret(omega) = N_Q [ 3/4 + (1/4)/(1 - omega^2/Omega_Q^2 - i chi_Q sigma_Q^can omega^5) ] + O(omega^6)`.

The odd coefficient is therefore

`Gamma5 = N_Q chi_Q sigma_Q^can / 4`.

With the canonical target normalization

`K0_target = 54 G c_s^5 / (5 a^5 c^5)`,

`sigma_Q^can = 4 a^5 / (27 c_s^5)`,

one gets the exact bridge

`K0_target sigma_Q^can / 4 = 2 G / (5 c^5)`.

So the full odd normalization condition is

`mhat_0^2 chi_Q N_Q = 1`.

On the natural source-map branch `mhat_0 -> 1`, the remaining reduced obstruction is purely outgoing:

`N_Q = 1 / chi_Q`.

So after Stage 272 the live reduced theorem gap is no longer a mixed support/source/outgoing problem. It is one outgoing normalization number.

---

## Stage 273 — canonical outgoing `l=2` DtN match

Using the exact compact outgoing spherical `l=2` DtN fingerprint

`Lambda_2^out(z) = -3 + z^2/3 + z^4/9 + i z^5/9 - 2 z^6/27 - i z^7/27 + O(z^8)`,

and the exact normalization relation

`Yhat_2^out(z) = -3 / Lambda_2^out(z)`,

one gets

`Yhat_2^out(z) = 1 + z^2/9 + 4 z^4/81 + i z^5/27 - 11 z^6/729 - i z^7/243 + O(z^8)`.

Now match this to the grouped-`P2` one-pole-plus-contact ansatz in the dimensionless variable `z = omega a / c_s`:

`Yhat_Q^ret(z) = 3/4 + (1/4)/(1 - alpha z^2 - i chi_Q B z^5) + O(z^6)`.

Exact coefficient matching gives

`alpha = 4/9`,

`B chi_Q = 4/27`.

Since

`alpha = c_s^2 / (a^2 Omega_Q^2)`,

the canonical pole scale is fixed to

`Omega_Q = 3 c_s / (2 a)`.

Then

`sigma_Q^can = 9 / (8 Omega_Q^5) = 4 a^5 / (27 c_s^5)`,

so the dimensionless odd coefficient is exactly `B = 4/27`, and therefore

`chi_Q = 1`.

This is the cleanest current proof that the canonical compact passive/outgoing `l=2` branch is the canonical reduced GR branch.

---

## Stage 274 — isotropic DtN deformation algebra

Now deform the isotropic outgoing DtN branch by

`Lambda_2^def(z)`
`= S Lambda_2^out(beta z) + Sigma0 + Sigma2 z^2 + Sigma4 z^4 + i Sigma5 z^5 + O(z^6)`.

Write again

`Yhat_2^def(z) = -3 / Lambda_2^def(z)`.

If the **normalized even moments** are required to stay canonical, i.e.

`y2 / y0 = 1/9`,

`y4 / y0 = 4/81`,

then the even deformation coefficients are fixed exactly:

`Sigma2 = - S beta^2 / 3 + S / 3 - Sigma0 / 9`,

`Sigma4 = - S beta^4 / 9 + S / 9 - Sigma0 / 27`.

After that elimination, the conservative and outgoing scalar defects are

`N_Q = y0 = 3 / (3 S - Sigma0)`,

`chi_Q = 27 y5 / y0 = 3 (S beta^5 + 9 Sigma5) / (3 S - Sigma0)`.

So once the even moments are pinned, only three isotropic DtN deformations can move the outgoing theorem:

- the scale deformation `beta`,
- the static DtN shift `Sigma0`,
- the odd `z^5` DtN shift `Sigma5`.

Everything else is already absorbed into maintaining the canonical even branch.

---

## Stage 275 — outgoing-defect linearization and exact no-shift family

The exact no-shift condition `chi_Q = 1` can be solved directly:

`Sigma0 = 3 S (1 - beta^5) - 27 Sigma5`.

So the canonical outgoing branch is **not unique**. Once the even moments are preserved, there is an exact two-parameter isotropic deformation family that still leaves `chi_Q = 1`.

Linearizing around the canonical compact branch,

`S = 1 + delta_S`,

`beta = 1 + delta_beta`,

`Sigma0 = delta_Sigma0`,

`Sigma5 = delta_Sigma5`,

one gets the first-order outgoing defect

`Delta_Q := chi_Q - 1`

`= 5 delta_beta + delta_Sigma0 / 3 + 9 delta_Sigma5 + O(2)`.

The overall amplitude deformation `delta_S` cancels out at first order. So at linear order the outgoing theorem is controlled only by:

1. the pole rescaling `beta`,
2. the static DtN shift `Sigma0`,
3. the odd DtN shift `Sigma5`.

A useful special slice is the already conservative-normalized branch `N_Q = 1`, which forces

`Sigma0 = 3 (S - 1)`.

On that slice the remaining outgoing defect reduces to

`chi_Q = S beta^5 + 9 Sigma5`.

So once the conservative isotropic branch is clean, the remaining reduced obstruction is entirely retarded.

---

## Best current reading after Stage 275

The reduced theorem gap is now smaller than it was at Stage 271.

1. The actual isotropic grouped-`P2` conservative branch carries one scalar defect `N_Q`.
2. On the natural source map, the remaining reduced odd obstruction is purely outgoing and is measured by `chi_Q`.
3. The canonical compact passive/outgoing `l=2` DtN branch gives `chi_Q = 1` exactly.
4. The most general isotropic DtN deformation that preserves the canonical even moments moves `chi_Q` only through `beta`, `Sigma0`, and `Sigma5`.
5. Near the canonical branch, the exact first-order outgoing defect is

   `Delta_Q = 5 delta_beta + delta_Sigma0 / 3 + 9 delta_Sigma5 + O(2)`.

So the next honest theorem gate is no longer on the scalar/geometry or explicit support/source side. It is to determine whether the actual passive/outgoing moving-throat DtN branch is the canonical one (`chi_Q = 1`) or, if not, which concrete isotropic DtN deformation data `beta`, `Sigma0`, `Sigma5` shift it away.
