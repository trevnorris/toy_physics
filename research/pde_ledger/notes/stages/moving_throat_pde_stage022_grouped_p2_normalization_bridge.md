# Moving-Throat PDE — Stage 022: Grouped Real `P2` Bundle and the Normalized Outgoing-Quadrupole Bridge

## Purpose

Stage 021 gave the matched one-lane bridge from

- the conservative wall/BdG/Maxwell/mixed block,
- to a passive/outgoing `l=2` fingerprint,
- with a positive static transfer factor.

That was enough to show **where** the 2.5PN quadrupole route can live, but not yet enough to speak the same language as the grouped real `P2` conservative packages or the 2.5PN normalization stack.

The next missing bridge is therefore not “solve the whole PDE.”
It is much sharper:

> lift the one-lane Stage-021 result to the grouped real `20/21/22` bundle, convert the Stage-003 and Stages-004--021 conservative operator moments into the normalized grouped-response moments used by the 2.5PN package, and isolate the exact normalization product that still has to hit the universal quadrupole target.

That is the whole point of this stage.

---

## 1. Grouped real `P2` lane setup

Work lane by lane with

`A in {20,21,22}`.

After the conservative Stage-003 + Stages-004--021 reductions, each grouped lane is described by a conservative wall/worldtube operator

`D_A^(cons)(omega) = D_{A0} + D_{A2} omega^2 + D_{A4} omega^4 + O(omega^6)`.

This notation is intentionally broad.
All of the following are already absorbed into the coefficients `D_{A0}, D_{A2}, D_{A4}`:

- the geometry-only finite-throat wall backbone,
- the stable BdG support self-energy,
- the conservative localized-Maxwell / mixed-sector self-energy.

So Stage 022 does **not** reopen the earlier reductions. It takes them as the conservative front end.

The outgoing mixed-sector dressing from Stage 021 is encoded in a transfer factor

`N_A(omega) = N_{A0} + N_{A2} omega^2 + N_{A4} omega^4 + O(omega^6)`.

The key new point is that the grouped real `P2` bundle is now described by **two** low-frequency objects per lane:

1. the conservative operator `D_A^(cons)`,
2. the outgoing transfer factor `N_A`.

---

## 2. Exact bridge from conservative operator moments to grouped-response moments

Stage 003 normalized the conservative operator itself.
The 2.5PN grouped-quadrupole package, however, is phrased most naturally in terms of the normalized **response**

`Y_A(omega) = D_{A0} / D_A^(cons)(omega)`.

Expanding,

`Y_A(omega) = 1 + u_2^(A) omega^2 + u_4^(A) omega^4 + O(omega^6)`.

The exact conversion formulas are

`u_2^(A) = - D_{A2} / D_{A0}`,

`u_4^(A) = ( D_{A2}^2 - D_{A0} D_{A4} ) / D_{A0}^2`.

This is the first important bookkeeping bridge of the stage.
It tells us exactly how the conservative wall/BdG/Maxwell data from Stages 3–4 enter the grouped real `P2` response moments used by the 2.5PN summary.

So from this point onward:

- `D_{A0}, D_{A2}, D_{A4}` are the conservative operator moments,
- `u_2^(A), u_4^(A)` are the normalized grouped-response moments.

---

## 3. Exact grouped inverse map and isotropy defects

For any grouped triple `(x_20, x_21, x_22)`, define

`xbar = (x_20 + 2 x_21 + 2 x_22) / 5`,

`a_x = (2 x_20 - x_21 - x_22) / 10`,

`b_x = (x_21 - x_22) / 2`.

Then the exact inverse map is

`x_20 = xbar + 4 a_x`,

`x_21 = xbar - a_x + b_x`,

`x_22 = xbar - a_x - b_x`.

Applied to the grouped response moments, this gives

`ubar_2 = (u_2^(20) + 2 u_2^(21) + 2 u_2^(22)) / 5`,

`a_2 = (2 u_2^(20) - u_2^(21) - u_2^(22)) / 10`,

`b_2 = (u_2^(21) - u_2^(22)) / 2`,

and similarly for `(ubar_4, a_4, b_4)`.

The corresponding anisotropy norm is

`A_2^2 = 4 a_2^2 + (4/5) b_2^2`,

and similarly for `A_4^2`.

So the grouped real `P2` isotropy gate is once again completely sharp:

`a_2 = b_2 = a_4 = b_4 = 0`.

If the three grouped lanes share the same conservative moments, then the anisotropy defects vanish automatically.

---

## 4. The outgoing bridge is controlled by an exact internal prefactor

To first order in the outgoing branch, the response correction is controlled by the internal prefactor

`Pref_A(omega) = D_{A0} N_A(omega) / D_A^(cons)(omega)^2`.

Expand

`Pref_A(omega) = P_{A0} + P_{A2} omega^2 + P_{A4} omega^4 + O(omega^6)`.

Then the exact coefficients are

`P_{A0} = N_{A0} / D_{A0}`,

`P_{A2} = ( D_{A0} N_{A2} - 2 D_{A2} N_{A0} ) / D_{A0}^2`,

`P_{A4}`
`= [ D_{A0}^2 N_{A4} - 2 D_{A0}( D_{A2} N_{A2} + D_{A4} N_{A0} ) + 3 D_{A2}^2 N_{A0} ] / D_{A0}^3`.

This is the central algebraic result of the stage.

It cleanly separates the job of the conservative moving-throat block into two pieces:

- `D_A^(cons)` tells us the grouped response moments `u_2^(A), u_4^(A)`,
- `Pref_A` tells us how strongly the outgoing `l=2` branch is transferred into each grouped lane.

---

## 5. Multiplying by the compact outgoing `l=2` fingerprint

The compact outgoing `l=2` branch already has the normalized fingerprint

`Yhat_2^(out)(omega)`
`= 1 + A omega^2 + B omega^4 + i G5 omega^5 + O(omega^6)`

with

`A = a^2 / (9 c_s^2)`,

`B = 4 a^4 / (81 c_s^4)`,

`G5 = a^5 / (27 c_s^5)`.

So the outgoing contribution in grouped lane `A` is

`delta Y_A^(out)(omega) = Pref_A(omega) * Yhat_2^(out)(omega)`.

Write the branch coefficients as

`delta Y_A^(out)(omega)`
`= K_{A0} + K_{A2} omega^2 + K_{A4} omega^4 + i Gamma_5^(A) omega^5 + ...`.

Then

`K_{A0} = P_{A0}`,

`K_{A2} = P_{A2} + A P_{A0}`,

`K_{A4} = P_{A4} + A P_{A2} + B P_{A0}`,

`Gamma_5^(A) = G5 P_{A0}`.

This formula is extremely informative.
It shows that only the **static prefactor** `P_{A0}` enters the first odd `i omega^5` coefficient.
The higher prefactor moments `P_{A2}, P_{A4}` matter for the even branch bookkeeping, but not for the leading 2.5PN odd coefficient.

So the 2.5PN normalization problem is narrower than the full grouped conservative problem.
It only needs the correct isotropic value of `P_{A0}` on the true moving-throat branch.

---

## 6. Stage-021 one-lane prototype gives explicit `N0/N2/N4`

The Stage-021 one-lane Maxwell/mixed model has

`N(omega) = (P0_proto - g_W omega^2)^2 / (Delta0 - S2 omega^2 + omega^4)^2`,

where

`Delta0 = Omega_A^2 Omega_W^2 - R^2`,

`S2 = Omega_A^2 + Omega_W^2`,

`P0_proto = Omega_A^2 g_W + R g_A`.

Its exact low-frequency coefficients are

`N0 = P0_proto^2 / Delta0^2`,

`N2 = 2 P0_proto ( P0_proto S2 - Delta0 g_W ) / Delta0^3`,

`N4`
`= [ Delta0^2 g_W^2 - 2 Delta0 P0_proto^2 - 4 Delta0 P0_proto S2 g_W + 3 P0_proto^2 S2^2 ] / Delta0^4`.

So the one-lane Maxwell/mixed prototype already gives explicit microscopic data for the first three outgoing-transfer moments.

These are not yet the full coupled moving-throat bundle values, but they are the exact reduced prototype formulas that the true grouped-lane computation must generalize.

---

## 7. Isotropic grouped-lane collapse

If the full grouped real `P2` bundle is isotropic at the conservative level, then

`D_{20,n} = D_{21,n} = D_{22,n}`,

and if the outgoing transfer is isotropic as well, then

`N_{20,n} = N_{21,n} = N_{22,n}`

for each low-frequency coefficient we keep.

In that case,

`u_2^(20) = u_2^(21) = u_2^(22) = ubar_2`,

`u_4^(20) = u_4^(21) = u_4^(22) = ubar_4`,

`P_{20,n} = P_{21,n} = P_{22,n} = P_n`,

`K_{20,n} = K_{21,n} = K_{22,n} = K_n`,

and therefore all grouped anisotropy defects vanish:

`a_2 = b_2 = a_4 = b_4 = 0`.

So Stage 022 makes the 3PN and 2.5PN grouped stories meet cleanly:

- the 3PN grouped real `P2` isotropy gate is a statement about the conservative operator moments,
- the 2.5PN grouped quadrupole normalization is a statement about the isotropic static outgoing prefactor `P_0`.

---

## 8. The exact normalization product still to be hit

The universal GR quadrupole target is

`gamma_GR = 2 G / (5 c^5)`.

From the Stage-022 branch formula,

`Gamma_5 = P_0 * a^5 / (27 c_s^5)`.

Including the source-map factor `mhat_0`, the invariant normalization condition is therefore

`mhat_0^2 * P_0 * a^5 / (27 c_s^5) = 2 G / (5 c^5)`.

Equivalently,

`mhat_0^2 * P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

This is the sharpest expression of the remaining theorem gap I can give at this stage.

Everything that came before was needed to make this quantity well defined.
Now the problem is no longer diffuse:

> compute the isotropic moving-throat value of the invariant product `mhat_0^2 P_0`.

On the natural source-map branch `mhat_0 = 1 + O(a^2/r^2)`, this reduces to the familiar 2.5PN target

`P_0 = 54 G c_s^5 / (5 a^5 c^5)`

at leading point-particle order.

If, in addition, the outgoing branch is carried by a **constant isotropic prefactor** (`P_2 = P_4 = 0`), then the even branch coefficients automatically become

`K_2 = P_0 a^2 / (9 c_s^2)`,

`K_4 = 4 P_0 a^4 / (81 c_s^4)`.

That is exactly the constant-prefactor branch isolated in the 2.5PN normalization package.

---

## 9. What this stage accomplishes physically

This stage does three things that were missing before.

### 9.1 It unifies the Stage-003 and Stages-004--021 conservative data with the 2.5PN grouped language

Before this note, the conservative moving-throat derivation was naturally phrased in terms of operator moments.
The 2.5PN grouped summary, however, is phrased in terms of normalized grouped-response coefficients and branch-normalization data.

Stage 022 now gives the exact conversion formulas between those languages.

### 9.2 It identifies the true internal quantity that controls 2.5PN

The leading odd coefficient does **not** depend on every conservative detail equally.
At the universal point-particle 2.5PN level, the only internal quantity that matters is the isotropic static outgoing prefactor `P_0 = N_0 / D_0` (or `mhat_0^2 P_0` in invariant form).

That is a much sharper theorem target than “solve the whole quadrupole branch.”

### 9.3 It explains exactly what the full PDE still has to provide

The full moving-throat PDE still has to do two hard things:

1. produce the true grouped-lane conservative coefficients `D_{A0}, D_{A2}, D_{A4}` and outgoing transfer coefficients `N_{A0}, N_{A2}, N_{A4}` on the coupled isotropic branch,
2. and land the invariant product `mhat_0^2 P_0` on the target value above.

That is now a concrete job description for the missing PDE.

---

## 10. Immediate next step

The right next derivation is now narrower than before.

Do **not** jump yet to the full nonlinear theorem.
Instead:

1. keep the grouped-lane formulas from this stage,
2. compute the actual coupled grouped real `20/21/22` conservative coefficients for the wall/BdG/Maxwell/mixed system,
3. evaluate the corresponding `P_0, P_2, P_4` prefactor data on the isotropic branch,
4. and test whether the invariant normalization product `mhat_0^2 P_0` lands on the universal target.

That is the smallest next theorem gate that directly attacks the remaining 2.5PN/4PN bottleneck.
