# 5PN continuation notes — stages 18 through 20

These stages do two things in parallel.

First, they sharpen the **isotropic full-bundle target surface** that the actual
moving-throat PDE must hit if the calibrated 5PN / 2.5PN / 4PN branch is real.
Second, they reinstate the omitted conservative Maxwell/mixed `Z_2,Z_4` sector in
the weak-axisymmetric even-gate problem and show exactly how it removes the fake
freedom left in the Stage-17 lower-bound picture.

## Stage 18 — exact isotropic full-bundle target surface

Work on the isotropic grouped-lane branch with

a) conservative operator moments

a) `D0 = K - B0 - Z0`,

b) `D2 = -(M + B2 + Z2)`,

c) `D4 = -(B4 + Z4)`.

Then the normalized grouped-response and outgoing-prefactor moments are

`u2 = -D2 / D0`,

`u4 = (D2^2 - D0 D4) / D0^2`,

`P0 = N0 / D0`,

`P2 = (D0 N2 - 2 D2 N0) / D0^2`,

`P4 = (D0^2 N4 - 2 D0 (D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3`.

The exact isotropic 5PN conservative one-pole defect is

`u4 - 4 u2^2 = [ D0 (B4 + Z4) - 3 (M + B2 + Z2)^2 ] / D0^2`.

So the isotropic one-pole surface is exactly

`D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`.

The constant-prefactor outgoing branch is also exact:

`P2 = 0  ->  N2 = 2 D2 N0 / D0`,

`P4 = 0  ->  N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

Finally, the universal 2.5PN / 4PN normalization target is

`mhat0^2 P0 = 54 G c_s^5 / (5 a^5 c^5)`

or equivalently

`mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`.

So the full isotropic moving-throat bundle must land on one exact combined target
surface:

1. `D0 = K - B0 - Z0`,
2. `D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`,
3. `mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`,
4. `N2 = 2 D2 N0 / D0`,
5. `N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

That is the sharpest isotropic full-bundle statement we have so far.

## Stage 19 — exact `Z`-sector bridge back into the even gates

The missing conservative Maxwell/mixed moments are

`Z0 = Q / Delta`,

`Z2 = (Q S2 - H Delta) / Delta^2`,

`Z4 = [ Q (S2^2 - Delta) - S2 H Delta ] / Delta^3`.

Their exact first variations are

`dZ0 = (Delta dQ - Q dDelta) / Delta^2`,

`dZ2 = [ Delta (-Delta dH - H dDelta + Q dS2 + S2 dQ) + 2 dDelta (Delta H - Q S2) ] / Delta^3`,

`dZ4 = -[ Delta^2 H dS2 + Delta^2 S2 dH + Delta^2 dQ - 2 Delta H S2 dDelta`
`         - 2 Delta Q S2 dS2 - 2 Delta Q dDelta - Delta S2^2 dQ + 3 Q S2^2 dDelta ] / Delta^4`.

Therefore the conservative `Z`-sector contributions to the two live 5PN even gates are

`K1_Z = -dZ2 - dZ0/9`,

`H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27`.

After inserting the exact Stage-13 normalized similarity kernel, both become honest
functions of the mixed-sector similarity directions `alpha_W` and `alpha_R`.

On the positive constructive slice

- `G_U = 5`, `G_W = 7`, `R = 2`,
- `Omega_U = 3`, `Omega_W = 4`,
- `chi_0 = 3/2`, `delta_U = 2/3`,
- `E_* = 1/4`, `F_* = 5/6`,

we get

`K1_Z = (78623501/25004700) alpha_OmegaU + (733046/6251175) alpha_R`
`       - (59010631/25004700) alpha_U - (32134513/50009400) alpha_W`,

`H_even,Z = -(28906377971/21003948000) alpha_OmegaU`
`           - (1174937411/21003948000) alpha_R`
`           + (11102468471/10501974000) alpha_U`
`           + (5617869293/21003948000) alpha_W`.

In particular, the pure directions give

`alpha_W:  K1_Z = -32134513/50009400,   H_even,Z = 5617869293/21003948000`,

`alpha_R:  K1_Z = 733046/6251175,       H_even,Z = -1174937411/21003948000`.

So the omitted `Z_2,Z_4` block does exactly what Stage 17 said it still had to do:
it activates the previously untouched mixed directions.

## Stage 20 — full even-gate solve on the constructive branch

Once wall-only, pure-BdG, and the reinstated `Z`-sector are combined, the full
constructive-slice even-gate matrix is

`Gate_full(slice) =`

`[ -25/9,  -32134513/50009400,   91017569/25004700,    733046/6251175,`
`  -71404699/25004700,  -8/9,  4/3 ]`

`[  52/27,  5617869293/21003948000, -30905427529/10501974000,`
`  -1174937411/21003948000, 55109414029/21003948000, 32/81, -16/27 ]`.

The matrix still has rank `2`, so the full even-gate intersection is still
5-dimensional. That part is unsurprising: there are still only two even equations.

The new structural fact is the mixed-sector minor:

`det Gate_(alpha_W, alpha_R) = 942737330573 / 205838690400000 != 0`.

So on this branch the full even system can solve directly for the previously
untouched mixed directions:

`alpha_W =`
`  (14503089433000/942737330573) alpha_K`
`+ (30450672110098/942737330573) alpha_OmegaU`
`- (29120459867142/942737330573) alpha_U`
`- (18876066395200/25453907925471) beta_B`
`+ (9438033197600/8484635975157) beta_varpi`,

`alpha_R =`
`  (101802968743000/942737330573) alpha_K`
`+ (189815725996721/942737330573) alpha_OmegaU`
`- (188832473718440/942737330573) alpha_U`
`+ (89510801038400/25453907925471) beta_B`
`- (44755400519200/8484635975157) beta_varpi`.

So there are no longer pure `alpha_W` or pure `alpha_R` null directions on the full
constructive branch. The mixed-sector freedom was never genuinely unconstrained; it
was just hidden in the omitted `Z_2,Z_4` block.

## Net result after Stage 20

The continuation point is sharper again.

1. The isotropic full-bundle target surface is now exact and explicit.
2. The omitted conservative `Z_2,Z_4` sector has been reinstated exactly.
3. On a clean constructive branch it does the precise job Stage 17 predicted:
   it removes the fake freedom in the mixed directions `alpha_W, alpha_R`.
4. So the remaining work is no longer “some mixed-sector freedom somewhere.”
   It is to compute the actual overlap data on the true moving-throat branch and see
   whether they land on the isotropic full-bundle target surface from Stage 18.
