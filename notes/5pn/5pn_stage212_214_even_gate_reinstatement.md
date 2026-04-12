# 5PN continuation notes — stages 212 through 214

This session pushed the weak-axisymmetric program one step past the Stage-209–211 sieve.

The open question after Stage 211 was no longer whether the normalized similarity-orbit
package controls `Xi_load = P1/P0`. That part was already clean. The real next gate was the one
flagged in the notes after Stage 17:

> reinstate the conservative Maxwell/mixed `Z_2,Z_4` sector and see whether it removes the fake
> freedom left in the lower-bound even-gate picture.

These three stages do exactly that.

---

## Stage 212 — exact isotropic full-bundle target surface

Files:
- `5pn_stage212_isotropic_target_surface.py`
- `5pn_stage212_isotropic_target_surface_output.txt`

### What was frozen

Work on the isotropic grouped-lane branch with

```text
D0 = K - B0 - Z0
D2 = -(M + B2 + Z2)
D4 = -(B4 + Z4)
```

and compile the normalized grouped-response and outgoing-prefactor moments

```text
u2 = -D2 / D0
u4 = (D2^2 - D0 D4) / D0^2

P0 = N0 / D0
P2 = (D0 N2 - 2 D2 N0) / D0^2
P4 = (D0^2 N4 - 2 D0(D2 N2 + D4 N0) + 3 D2^2 N0) / D0^3.
```

### Exact one-pole defect

The script proves

```text
u4 - 4 u2^2 = [ D0 (B4 + Z4) - 3 (M + B2 + Z2)^2 ] / D0^2.
```

So the isotropic conservative one-pole surface is exactly

```text
D0 (B4 + Z4) = 3 (M + B2 + Z2)^2.
```

### Exact constant-prefactor branch

Imposing

```text
P2 = 0
P4 = 0
```

gives the exact outgoing conditions

```text
N2 = 2 D2 N0 / D0

N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2.
```

### Exact 2.5PN / 4PN normalization hit

The universal normalization condition is

```text
mhat0^2 P0 = 54 G c_s^5 / (5 a^5 c^5),
```

i.e.

```text
mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5).
```

### Most useful Stage-212 result

The completed isotropic moving-throat bundle must land on one exact algebraic target surface:

1. `D0 = K - B0 - Z0`,
2. `D0 (B4 + Z4) = 3 (M + B2 + Z2)^2`,
3. `mhat0^2 N0 / D0 = 54 G c_s^5 / (5 a^5 c^5)`,
4. `N2 = 2 D2 N0 / D0`,
5. `N4 = [2 D0 (D2 N2 + D4 N0) - 3 D2^2 N0] / D0^2`.

So the actual PDE does not need to “show 5PN” in any vague sense. It needs to hit that surface.

---

## Stage 213 — exact `Z`-sector bridge back into the live even gates

Files:
- `5pn_stage213_z_sector_even_gate_bridge.py`
- `5pn_stage213_z_sector_even_gate_bridge_output.txt`

### Conservative Maxwell/mixed moments reinstated exactly

The omitted conservative block is

```text
Z0 = Q / Delta
Z2 = (Q S2 - H Delta) / Delta^2
Z4 = [ Q (S2^2 - Delta) - S2 H Delta ] / Delta^3
```

with

```text
Delta = Omega_U^2 Omega_W^2 - R^2
Q     = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2
H     = G_U^2 + G_W^2
S2    = Omega_U^2 + Omega_W^2.
```

The script differentiates these exactly and verifies

```text
dZ0 = (Delta dQ - Q dDelta) / Delta^2

dZ2 = [ Delta(-Delta dH - H dDelta + Q dS2 + S2 dQ)
        + 2 dDelta (Delta H - Q S2) ] / Delta^3

dZ4 = -[ Delta^2 H dS2 + Delta^2 S2 dH + Delta^2 dQ
          - 2 Delta H S2 dDelta - 2 Delta Q S2 dS2
          - 2 Delta Q dDelta - Delta S2^2 dQ + 3 Q S2^2 dDelta ] / Delta^4.
```

### Exact even-gate contribution

The missing `Z`-sector contribution to the two live even gates is

```text
K1_Z     = -dZ2 - dZ0/9
H_even,Z = -dZ4 + (2/3) dZ2 + dZ0/27.
```

### Constructive-slice evaluation

On the same constructive slice recorded in the notes,

```text
G_U = 5,  G_W = 7,  R = 2,
Omega_U = 3,  Omega_W = 4,
chi_0 = 3/2,  delta_U = 2/3,
E_* = 1/4,  F_* = 5/6,
```

the script gives

```text
K1_Z = (78623501/25004700) alpha_OmegaU + (733046/6251175) alpha_R
       - (59010631/25004700) alpha_U - (32134513/50009400) alpha_W

H_even,Z = -(28906377971/21003948000) alpha_OmegaU
           - (1174937411/21003948000) alpha_R
           + (11102468471/10501974000) alpha_U
           + (5617869293/21003948000) alpha_W.
```

In particular, the pure mixed directions are now active:

```text
alpha_W :  K1_Z = -32134513/50009400,
           H_even,Z = 5617869293/21003948000

alpha_R :  K1_Z = 733046/6251175,
           H_even,Z = -1174937411/21003948000.
```

### Most useful Stage-213 result

This is the first exact executable proof that the omitted `Z_2,Z_4` block does the precise job the
Stage-17 lower-bound picture predicted: it activates the previously untouched mixed directions
`alpha_W` and `alpha_R` in the even-gate problem.

---

## Stage 214 — full constructive-slice even-gate solve

Files:
- `5pn_stage214_full_even_gate_constructive_slice.py`
- `5pn_stage214_full_even_gate_constructive_slice_output.txt`

### Lower-bound matrix before `Z_2,Z_4`

Using only the matched wall block and the explicit one-mode BdG block with

```text
K = 2,
M = 3,
B0 = 2,
varpi = 3,
```

the lower-bound gate matrix in the seven directions

```text
(alpha_K, alpha_W, alpha_U, alpha_R, alpha_OmegaU, beta_B, beta_varpi)
```

has rank `2`, but its `alpha_W` and `alpha_R` columns vanish exactly.
So the lower-bound picture leaves those mixed directions completely untouched.

### Full matrix after exact `Z`-sector reinstatement

Adding the Stage-213 `Z`-sector contributions gives the full constructive-slice even-gate matrix

```text
[ -25/9,  -32134513/50009400,   91017569/25004700,    733046/6251175,
  -71404699/25004700,  -8/9,  4/3 ]

[  52/27,  5617869293/21003948000, -30905427529/10501974000,
  -1174937411/21003948000, 55109414029/21003948000, 32/81, -16/27 ]
```

with rank still `2`, so the full even-gate intersection remains five-dimensional.

### The structural fact that matters

The mixed-sector minor is now nonzero:

```text
det Gate_(alpha_W, alpha_R) = 942737330573 / 205838690400000 != 0.
```

So the apparently free mixed directions of the lower-bound picture were never genuine null
modes. They were only hidden in the omitted `Z_2,Z_4` block.

### Exact solve for the mixed directions

The full constructive-slice even-gate system solves directly for
`alpha_W` and `alpha_R`:

```text
alpha_W =
  (14503089433000/942737330573) alpha_K
+ (30450672110098/942737330573) alpha_OmegaU
- (29120459867142/942737330573) alpha_U
- (18876066395200/25453907925471) beta_B
+ (9438033197600/8484635975157) beta_varpi

alpha_R =
  (101802968743000/942737330573) alpha_K
+ (189815725996721/942737330573) alpha_OmegaU
- (188832473718440/942737330573) alpha_U
+ (89510801038400/25453907925471) beta_B
- (44755400519200/8484635975157) beta_varpi.
```

The script verifies by exact back-substitution that these kill both live even gates.

### Most useful Stage-214 result

There are no longer pure `alpha_W` or pure `alpha_R` null directions on the full constructive
branch. Once the conservative `Z_2,Z_4` sector is restored, the mixed-sector freedom is no longer
an open qualitative question. It is algebraically slaved to the remaining five directions.

---

## Best current status after stages 212–214

The next theorem gate is sharper again.

1. The exact isotropic full-bundle target surface is now explicit.
2. The omitted conservative `Z_2,Z_4` sector has been reinstated exactly.
3. On a clean constructive branch, that block removes the fake mixed-sector freedom left in the
   Stage-17 lower-bound picture.
4. So the remaining work is no longer “some mixed-sector freedom somewhere.” It is:

> compute the actual overlap data on the true moving-throat branch and test whether they land on
> the isotropic full-bundle target surface from Stage 212.

That is the cleanest next continuation point.
