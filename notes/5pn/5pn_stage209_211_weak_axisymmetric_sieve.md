# 5PN continuation notes — stages 209 and 210

This session pushed the explicit isotropic overlap prototype into the first runnable
weak-axisymmetric mechanism sieve.

The clean continuation point after Stages 207–208 was to stop speaking abstractly
about “the weak-axisymmetric defect” and instead compile it directly from the same
finite-throat overlap model that already feeds Packet A and Packet B on the isotropic
branch.

The two new stages do exactly that.

---

## Stage 209 — primitive weak-axisymmetric extractor from the explicit overlap model

Files:
- `5pn_stage209_primitive_weak_axisymmetric_extractor.py`
- `5pn_stage209_primitive_weak_axisymmetric_extractor_output.txt`

### What was fixed

Start from the explicit isotropic overlap prototype already used in Stages 207–208:

- `I_cross = 8/(3 pi)`
- `C   = lambda_B I_cross`
- `G_U = lambda_U`
- `G_W = lambda_W I_cross`
- `R   = lambda_R I_cross`

and add the primitive weak-axisymmetric absolute slopes

- `dK`, `dM`
- `d(lambda_B)`, `d(varpi)`
- `d(lambda_U)`, `d(lambda_W)`, `d(lambda_R)`
- `d(Omega_U)`, `d(Omega_W)`.

From the isotropic bundle data

- `Delta = Omega_U^2 Omega_W^2 - R^2`
- `P = Omega_U^2 G_W + R G_U`
- `Q = G_U^2 Omega_W^2 + 2 G_U G_W R + G_W^2 Omega_U^2`
- `H = G_U^2 + G_W^2`
- `S2 = Omega_U^2 + Omega_W^2`

this stage derives the exact first-order primitive variations

- `Delta01`
- `P01`
- `Q01`
- `H01`
- `S201`

and then the full slope packet

- `B01,B21,B41`
- `Z01,Z21,Z41`
- `N01`
- `D01,D21,D41`.

### Exact bundle decomposition

The script proves the exact weak-axisymmetric bundle identities

```text
D01 = dK - B01 - Z01
D21 = -(dM + B21 + Z21)
D41 = -(B41 + Z41)
Xi_load = N01/N0 - D01/D0 = P1/P0
```

so the full primitive deformation only reaches the grouped first-order obstruction triplet
through the four numbers

```text
(D01, D21, D41, N01).
```

### Compact BdG slope formulas

The explicit BdG support slopes reduce to the clean logarithmic forms

```text
B01 = 2 B0 (d(lambda_B)/lambda_B - d(varpi)/varpi)
B21 = 2 B2 (d(lambda_B)/lambda_B - 2 d(varpi)/varpi)
B41 = 2 B4 (d(lambda_B)/lambda_B - 3 d(varpi)/varpi).
```

So the support part of the primitive deformation is already exactly readable.

### Exact compensation surfaces

On the canonical compensated isotropic branch, the first two first-order gates are
algebraic.

1. **Even-preserving compensation** `K1 = 0` fixes the inertia-side slope exactly:


```text
K1 = D21 + D01/9 = 0
=> dM = D01/9 - B21 - Z21 = (dK - B01 - Z01)/9 - B21 - Z21.
```

2. **Odd/normalization-preserving compensation** `Xi_load = 0` fixes the static
   wall-loading slope exactly:

```text
Xi_load = N01/N0 - D01/D0 = 0
=> dK = B01 + Z01 + D0 (N01/N0).
```

So after Stage 209 the next live first-order gate is no longer vague. It is the
remaining hidden-even residual

```text
H_even = D41 - (2/3) D21 - D01/27.
```

That is exactly the explicit first-order 5PN gate left after the two compensation
surfaces are imposed.

---

## Stage 210 — mechanism sieve and the surviving outgoing corridor

Files:
- `5pn_stage210_mechanism_sieve_and_outgoing_corridor.py`
- `5pn_stage210_mechanism_sieve_and_outgoing_corridor_output.txt`

This stage turns the Stage-209 obstruction triplet into an actual sieve.

The three first-order gates are

```text
K1      = D21 + D01/9
Xi_load = N01/N0 - D01/D0
H_even  = D41 - (2/3) D21 - D01/27.
```

### Wall-only weak-axisymmetric anisotropy is dead

With only `(dK,dM)` active,

```text
D01 = dK
D21 = -dM
D41 = 0
N01 = 0
```

so

```text
K1_wall  = dK/9 - dM
Xi_wall  = -dK/D0
H_wall   = 2 dM/3 - dK/27.
```

The exact solve of

```text
(K1_wall, Xi_wall, H_wall) = 0
```

gives only the trivial branch

```text
dK = 0,
dM = 0.
```

So wall-only weak-axisymmetric anisotropy cannot carry the first-order 5PN closure.

### Pure BdG weak-axisymmetric anisotropy is also dead

Using logarithmic support drifts

```text
x_c     = delta ln C
x_varpi = delta ln varpi
```

the exact support slopes are

```text
B01 = 2 B0 (x_c - x_varpi)
B21 = 2 B2 (x_c - 2 x_varpi)
B41 = 2 B4 (x_c - 3 x_varpi)
```

and with wall/gauge/mixed data fixed,

```text
D01 = -B01
D21 = -B21
D41 = -B41
N01 = 0.
```

The exact solve of the full first-order system

```text
(K1, Xi_load, H_even) = 0
```

again gives only the trivial branch

```text
x_c = 0,
x_varpi = 0.
```

So pure BdG weak-axisymmetric anisotropy is also ruled out.

### BdG self-similarity kills only the load defect

On the pure-BdG branch, impose the natural self-similar condition

```text
x_c = x_varpi.
```

Then the script verifies

```text
Xi_load = 0,
K1      = 2 B2 x_varpi,
H_even  = (4/3) x_varpi (-B2 + 3 B4).
```

So BdG self-similarity kills only the normalization defect. It does **not** kill the
full first-order 5PN triplet unless the branch is trivial.

### Exact one-port outgoing-load factorization

The surviving nontrivial corridor is the outgoing mixed-sector load factor.
For one port,

```text
N0 = P^2 / Delta^2
```

with

```text
P     = Omega_U^2 G_W + R G_U
Delta = Omega_U^2 Omega_W^2 - R^2.
```

The script proves the exact factorization

```text
N0 / K = M_leg^2 (1 + I)^2 / (1 - H)^2
```

where

```text
M_leg = G_W / (Omega_W^2 sqrt(K))
I     = R G_U / (Omega_U^2 G_W)
H     = R^2 / (Omega_U^2 Omega_W^2).
```

### Exact one-port outgoing defect

The corresponding outgoing defect is

```text
Sigma^(N) = 2 m + 2 I/(1+I) i + 2 H/(1-H) h
```

with

- `m` the raw mixed-leg drift,
- `i` the interference-ratio drift,
- `h` the hybridization-ratio drift.

Under rigid interference and hybridization ratios,

```text
i = 0,
h = 0,
```

this collapses to

```text
Sigma^(N) = 2 m.
```

Writing the raw mixed-leg drift as

```text
m = g_W - 2 o_W - kappa_1/2,
```

the zero-defect condition becomes

```text
2 g_W - 4 o_W - kappa_1 = 0.
```

Equivalently,

```text
G_W / Omega_W^2 ∝ sqrt(K).
```

This is the exact **square-root mixed-leg law**. It is the first nontrivial
surviving first-order cancellation condition after the wall-only and BdG-only
primitive sectors are killed.

---

## Best current status after Stages 209–210

The weak-axisymmetric continuation is now much sharper.

The live first-order problem is no longer “find some anisotropy that works.”
What survives is a very short list:

1. the full primitive deformation only enters through
   ```text
   (D01, D21, D41, N01),
   ```
2. `K1 = 0` and `Xi_load = 0` are exact algebraic compensation surfaces,
3. wall-only and pure-BdG primitive sectors are dead,
4. BdG self-similarity kills only `Xi_load`, not the full 5PN triplet,
5. the surviving nontrivial corridor is outgoing mixed-sector co-loading,
6. and under rigid interference/hybridization the remaining first-order condition is
   exactly the square-root mixed-leg law
   ```text
   G_W / Omega_W^2 ∝ sqrt(K).
   ```

So the next honest theorem gate is now:

> choose a concrete mixed-sector weak-axisymmetric mechanism, substitute it into the
> Stage-209 compiler, and test whether it realizes the outgoing co-loading law while
> also killing the hidden-even residual.

That is the clean continuation point.


---

## Stage 211 — normalized monomial bridge and exact similarity kernel

Files:
- `5pn_stage211_normalized_similarity_kernel.py`
- `5pn_stage211_normalized_similarity_kernel_output.txt`

After Stage 210, the weak-axisymmetric continuation still needed the exact bridge
back to the later monomial/similarity-orbit package. Stage 211 now makes that bridge
runnable.

### Exact normalized monomial-drift matrix

In the normalized variable set

- `dln G_W`
- `dln G_U`
- `dln R`
- `dln K`
- `dln M`
- `dln Omega_U`
- `dln Omega_W`
- `dln delta_U`

and with constructive-branch parameters `(chi_0, delta_U, E_*, F_*)`, the direct
monomial-drift matrix is exactly

```text
M_norm =
[ -(1+delta_U),  1+delta_U,  1+delta_U,  0, 0, -2(1+delta_U),        0,      1+chi_0 ]
[            2,          0,      2E_*, -1, 1,        -2E_*, -(4+2E_*),        -F_* ]
[            0,          2,          0, -1, 1,           -2,         0,           0 ]
```

The script verifies this matrix exactly and shows

```text
rank(M_norm) = 3.
```

So the normalized zero-defect tangent space has dimension `5`.

### Exact triangular zero-defect solution

The zero-defect equations solve triangularly.

1. **Tracking** fixes `dln(delta_U)`:

```text
dln(delta_U)
= - (1+delta_U)/(1+chi_0)
  [ dln R + dln G_U - dln G_W - 2 dln Omega_U ].
```

2. **Dressing** fixes `dln(M)`:

```text
dln(M) = dln(K) - 2 dln(G_U) + 2 dln(Omega_U).
```

3. **Nontracking** fixes `dln(Omega_W)`:

```text
dln(Omega_W)
= [ dln(G_W) - dln(G_U) + (1-E_*) dln(Omega_U) + E_* dln(R)
    - (F_*/2) dln(delta_U) ] / (E_* + 2).
```

So once a candidate weak-axisymmetric branch supplies the five drifts

```text
(dln K, dln G_U, dln G_W, dln R, dln Omega_U),
```

the exact similarity-orbit theorem predicts the three co-drifts required to stay on
zero defect:

```text
(dln delta_U, dln M, dln Omega_W).
```

### Support-blind extension back to the explicit Stage-5 primitive space

Extending the primitive space back to

- `dln lambda_B`
- `dln varpi`
- plus the eight normalized drifts above,

the script proves that the first two columns of the extended monomial matrix are
exact zeros.

So the direct monomial package is exactly support-blind in the two explicit BdG
primitive directions:

```text
partial_{ln lambda_B} (dln C_tr, dln C_nt, dln epsilon_eta) = 0,
partial_{ln varpi}    (dln C_tr, dln C_nt, dln epsilon_eta) = 0.
```

The extended primitive monomial matrix still has rank `3`, so its nullity is

```text
10 - 3 = 7.
```

This splits into

1. the five normalized similarity directions, and
2. the two support-blind directions `(dln lambda_B, dln varpi)`.

### Why Stage 211 matters

This is the exact splice between the two halves of the 5PN continuation:

- Stages 209–210 now isolate which primitive first-order sectors are dead and which
  mixed/outgoing corridor survives,
- while Stage 211 shows exactly how the surviving mixed/wall/U normalization problem
  is constrained by monomial rigidity.

So the current weak-axisymmetric picture is now structurally clean:

- `Xi_load = P1/P0` lives on the same mixed/wall/U normalization problem controlled by
  the similarity-orbit theorem,
- but the explicit BdG-support directions remain outside that theorem because they are
  exact zero columns of the monomial map.

That is the sharpest current continuation point.
