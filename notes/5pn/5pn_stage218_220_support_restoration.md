# 5PN continuation notes — stages 218 through 220

These stages answer the fork opened by the Stage-217 support-blind obstruction.

The question was:

> does the failure of the support-blind mixed corridor force an isotropic-baseline pivot,
> or does it only mean that the next honest move is to restore the support carriers and
> solve the true first-order Packet-A/B master matrix on the same branch?

The answer is now sharp.

It does **not** force a pivot yet.

The support-blind corridor was killed, but the isotropic baseline itself survives.
Once the support-carrying directions are restored, the true first-order Packet-A/B master
matrix immediately regains nontrivial nullspace on the same explicit isotropic prototype.

So the correct next roadmap is:

1. restore the support variables,
2. solve the full first-order packet master matrix,
3. then ask which of those algebraic null directions can actually be realized by the
   moving-throat branch.

A baseline pivot is now the **backup** path, not the default one.

---

## Stage 218 — support-restored master matrix theorem

Files:
- `5pn_stage218_support_restored_master_matrix.py`
- `5pn_stage218_support_restored_master_matrix_output.txt`

### 1. Full true first-order packet system

Keep the same explicit positive isotropic finite-throat prototype used in Stages 215–217.
The true first-order packet-tangency system is still

```text
u2^(1) = 0,
u4^(1) = 0,
Xi_load = 0.
```

Restore the full support-carrying direction set

```text
(alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi),
```

with the orbit-lock-dependent drifts imposed exactly:

```text
alpha_deltaU = -(1+deltaU)/(1+chi0) (alpha_R + alpha_GU - alpha_GW - 2 alpha_OU),
```

```text
alpha_OW = [alpha_GW - alpha_GU + (1-E_*) alpha_OU + E_* alpha_R - (F_*/2) alpha_deltaU] / (E_* + 2),
```

```text
alpha_M = alpha_K - 2 alpha_GU + 2 alpha_OU.
```

So the full first-order Packet-A/B master matrix is the `3 x 7` Jacobian of

```text
(u2^(1), u4^(1), Xi_load)
```

with respect to those seven free directions.

### 2. Support-only 3x3 minor is already nonzero in the positive constructive quadrant

The key new theorem is that the support-only `3 x 3` minor in

```text
(alpha_K, alpha_GU, alpha_OU)
```

has exact determinant

```text
det M_support
= - pi^2 (8575 + 12717 pi^2)^2 [A_E E_* + A_F F_* + A_0]
  / [3961650000 (490 + 1503 pi^2)^6 (E_* + 2)],
```

with

```text
A_E = 36399941765600000
    + 252110109946857750 pi^2
    + 578056501787708475 pi^4
    + 433070868615970095 pi^6,
```

```text
A_F = 9400136276160000
    + 48552477588959400 pi^2
    + 110238264175381020 pi^4
    + 85279016350877148 pi^6,
```

```text
A_0 = 22826428394560000
    + 216604585787832900 pi^2
    + 597245590420119270 pi^4
    + 566382238265309238 pi^6.
```

All three coefficients are strictly positive.
So for every

```text
E_* >= 0,
F_* >= 0,
```

and in particular on the positive constructive branch `F_* > 0`,

```text
det M_support != 0.
```

### 3. Rank/nullity verdict

Because that support-only minor is already nonzero, the exact verdict is:

```text
rank(M_support) = 3,
nullity(M_support) = 2,
```

and therefore the full support-restored master matrix also has

```text
rank(M_master) = 3,
nullity(M_master) = 4.
```

So the Stage-217 obstruction has a very specific meaning:

- it kills the **support-blind** corridor,
- but it does **not** kill the current isotropic baseline.

That is the main decision result of this session.

---

## Stage 219 — support-only constructive rescue

Files:
- `5pn_stage219_support_only_constructive_rescue.py`
- `5pn_stage219_support_only_constructive_rescue_output.txt`

To make the Stage-218 theorem concrete, I then froze the representative constructive point

```text
(E_*, F_*) = (1/4, 5/6)
```

and restricted to the support-only slice

```text
alpha_GW = 0,
alpha_R  = 0.
```

So the remaining free support carriers are

```text
(beta_B, beta_varpi).
```

The exact solve of

```text
u2^(1) = 0,
u4^(1) = 0,
Xi_load = 0
```

for

```text
(alpha_K, alpha_GU, alpha_OU)
```

gives a two-parameter packet-null family:

```text
alpha_K  =  0.0898376063372746 beta_B - 0.183439324158911 beta_varpi,
alpha_GU =  0.0736940106940628 beta_B - 0.154862196153272 beta_varpi,
alpha_OU =  0.0385380345483458 beta_B - 0.105749511712718 beta_varpi.
```

So the constructive branch already contains two independent exact support-only null directions.
A convenient basis is:

### `beta_B` basis

```text
(alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi)
≈
( 0.0898376063,
  0,
  0.0736940107,
  0,
  0.0385380345,
  1,
  0).
```

### `beta_varpi` basis

```text
(alpha_K, alpha_GW, alpha_GU, alpha_R, alpha_OU, beta_B, beta_varpi)
≈
(-0.1834393242,
  0,
 -0.1548621962,
  0,
 -0.1057495117,
  0,
  1).
```

Both were checked directly in the script and satisfy

```text
u2^(1) = 0,
u4^(1) = 0,
Xi_load = 0
```

exactly.

So on the same explicit isotropic prototype that killed the support-blind corridor,
restoring just the BdG support carriers is already enough to recover a true first-order
Packet-A/B null family.

---

## Stage 220 — full constructive support-restored master solve

Files:
- `5pn_stage220_full_constructive_master_solver.py`
- `5pn_stage220_full_constructive_master_solver_output.txt`

Finally, I solved the **full** support-restored constructive master system at

```text
(E_*, F_*) = (1/4, 5/6)
```

without freezing the mixed carriers.

A convenient free-carrier choice is

```text
(alpha_GW, alpha_R, beta_B, beta_varpi),
```

with the wall/support drifts

```text
(alpha_K, alpha_GU, alpha_OU)
```

fixed linearly by the true packet-null conditions.

Numerically, the exact constructive solve is

```text
alpha_K
=  0.729738128983757 alpha_GW
 - 0.779773212960421 alpha_R
 + 0.0898376063372746 beta_B
 - 0.183439324158911 beta_varpi,
```

```text
alpha_GU
= -0.269396497721430 alpha_GW
 + 0.367875865936610 alpha_R
 + 0.0736940106940628 beta_B
 - 0.154862196153272 beta_varpi,
```

```text
alpha_OU
= -0.216658637559893 alpha_GW
 + 0.268624898608163 alpha_R
 + 0.0385380345483458 beta_B
 - 0.105749511712718 beta_varpi.
```

So the full constructive branch carries a `4`-parameter exact packet-null family.
A canonical numerical basis is:

### `alpha_GW` basis

```text
( 0.729738129,
  1,
 -0.269396498,
  0,
 -0.216658638,
  0,
  0).
```

### `alpha_R` basis

```text
(-0.779773213,
  0,
  0.367875866,
  1,
  0.268624899,
  0,
  0).
```

### `beta_B` basis

```text
( 0.0898376063,
  0,
  0.0736940107,
  0,
  0.0385380345,
  1,
  0).
```

### `beta_varpi` basis

```text
(-0.1834393242,
  0,
 -0.1548621962,
  0,
 -0.1057495117,
  0,
  1).
```

The Stage-219 support-only rescue is recovered exactly by setting

```text
alpha_GW = alpha_R = 0.
```

So the full master solve confirms the roadmap verdict from Stage 218.

---

## What changed in the roadmap

Before these stages, the honest continuation fork was:

1. restore the support carriers and solve the master matrix, or
2. pivot the isotropic baseline.

After these stages, the updated verdict is:

### Not yet time to pivot

The current isotropic baseline remains live.
The support-blind obstruction was real, but it was only a slice obstruction.

### The correct next move is now explicit

The next theorem gate is:

> determine which of the support-restored first-order packet-null directions can actually
> be realized by the moving-throat branch.

That means the next work should be about **realizability**, not about reopening the baseline choice.

Concretely, the next continuation should now attack one of two closely related tasks:

1. derive the actual branch-induced weak-axisymmetric support drifts
   `(beta_B, beta_varpi)` and mixed drifts `(alpha_GW, alpha_R)` from the moving-throat overlap data,
   then project them onto the Stage-220 null basis;
2. or derive the branch tangent map from the moving-throat reduced PDE into
   `(alpha_K, alpha_GU, alpha_GW, alpha_R, alpha_OU, beta_B, beta_varpi)`
   and test whether its image intersects the exact packet-null family.

That is a much sharper continuation point than where Stage 217 left us.

---

## Honest caveat

This is a real rescue of the current isotropic prototype, but it is still only a
**first-order reduced packet-null result**.

What is still not proved is that the actual moving-throat branch realizes any of these null directions.
So the serious remaining problem is now:

> algebraic packet-null existence is established,
> but branch realization is still open.

That is the right place for the next push.
