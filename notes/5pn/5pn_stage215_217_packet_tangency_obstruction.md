# 5PN continuation notes — stages 215 through 217

This session tightened the weak-axisymmetric 5PN program in a way that changes the next honest roadmap.

The key issue was that the recent mixed-sector continuation had still been phrased in the surrogate even-gate variables

- `K1 = D21 + D01/9`,
- `H_even = D41 - (2/3) D21 - D01/27`,
- `Xi_load = N01/N0 - D01/D0 = P1/P0`,

whereas the true reduced endgame from Stages 199–201 is formulated in the grouped Packet-A / Packet-B data.

So the first question for this session was:

> when does a weak-axisymmetric microscopic deformation actually remain tangent to the *true* first-order Packet-A/B endgame, rather than only solving a lower-bound surrogate gate?

That question is now answered exactly.

---

## Stage 215 — exact first-order Packet-A tangency theorem

Files:
- `5pn_stage215_first_order_packet_tangency_theorem.py`
- `5pn_stage215_first_order_packet_tangency_theorem_output.txt`

### Setup

Take an isotropic grouped-lane baseline with

```text
D_A0 = D0,
D_A2 = D2,
D_A4 = D4,
N_A0 = N0,
```

and add a weak-axisymmetric grouped perturbation with the exact signature

```text
lambda_(20) = 1,
lambda_(21) = 1/2,
lambda_(22) = -1.
```

So the perturbed grouped-lane data are

```text
D_A0 = D0 + epsilon lambda_A D01
D_A2 = D2 + epsilon lambda_A D21
D_A4 = D4 + epsilon lambda_A D41
N_A0 = N0 + epsilon lambda_A N01.
```

### Exact first-order grouped anomalies

Compiling the grouped response moments

```text
u2^(A) = -D_A2 / D_A0
u4^(A) = (D_A2^2 - D_A0 D_A4) / D_A0^2
P0^(A) = N_A0 / D_A0,
```

the script proves two exact facts.

First, the weighted grouped trace is invisible at first order:

```text
d/d epsilon ubar(u2)|0 = 0

d/d epsilon ubar(u4)|0 = 0

d/d epsilon ubar(P0)|0 = 0.
```

So on a genuine weak-axisymmetric grouped branch, the scalar trace defects
`Delta_pole` and `Delta_norm` are automatically invisible at `O(epsilon)`.

Second, the live first-order grouped anisotropy slots are one-scalar each:

```text
a(u2) =  epsilon u2^(1) / 4,
 b(u2) = 3 epsilon u2^(1) / 4,

a(u4) =  epsilon u4^(1) / 4,
 b(u4) = 3 epsilon u4^(1) / 4,

a(P0) =  epsilon P1 / 4,
 b(P0) = 3 epsilon P1 / 4.
```

So the entire first-order Packet-A tangency problem collapses to the three scalars

```text
u2^(1),  u4^(1),  P1.
```

### Exact operator formulas

The script derives the exact operator-slope formulas

```text
u2^(1) = (-D0 D21 + D2 D01) / D0^2
       = -(D21 + u2 D01) / D0
```

```text
u4^(1)
= [ D0(-D0 D41 + 2 D2 D21 - D4 D01) + 2 D01(D0 D4 - D2^2) ] / D0^3
= -(D41 + 2 u2 D21 + (u2^2 + u4) D01) / D0
```

and

```text
P1/P0 = N01/N0 - D01/D0 = Xi_load.
```

So the **true** first-order packet tangency conditions are exactly

```text
u2^(1) = 0,
nu4^(1) = 0,
P1 = 0.
```

### One-pole specialization

On an isotropic one-pole baseline `u4 = 4 u2^2`,

```text
u2^(1) = 0  ->  D21 = -u2 D01.
```

If this holds, then the second condition reduces to

```text
u4^(1) = 0  ->  D41 = -3 u2^2 D01.
```

So the exact one-pole tangent conditions are

```text
D21 = -u2 D01,
D41 = -3 u2^2 D01,
N01/N0 = D01/D0.
```

### Canonical compensated branch and the surrogate gates

On the **canonical compensated branch** used in the earlier notes,

```text
u2 = 1/9,
u4 = 4/81,
```

so the exact true tangency conditions become

```text
D21 = -D01/9,
D41 = -D01/27,
N01/N0 = D01/D0.
```

The script then proves the exact translation

```text
u2^(1) = -K1 / D0
```

and

```text
nu4^(1) = -(H_even + (8/9) K1) / D0.
```

Therefore, on the canonical compensated branch,

```text
u2^(1) = 0,
nu4^(1) = 0,
P1 = 0
```

is exactly equivalent to

```text
K1 = 0,
H_even = 0,
Xi_load = 0.
```

### Most useful Stage-215 result

This stage cleanly separates two things that had been easy to conflate:

1. **true first-order Packet-A tangency**, and
2. the earlier **surrogate even-gate variables**.

They coincide exactly on the canonical compensated branch, but not on a generic isotropic prototype with a different baseline `u2`.

That turns out to matter immediately in the next stage.

---

## Stage 216 — the support-blind mixed sector and the exact packet-null line

Files:
- `5pn_stage216_support_blind_packet_null_condition.py`
- `5pn_stage216_support_blind_packet_null_condition_output.txt`

### Explicit isotropic prototype branch

I then returned to the explicit positive isotropic finite-throat overlap prototype built from

```text
I_cross = 8 / (3 pi)
lambda_B = 1
varpi    = 2
G_U      = 1
G_W      = 1
R        = 1/2
Omega_U  = 3/2
Omega_W  = 2
M        = 1,
```

with `K` fixed on the exact isotropic one-pole surface

```text
D0 (B4 + Z4) = 3 (M + B2 + Z2)^2.
```

The first critical fact is that this explicit prototype does **not** sit on the canonical compensated branch. Its actual baseline value is

```text
u2 = (8575 + 12717 pi^2) / [210 (490 + 1503 pi^2)] ≈ 0.0416671714,
```

so

```text
u2 - 1/9 != 0.
```

Therefore the correct first-order packet test on this prototype is not the surrogate pair `(K1,H_even)`. It is the true Stage-215 packet system

```text
u2^(1) = 0,
nu4^(1) = 0,
Xi_load = 0.
```

### Support-blind mixed sector

I then imposed the same support-blind mixed-sector restriction at the input level:

- free drifts:
  ```text
  (alpha_K, alpha_GW, alpha_R),
  ```
- frozen support / upper-leg directions:
  ```text
  alpha_GU = 0,
  alpha_OU = 0,
  beta_B = 0,
  beta_varpi = 0,
  ```
- dependent co-drifts from the normalized orbit relations:
  ```text
  alpha_deltaU = -(1+deltaU)/(1+chi0) (alpha_R - alpha_GW),
  alpha_OW     = [alpha_GW + E_* alpha_R - (F_*/2) alpha_deltaU] / (E_* + 2).
  ```

So the support-blind mixed family is still parameterized by the three free directions

```text
(alpha_K, alpha_GW, alpha_R),
```

but now the **true** packet-null problem is the homogeneous `3 x 3` linear system generated by

```text
u2^(1),
nu4^(1),
Xi_load.
```

### Exact determinant and the packet-null line

The script computes the exact determinant of that `3 x 3` matrix and proves that it factors to a single affine line in the orbit constants `(E_*,F_*)`:

```text
det(M_packet) = 0
iff
F_* = F_crit(E_*),
```

with exact affine law

```text
F_crit(E_*) = -(A_E E_* + A_0) / A_F,
```

where

```text
A_E = 263797293760000
    + 1757766806455275 pi^2
    + 3339557838723645 pi^4
    + 1551622258297188 pi^6,

A_F = 48655861632000
    + 389171318788980 pi^2
    + 930178880126748 pi^4
    + 694451446430976 pi^6,

A_0 = 102703468791960 pi^6
    - 155911749769062 pi^4
    - 147002028439770 pi^2
    - 25886544768000.
```

So the support-blind mixed sector does **not** contain a generic true first-order Packet-A/B null direction.
It contains one only on the affine codimension-1 line

```text
F_* = F_crit(E_*).
```

### Concrete rank checks

The script then evaluates two concrete sample points.

At the representative constructive point already used in the notes,

```text
(E_*, F_*) = (1/4, 5/6),
```

the exact packet-null matrix has

```text
rank = 3,
```

so the support-blind packet-null system is **trivial** there.

At the illustrative point

```text
(E_*, F_*) = (1, 1),
```

the matrix again has

```text
rank = 3.
```

So the tempting support-blind mixed corridor that survived the lower-bound surrogate even gates does **not** survive the true first-order packet-null test on the explicit isotropic prototype at these positive orbit values.

### Most useful Stage-216 result

This is the first exact executable proof that the support-blind mixed corridor is much narrower than the lower-bound even-gate picture suggested.

It is not generically a first-order Packet-A/B null direction.
It exists only on one affine line in `(E_*,F_*)`.

---

## Stage 217 — positive-quadrant obstruction and the updated roadmap

Files:
- `5pn_stage217_positive_quadrant_obstruction.py`
- `5pn_stage217_positive_quadrant_obstruction_output.txt`

### Exact sign theorem for the packet-null line

Stage 217 turns the affine line from Stage 216 into a real obstruction theorem.

The exact line is

```text
F_* = F_crit(E_*),
```

and the script proves that both its slope and intercept are negative:

```text
dF_crit/dE_* < 0,
F_crit(0) < 0.
```

Numerically,

```text
F_crit(E_*) ≈ -0.1076895540 - 2.4072204236 E_*.
```

In particular,

```text
F_crit(0)   ≈ -0.1076895540,
F_crit(1/4) ≈ -0.7094946599,
F_crit(1/2) ≈ -1.3112997658,
F_crit(1)   ≈ -2.5149099776.
```

So for every `E_* >= 0`,

```text
F_crit(E_*) < 0.
```

### Consequence

The support-blind packet-null line never enters the constructive positive quadrant

```text
E_* >= 0,
F_* > 0.
```

So the support-blind mixed sector does admit a **mathematical** first-order packet-null corridor,
but only outside the positive constructive-orbit branch.

### Explicit mathematical null vector outside the constructive quadrant

To show the obstruction is physical rather than algebraic, the script evaluates one point on the line, namely

```text
E_* = 1/4,
F_* = F_crit(1/4) ≈ -0.7094946599.
```

There the packet-null matrix has

```text
rank = 2,
```

so a genuine one-parameter null corridor exists. A normalized null vector is

```text
(alpha_K, alpha_GW, alpha_R)
≈ (0.2698838731, 0.6492432127, 1).
```

So the support-blind mixed corridor is not algebraically impossible.
It is simply pushed out of the positive constructive-orbit branch.

### Most useful Stage-217 result

This is the new roadmap verdict:

1. the earlier support-blind mixed corridor was only a lower-bound even-gate survivor,
2. once the **true** first-order Packet-A/B tangency conditions are imposed on the explicit isotropic prototype,
   that corridor survives only on a negative-`F_*` affine line,
3. so it does **not** survive on the positive constructive-orbit branch.

That means the next honest 5PN continuation is now much sharper:

> restore the support-carrying directions
> `(alpha_GU, alpha_OU, beta_B, beta_varpi)`
> and recompute the full first-order packet-null matrix,
> or else move to a different isotropic baseline branch before re-testing the endgame.

This is a real change in the continuation point.
The support-blind mixed sector is no longer the default favorite corridor once the actual first-order Packet-A/B endgame is enforced.
