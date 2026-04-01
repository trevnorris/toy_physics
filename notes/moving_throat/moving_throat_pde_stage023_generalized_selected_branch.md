# Moving-Throat PDE — Stage 23: Generalized Selected-Branch Normalization with Source/Loading Mismatch

## Purpose

Stage 22 showed that the first axial structure of the internal `U` sector does **not** destroy the scalar continuum placement map, but it **does** generically rotate the mixed loading vector away from the source/support direction.

That means the old Stage-18/19 branch functions cannot be assumed blindly once `delta_U != 0`.

The right next question is therefore:

> how does the selected-branch normalization law deform when the source vector and the loading vector are no longer the same?

This stage answers that exactly.

The main result is that the selected-branch geometry survives, but it becomes a one-parameter deformation of the flat-`U` branch.

More concretely:

1. the selected lower wall mode can still be solved exactly,
2. the old normalization function `F(xi,delta)` is replaced by a two-vector function `F_(q,eta)(xi,delta)`,
3. the required baseline loading remains one-dimensional through a deformed function `G_q(xi,delta)`,
4. for the split-`U` continuum the whole deformation collapses to one exact ratio `R_U`,
5. and setting `R_U = 1` reproduces the Stage-18/19 flat-`U` branch exactly.

So the first non-flat `U` structure does not kill the theorem geometry — it deforms it in a controlled way.

---

## 1. Exact rank-1 loaded 2x2 branch solve

Consider the selected wall basis after the Stage-22 direct `U` softening has already been absorbed into the diagonal baseline matrix

`K_base = diag(A_0, A_0(1 + delta)).`

Now add a single directional loading

`K_loaded = K_base - alpha z z^T,`

with loading vector

`z = (z_0, z_1)^T.`

Write the lower selected eigenvalue as

`lambda_- = A_0 (1 - xi),`

with stable branch parameter

`0 <= xi < 1.`

Then the exact required loading is

`alpha_req = A_0 xi (xi + delta) / [ z_0^2 ( delta + (1 + q^2) xi ) ],`

where

`q := z_1 / z_0.`

So the basic one-parameter selected-branch geometry survives exactly even before any source map is specified.

The exact lower-mode eigenvector ratio is

`e_1 / e_0 = q xi / (delta + xi).`

That is the universal 2x2 rank-1-loaded backbone behind everything that follows.

---

## 2. Exact overlap formulas with distinct source and loading directions

Now separate two vectors that coincided in the flat-`U` stages:

- the **loading vector** `z`,
- the **source vector** `s`.

Define the signed mismatch parameter

`eta := (s_1 z_1) / (s_0 z_0).`

Then the exact normalized overlaps of the selected lower mode are

`(z.e_-)^2 / z_0^2`
`= [ delta + (1 + q^2) xi ]^2`
`  / [ (delta + xi)^2 + q^2 xi^2 ],`

`(s.e_-)^2 / s_0^2`
`= [ delta + (1 + eta) xi ]^2`
`  / [ (delta + xi)^2 + q^2 xi^2 ].`

So the selected-branch normalization product is no longer governed by the old single-vector function.
It is governed by the exact two-vector shape function

`F_(q,eta)(xi,delta)`
`= [ delta + (1 + q^2) xi ]^2 [ delta + (1 + eta) xi ]^2`
`  / [ (1 - xi) ( (delta + xi)^2 + q^2 xi^2 )^2 ].`

At the same time, the required baseline loading stays one-dimensional:

`G_q(xi,delta) = xi (xi + delta) / [ delta + (1 + q^2) xi ].`

This is the exact generalization of the Stage-18/19 pair `(F,G)`.

---

## 3. Specialization to the split-`U` continuum

For the split-`U` continuum of Stage 22, the source vector is still the original D/N overlap direction

`s = v = (kappa_0, kappa_1)^T,`

while the loading vector obeys

`z_1 / z_0 = (kappa_1 / kappa_0) R_U.`

Using

`lambda_0 := kappa_1^2 / kappa_0^2 = 2/9,`

this means

`q = - sqrt(lambda_0) R_U = - (sqrt(2)/3) R_U,`

and therefore

`eta = (s_1 z_1)/(s_0 z_0) = lambda_0 R_U = (2/9) R_U.`

So the full source/loading mismatch collapses to one exact parameter `R_U`.

Substituting that into the general formulas gives the exact split-`U` branch functions

`F_U(xi,delta;R_U)`
`= [ 9 delta + (9 + 2 R_U^2) xi ]^2 [ 9 delta + (9 + 2 R_U) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R_U^2) xi^2 )^2 ],`

`G_U(xi,delta;R_U)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R_U^2) xi ].`

This is the exact one-parameter deformation of the old flat-`U` selected-branch geometry.

---

## 4. Exact recovery of the Stage-18/19 flat-`U` branch

If `R_U = 1`, then `q = kappa_1/kappa_0` and the source/loading mismatch disappears.
The exact branch functions collapse back to the Stage-18/19 formulas:

`F_U(xi,delta;1)`
`= (9 delta + 11 xi)^4 / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2 ],`

`G_U(xi,delta;1)`
`= 9 xi (xi + delta) / (9 delta + 11 xi ).`

So the flat-`U` selected-branch theorem is not replaced. It is recovered as the exact collinear limit of the more general split-`U` theory.

---

## 5. Exact small-deformation expansion about the flat-`U` branch

Write

`R_U = 1 + eps,`

with `eps` small.

Then the exact normalized deformation of the selected-branch functions is

`F_U / F_flat = 1 + eps H_F + O(eps^2),`

`G_U / G_flat = 1 + eps H_G + O(eps^2),`

with

`H_F`
`= 4 xi ( 27 delta^2 + 36 delta xi + 11 xi^2 )`
`  / [ (9 delta + 11 xi) (9 delta^2 + 18 delta xi + 11 xi^2) ],`

`H_G = - 4 xi / (9 delta + 11 xi).`

So the deformation is smooth and exact.

Now combine this with the Stage-22 result

`R_U = 1 - [ rho_0 / (1 + rho_0) ] delta_U + O(delta_U^2).`

On the natural constructive branch `rho_0 > 0`, positive internal axial splitting implies `R_U < 1`, so:

- `F_U` is shifted **below** the flat-`U` normalization function at fixed `(xi,delta)`,
- `G_U` is shifted **above** the flat-`U` baseline-loading function at fixed `(xi,delta)`.

That is the first exact robustness statement for the selected branch once source/loading collinearity is relaxed.

---

## 6. Best current theorem statement after Stage 23

The exact theorem picture is now much sharper.

### What survives intact

- the Stage-22 scalar continuum placement map,
- the existence of a one-parameter selected branch `xi`,
- the exact lower-eigenvalue softening law,
- and the exact recovery of the Stage-18/19 functions when the split-`U` deformation is removed.

### What changes structurally

- the selected-branch normalization no longer depends on one vector only,
- the Stage-18 function `F(xi,delta)` is replaced by `F_(q,eta)(xi,delta)`,
- and in the physical split-`U` continuum this becomes a one-parameter deformation `F_U(xi,delta;R_U)`.

### The new bottleneck

The flat-`U` program reduced the selected branch to one collinear direction shared by

- source,
- support,
- and mixed loading.

Stages 22–23 show that the first non-flat `U` structure breaks that coincidence generically.

So the next honest derivation step is no longer “compute more scalar placement ratios.”
Those are already exact.

The next step is:

> determine how the additional support/BdG loading enters the now noncollinear selected-branch geometry, and whether the physical support direction tracks the deformed loading vector or remains tied to the original source vector.

That is the first place where the old Stage-19 support-feasibility theorem may need a true rank-2 extension.
