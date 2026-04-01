# Moving-Throat PDE — Stage 25: Selected-Mode Normalization Under Rank-2 Support Completion

## Purpose

Stage 24 closed the rank-2 *existence* problem for the selected wall branch: with mixed loading `z` and support loading `y`, the exact support load needed to reach a chosen softening depth is known.

But the 2.5PN bridge does not depend only on branch existence. It depends on the **selected-mode normalization product** that combines

- the outgoing mixed overlap,
- the source-map overlap,
- and the selected stiffness.

So the right next question is:

> how does the selected-mode normalization law deform once the branch is supported by two different directions?

This stage answers that exactly.

The main result is that the selected-mode geometry remains closed-form, but the normalization now depends on **three** directions:

- the mixed/outgoing direction `z`,
- the support direction `y`,
- the source direction `s`.

That gives an exact generalized normalization function `F_(q,r,t)(xi,delta;m)`. The old Stage-23 law is recovered as the special case `r=q`.

So the rank-2 bottleneck is no longer conceptual. It is now a sharply defined source/support/outgoing alignment problem.

---

## 1. Exact selected-mode geometry with two loading directions

Carry forward the Stage-24 wall operator

`K_loaded / A_0 = diag(1,1+delta) - m (1,q)^T(1,q) - n (1,r)^T(1,r),`

with the lower selected eigenvalue parameterized by

`lambda_- = A_0 (1 - xi),`

and with the support loading fixed to the exact Stage-24 value

`n = n_req(xi,delta;m,q,r).`

Then the lower selected eigenvector ratio is still exact:

`e_1 / e_0`
`= [ m(q-r) + r xi ] / [ delta + xi - m q(q-r) ].`

So the two-direction selected branch is still described by one scalar deformation variable `xi`, but the eigenvector now depends explicitly on the baseline mixed loading `m` whenever the support direction differs from the mixed direction.

That dependence was absent in the Stage-23 one-direction problem.

---

## 2. Exact overlap formulas

Now introduce a general source vector

`s = s_0 (1,t)^T.`

Using the exact selected eigenvector above, the normalized overlaps are

`(z.e_-)^2 / z_0^2`
`= [ delta + (1 + q r) xi ]^2`
`  / D_(q,r)(xi,delta;m),`

`(s.e_-)^2 / s_0^2`
`= [ delta + (1 + r t) xi - m (q-r)(q-t) ]^2`
`  / D_(q,r)(xi,delta;m),`

with the exact denominator

`D_(q,r)(xi,delta;m)`
`= [ delta + xi - m q(q-r) ]^2 + [ m(q-r) + r xi ]^2.`

So the normalization data no longer factor through one direction only. The mixed/outgoing overlap and the source overlap are deformed differently once `q != r`.

---

## 3. Exact rank-2 normalization function

The selected-branch normalization product is therefore

`F_(q,r,t)(xi,delta;m)`
`= [ delta + (1 + q r) xi ]^2`
`  [ delta + (1 + r t) xi - m(q-r)(q-t) ]^2`
`  / [ (1 - xi) D_(q,r)(xi,delta;m)^2 ].`

This is the exact rank-2 generalization of the Stage-23 selected-mode normalization law.

It is the object that the outgoing quadrupole bridge now depends on whenever the support direction and the mixed direction differ.

---

## 4. Exact collapse to Stage 23 when support tracks the mixed vector

If the support tracks the mixed direction,

`r = q,`

then

`D_(q,q) = (delta + xi)^2 + q^2 xi^2,`

and the exact rank-2 normalization function collapses to

`F_(q,q,t)(xi,delta;m)`
`= [ delta + (1 + q^2) xi ]^2 [ delta + (1 + q t) xi ]^2`
`  / [ (1 - xi) ((delta + xi)^2 + q^2 xi^2)^2 ].`

So the second sharp theorem is:

> **If the support/BdG loading follows the deformed mixed direction, the entire rank-2 normalization law collapses exactly to the Stage-23 two-vector function.**

In that case the rank-2 completion adds no new normalization geometry.

---

## 5. Source-tied specialization for the physical split-`U` continuum

Now take the physically interesting opposite hypothesis:

- support direction tied to the original source vector,
- source vector itself also equal to that original D/N direction.

Write

`t := kappa_1 / kappa_0,`

so that

`t^2 = lambda_0 = 2/9.`

From Stage 22 the mixed direction obeys

`q = t R_U,`

while the source-tied support hypothesis gives

`r = t.`

Then the exact normalization function becomes

`F_src(xi,delta;m,R_U)`
`= [ delta + (1 + lambda_0 R_U) xi ]^2`
`  [ delta + (1 + lambda_0) xi - m lambda_0 (R_U - 1)^2 ]^2`
`  / [ (1 - xi) D_src(xi,delta;m,R_U)^2 ],`

with

`D_src(xi,delta;m,R_U)`
`= [ delta + xi - m lambda_0 R_U (R_U - 1) ]^2`
`  + lambda_0 [ xi + m (R_U - 1) ]^2.`

This is the first exact selected-mode normalization law beyond Stage 23 that depends explicitly on the mixed baseline loading `m`.

That dependence is the clean reduced signature of the source/support/outgoing mismatch.

---

## 6. Exact flat-`U` recovery and first-order deformation

When `R_U = 1`, all three directions become collinear and the source-tied rank-2 completion collapses exactly to the old flat-`U` law:

`F_src(xi,delta;m,1) = F_flat(xi,delta),`

where

`F_flat(xi,delta)`
`= [ delta + (1 + lambda_0) xi ]^4`
`  / [ (1 - xi) ( (delta + xi)^2 + lambda_0 xi^2 )^2 ].`

So the source-tied rank-2 branch is a genuine deformation of Stage 23, not a different disconnected object.

Write

`R_U = 1 + eps.`

Then the exact first-order support-loading deformation is

`n_src = G_flat(xi,delta) - m + eps H_n^(src) + O(eps^2),`

with

`H_n^(src)`
`= - 2 lambda_0 m xi / [ delta + (1 + lambda_0) xi ].`

Likewise the exact first-order normalization deformation is

`F_src / F_flat = 1 + eps H_F^(src) + O(eps^2),`

with

`H_F^(src)`
`= 2 lambda_0`
`  [ xi ( (delta + xi)^2 + lambda_0 xi^2 )`
`    + 2 m delta ( delta + (1 + lambda_0) xi ) ]`
`  / [ ( delta + (1 + lambda_0) xi ) ( (delta + xi)^2 + lambda_0 xi^2 ) ].`

On the constructive split-`U` branch of Stage 22, one has `R_U < 1`, i.e. `eps < 0`, so:

- the source-tied support requirement is **raised** above the simple flat subtraction `G_flat - m`,
- the selected-mode normalization function is **lowered** below the flat-`U` value.

So the source-tied hypothesis makes the selected-branch normalization test strictly harder on the natural split-`U` branch.

---

## 7. Best current theorem statement after Stage 25

The rank-2 support bottleneck is now completely explicit.

### If support tracks the mixed vector

- Stage 24 gives `n_req = G_q - m`,
- Stage 25 gives exact collapse back to the Stage-23 normalization law.

So the old one-parameter deformation survives intact.

### If support remains tied to the original source direction

- Stage 24 gives the new exact source-tied support-feasibility formula,
- Stage 25 gives the new exact source-tied normalization function `F_src`.

So the selected branch becomes a true two-layer deformation:

1. Stage 22 deforms the mixed direction through `R_U`,
2. Stages 24–25 decide whether the support kernel follows that deformation or resists it.

The open PDE-side theorem question is therefore no longer diffuse:

> compute the actual projected support/BdG loading direction from the moving-throat operator and determine whether the physical kernel lands on the tracking closure or the source-tied closure.

That is now the sharp reduced theorem gate for the next stage.
