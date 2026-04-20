# Moving-Throat PDE — Stage 24: Rank-2 Support Completion and the Exact Support-Loading Theorem

## Purpose

Stage 23 showed that once the first axial structure of the internal `U` sector is turned on, the mixed loading vector `z` is no longer generically collinear with the source/support direction `v`.

That leaves one honest reduced-theorem bottleneck:

> what happens when the selected wall branch is loaded by **two** directions rather than one?

This stage answers that exactly.

The main result is that the rank-2 completion is still analytically tractable because the determinant stays **linear** in the support loading. That gives an exact support-loading theorem.

More concretely:

1. the selected lower branch with mixed loading `z` and support loading `y` is still solvable in closed form,
2. the exact required support loading is a rational function `n_req(xi,delta;m,q,r)`,
3. increasing the mixed baseline always lowers the support needed to hit a given selected branch,
4. if the support tracks the deformed mixed direction, the whole rank-2 problem collapses back to Stage 23,
5. if the support stays tied to the original source direction, the branch geometry becomes genuinely new, with an exact support-feasibility window.

So the flat-`U` simplification did not just hide a numerical correction. It hid the first place where the support direction itself becomes a real theorem variable.

---

## 1. Exact rank-2 loaded selected branch

Work in the 2D wall basis of Stages 22–23 and write the diagonal baseline as

`K_base = A_0 diag(1, 1 + delta).`

Now allow two independent rank-1 loadings:

- a mixed baseline loading carried by `z`,
- a support/BdG loading carried by `y`.

Normalize the first wall-basis component of each direction and define

`z = z_0 (1, q)^T,`

`y = y_0 (1, r)^T,`

with dimensionless loadings

`m := alpha_mix z_0^2 / A_0,`

`n := alpha_supp y_0^2 / A_0.`

The full selected-branch matrix is therefore

`K_loaded / A_0 = diag(1,1+delta) - m (1,q)^T(1,q) - n (1,r)^T(1,r).`

As before, parameterize the lower selected eigenvalue by the softening depth

`lambda_- = A_0 (1 - xi),`

with stable branch parameter

`0 <= xi < 1.`

Then the exact determinant condition is

`D_sel = xi (delta + xi)`
`        - m [ delta + (1 + q^2) xi ]`
`        - n [ delta + (1 + r^2) xi ]`
`        + m n (q - r)^2`
`      = 0.`

This is the structural point of the stage:

> even with two loading directions, the determinant remains **linear** in the support loading `n`.

So the rank-2 completion is not a high-dimensional numerical problem yet. It is still an exact reduced theorem.

---

## 2. Exact support-loading theorem

Solving the determinant condition for the support loading gives

`n_req(xi,delta;m,q,r)`
`= [ xi(delta + xi) - m( delta + (1 + q^2) xi ) ]`
`  / [ delta + (1 + r^2) xi - m (q - r)^2 ].`

This is the exact rank-2 support-loading theorem.

The denominator shows the first genuinely new effect of noncollinearity:

`delta + (1 + r^2) xi - m (q - r)^2.`

If the two directions differ, the mixed baseline modifies the *support* denominator directly. That effect was absent in the collinear Stage-18/19/23 branch.

A second exact theorem follows immediately by differentiation:

`d n_req / d m`
`= - [ delta + (1 + q r) xi ]^2`
`  / [ delta + (1 + r^2) xi - m (q - r)^2 ]^2`
`< 0.`

So, whenever the branch is regular, **increasing the mixed baseline always lowers the support loading needed to reach the same softening depth**.

That monotonicity is exact and independent of the detailed direction mismatch.

---

## 3. Tracking theorem: support follows the deformed mixed direction

The first natural hypothesis is

`y || z`,

so that

`r = q.`

Then the exact support theorem collapses immediately to

`n_req = xi (delta + xi) / [ delta + (1 + q^2) xi ] - m`
`      = G_q(xi,delta) - m.`

This is exactly the old Stage-19/23 geometry with a split between

- mixed baseline loading `m`,
- and additional support loading `n`.

So the first sharp theorem is:

> **If the physical support/BdG loading follows the deformed mixed direction, the rank-2 completion introduces no new selected-branch geometry.**

The whole problem collapses back to the Stage-23 one-parameter deformation.

That makes `y || z` a very strong closure hypothesis.

---

## 4. Source-tied theorem: support remains aligned with the original source direction

The second natural hypothesis is the opposite one:

`y || v,`

where `v = (kappa_0, kappa_1)^T` is the original D/N source/support direction.

Write the source ratio as

`t := kappa_1 / kappa_0,`

so that

`t^2 = lambda_0 = 2/9.`

From Stage 22, the mixed vector obeys

`q = t R_U,`

where `R_U` is the exact split-`U` direction ratio. Under the source-tied hypothesis,

`r = t.`

The exact required support loading becomes

`n_req^(src)`
`= [ xi(delta + xi) - m( delta + (1 + lambda_0 R_U^2) xi ) ]`
`  / [ delta + (1 + lambda_0) xi - m lambda_0 (R_U - 1)^2 ].`

This is the first genuinely new branch formula beyond Stage 23.

Two exact support-feasibility conditions follow.

### Regularity condition

The branch stays regular only if

`delta + (1 + lambda_0) xi - m lambda_0 (R_U - 1)^2 > 0,`

or equivalently

`m < [ delta + (1 + lambda_0) xi ] / [ lambda_0 (R_U - 1)^2 ]`

for `R_U != 1`.

### Positive-support condition

Assuming the regularity denominator is positive, nonnegative support loading requires

`m <= [ xi(delta + xi) ] / [ delta + (1 + lambda_0 R_U^2) xi ].`

So the source-tied rank-2 branch has a sharp mixed-loading ceiling before the support channel can even remain physical.

That ceiling disappears in the collinear flat-`U` limit `R_U = 1`, where the denominator correction vanishes and the branch collapses back to Stage 19.

---

## 5. Exact comparison of the two hypotheses

The rank-2 support completion therefore splits into two sharply testable reduced hypotheses.

### Hypothesis A — support tracks the mixed vector

`y || z`

Then

`n_req = G_q - m,`

and no new branch geometry appears.

### Hypothesis B — support stays tied to the original source vector

`y || v`

Then

`n_req = n_req^(src)`

with the exact denominator correction

`- m lambda_0 (R_U - 1)^2`.

This is the first place where the noncollinearity of Stage 22 becomes a true structural effect instead of just a deformation of Stage-23 constants.

So the exact theorem picture is now:

- **tracking support** preserves the Stage-23 geometry,
- **source-tied support** creates a new rank-2 support-feasibility problem.

That is the sharpest reduced statement yet of the support-direction bottleneck.

---

## 6. Best current theorem statement after Stage 24

The rank-2 support completion is no longer vague.

### What is now exact

- the selected-branch determinant with two loading directions,
- the exact required support loading `n_req`,
- the exact monotonicity theorem `d n_req / d m < 0`,
- the exact collapse to Stage 23 when `y || z`,
- and the exact source-tied branch formula when `y || v`.

### What remains open

The unresolved PDE-side question is now very specific:

> does the physical support/BdG kernel align with the deformed mixed vector `z`, or does it remain tied to the original source direction `v`?

If it aligns with `z`, the Stage-23 one-parameter deformation is already the whole selected-branch geometry.
If it stays tied to `v`, the branch acquires a new exact rank-2 support-feasibility window.

That is the next theorem gate the completed moving-throat operator has to decide.
