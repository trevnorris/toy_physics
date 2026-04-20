# Moving-Throat PDE — Stage 29: Exact Tracking-Branch Bounds, Monotonicity in `R_tr`, and the Residual Comparison Theorem

## Purpose

Stage 28 reduced the first coherent local D/N kernel to an exact one-parameter tracking branch:

`M_tr = G_tr(xi,delta;R_tr),`

`R_target = F_tr(xi,delta;R_tr).`

So the rank-2 closure ambiguity is gone on that branch.

The next honest question is not “what closure do we use?” It is:

> does the physical split-`U` deformation help or hurt the normalization test relative to the old flat branch?

This stage answers that exactly.

The main result is that the coherent local tracking branch is **ordered** by the tracking factor `R_tr`.

At fixed `(xi,delta)`:

- the required total loading `G_tr` decreases with `R_tr`,
- the normalized selected-branch response `F_tr` increases with `R_tr`.

Because the constructive split-`U` branch satisfies `R_tr < 1`, the physical local D/N branch therefore

- requires **more** total loading than the flat branch,
- but delivers **less** normalized response than the flat branch.

So the split-`U` deformation makes the Stage-18/19 normalization target harder, not easier.

---

## 1. Tracking-branch functions

On the coherent local D/N branch from Stage 28,

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ],`

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

The old flat branch is recovered at

`R = 1,`

so

`G_flat(xi,delta) = G_tr(xi,delta;1) = 9 xi (xi + delta) / (9 delta + 11 xi),`

`F_flat(xi,delta) = F_tr(xi,delta;1)`
`                 = (9 delta + 11 xi)^4`
`                   / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2 ].`

---

## 2. Exact monotonicity in the tracking factor `R`

Differentiate the branch functions at fixed `(xi,delta)`.

### 2.1 Loading monotonicity

The exact derivative is

`dG_tr/dR`
`= - 36 R xi^2 (delta + xi) / [ 2 R^2 xi + 9 delta + 9 xi ]^2.`

So for the physical branch `R > 0`,

`dG_tr/dR < 0.`

Thus lowering `R` increases the total loading required to realize the same softening depth `xi`.

### 2.2 Normalization monotonicity

The exact derivative is

`dF_tr/dR`
`= 4 xi (2 R xi + 9 delta + 9 xi) (2 R^2 xi + 9 delta + 9 xi) P_R`
`  / [ 81 (1 - xi) ( 2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2 )^3 ],`

with the positive polynomial

`P_R`
`= 4 R^4 xi^3`
`  + 54 R^2 delta^2 xi + 90 R^2 delta xi^2 + 36 R^2 xi^3`
`  + 162 R delta^3 + 324 R delta^2 xi + 162 R delta xi^2`
`  + 81 delta^3 + 243 delta^2 xi + 243 delta xi^2 + 81 xi^3.`

Every coefficient in `P_R` is positive, so for the physical branch `0 < xi < 1`, `delta > 0`, `R > 0`,

`dF_tr/dR > 0.`

Thus lowering `R` decreases the normalized selected-branch response at fixed `(xi,delta)`.

---

## 3. Exact comparison with the flat branch

Since the constructive split-`U` branch has `R_tr < 1`, the derivative theorems above already imply

`G_tr(xi,delta;R_tr) > G_flat(xi,delta),`

`F_tr(xi,delta;R_tr) < F_flat(xi,delta).`

But the comparison can be written in exact closed form.

### 3.1 Exact loading excess over the flat branch

`G_tr - G_flat`
`= 18 xi^2 (1 - R^2) (delta + xi)`
`  / [ (9 delta + 11 xi) (2 R^2 xi + 9 delta + 9 xi) ].`

So for `0 < R < 1`,

`G_tr - G_flat > 0.`

### 3.2 Exact normalization deficit relative to the flat branch

`F_flat - F_tr`
`= 4 xi (1 - R) P_1 P_2`
`  / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2`
`      (2 R^2 xi^2 + 9 delta^2 + 18 delta xi + 9 xi^2)^2 ],`

with

`P_1`
`= 18 R^2 delta^2 xi + 36 R^2 delta xi^2 + 22 R^2 xi^3`
`  + 81 R delta^3 + 180 R delta^2 xi + 99 R delta xi^2`
`  + 162 delta^3 + 423 delta^2 xi + 360 delta xi^2 + 99 xi^3,`

`P_2`
`= 18 R^3 delta^2 xi^2 + 36 R^3 delta xi^3 + 22 R^3 xi^4`
`  + 81 R^2 delta^3 xi + 324 R^2 delta^2 xi^2 + 459 R^2 delta xi^3 + 220 R^2 xi^4`
`  + 81 R delta^3 xi + 243 R delta^2 xi^2 + 261 R delta xi^3 + 99 R xi^4`
`  + 729 delta^4 + 3078 delta^3 xi + 4959 delta^2 xi^2 + 3600 delta xi^3 + 990 xi^4.`

Every coefficient in `P_1` and `P_2` is positive, so for `0 < R < 1`,

`F_flat - F_tr > 0.`

So the split-`U` tracking branch is **strictly below** the flat-branch normalization function at fixed `(xi,delta)`.

---

## 4. Exact endpoint formulas and bounds

The formal strong-split endpoint is `R -> 0`.

At that endpoint the branch functions simplify exactly to

`G_tr(xi,delta;0) = xi,`

`F_tr(xi,delta;0) = 1 / (1 - xi).`

So for `0 <= R <= 1` the exact bounds are

`G_flat(xi,delta) <= G_tr(xi,delta;R) <= xi,`

`1/(1 - xi) <= F_tr(xi,delta;R) <= F_flat(xi,delta).`

These are absolute tracking-branch bounds. The actual constructive branch has the tighter Stage-28 range

`1/(1+delta_U) < R_tr < 1,`

so the physical kernel sits strictly between the strong-split endpoint and the flat branch.

---

## 5. Exact residual comparison theorem

Define the tracking-branch normalization residual

`E_tr(xi) := R_target - F_tr(xi,delta;R_tr),`

and the flat-branch residual

`E_flat(xi) := R_target - F_flat(xi,delta).`

Then at fixed `(xi,delta)`,

`E_tr - E_flat = F_flat - F_tr > 0`

on the constructive split-`U` branch.

So the local split-`U` deformation worsens the normalization residual relative to the old flat branch.

Likewise the exact loading excess theorem says that the local split-`U` branch requires a larger total baseline to reach the same `xi`.

This is the sharpest exact comparison yet between the abstract Stage-18/19 branch and the first coherent concrete local kernel.

---

## 6. Best current theorem statement after Stage 29

### What is now exact

- the tracking-branch functions `G_tr` and `F_tr`,
- their exact derivatives with respect to `R`,
- the exact loading excess `G_tr - G_flat`,
- the exact normalization deficit `F_flat - F_tr`,
- the strong-split endpoint formulas,
- and the exact residual ordering `E_tr > E_flat` at fixed `(xi,delta)` on the constructive split-`U` branch.

### What this means physically

The first coherent local D/N kernel does not rescue the Stage-18/19 normalization target by hidden split-`U` assistance. It does the opposite.

The constructive split-`U` deformation

- lowers the tracking factor from `1` to `R_tr < 1`,
- increases the total loading required to reach a given softening depth,
- and decreases the normalized quadrupole response at the same point.

So the remaining theorem gap is now even sharper:

> the completed moving-throat PDE must supply kernel data strong enough to overcome this exact tracking-branch deficit, not merely to match the old flat-branch heuristic.
