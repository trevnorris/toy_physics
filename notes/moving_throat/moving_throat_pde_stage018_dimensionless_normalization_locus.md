# Moving-Throat PDE — Stage 18: Dimensionless D/N Shape Function, Unique Normalization Locus, and Exact Required Support Coupling

## Purpose

Stage 17 rewrote the selected-branch normalization problem as one scalar equation in the softening depth `x`.
That was already a real simplification, but the formulas still carried the dimensional parameters `(A, DeltaK_ax, beta_0)` in a way that partly obscured the geometry of the branch.

The next step is therefore to pass to the exact D/N finite-throat constants and reduce the entire problem to a **dimensionless shape function**.

This is worthwhile because it gives the strongest selected-branch result so far:

1. the whole normalization problem collapses to one equation
   `F(xi,delta) = R_target`,
2. the shape function `F` is strictly monotone on the stable branch,
3. it maps the stable branch exactly from `1` to `+infinity`,
4. so the normalization locus is unique whenever the onset condition is satisfied,
5. and the exact required total loading and required support loading then follow immediately.

So Stage 18 turns the microscopic normalization problem into a one-dimensional universal branch-geometry problem plus a final support-coupling feasibility check.

---

## 1. Exact D/N constants and dimensionless variables

Insert the exact finite-throat D/N overlap constants,

`kappa_0^2 = 8 / pi^2`,

`kappa_1^2 = 16 / (9 pi^2)`,

so the exact overlap ratio is

`eta := kappa_1^2 / kappa_0^2 = 2/9`.

Now define dimensionless selected-branch variables

`xi := x / A`,

`delta := DeltaK_ax / A`.

Then stable selected branches satisfy

`0 <= xi < 1`,

while `delta > 0` is the axial anisotropy of the bare wall block in units of the flat stiffness.

---

## 2. Exact dimensionless shape function

Using the Stage-17 softening-depth normal form, the exact normalization product may be written as

`N_-(x) = N_-(0) F(xi,delta)`,

with

`N_-(0) = beta_0 kappa_0^2 / A`.

After inserting the D/N constants, the exact dimensionless shape function is

`F(xi,delta)`
`= (9 delta + 11 xi)^4`
`  / [ 81 (1 - xi) (9 delta^2 + 18 delta xi + 11 xi^2)^2 ]`.

So the full selected-branch target becomes

`F(xi,delta) = R_target`,

where

`R_target := N_Q^(target) / N_-(0)`
`          = N_Q^(target) A / ( beta_0 kappa_0^2 )`.

This is the cleanest reduced theorem form yet.
All dependence on the microscopic transfer factor and the overall wall scale has collapsed into the single dimensionless target ratio `R_target`.

---

## 3. Exact monotonicity, onset, and softening limits

Differentiate the exact D/N shape function.
One obtains

`dF/dxi`
`= (9 delta + 11 xi)^3`
`  [ 81 delta^3 + 72 delta^2 + 189 delta^2 xi + 297 delta xi^2 + 121 xi^3 ]`
`  / [ 81 (1 - xi)^2 (9 delta^2 + 18 delta xi + 11 xi^2)^3 ]`
`> 0`.

So `F` is strictly monotone increasing on the entire stable branch `0 <= xi < 1`.

Its exact endpoint values are

`F(0,delta) = 1`,

`lim_{xi -> 1^-} F(xi,delta) = +infinity`.

More precisely,

`F(xi,delta) ~ C(delta) / (1 - xi)`

as `xi -> 1^-`, with

`C(delta) = (9 delta + 11)^4 / [ 81 (9 delta^2 + 18 delta + 11)^2 ]`.

This is a theorem-level improvement over Stage 16:

- if `R_target < 1`, the target is impossible on the stable branch,
- if `R_target = 1`, the only hit is the unloaded onset point `xi = 0`,
- if `R_target > 1`, there is one and only one stable selected-branch softening depth `xi_req` that hits the target.

So the Stage-16 onset inequality is now upgraded to a full uniqueness theorem for the selected normalization locus.

---

## 4. Exact required total loading

The Stage-17 total loading law becomes, in D/N dimensionless form,

`alpha_req(xi,delta)`
`= 9 pi^2 A xi (xi + delta) / [ 8 (9 delta + 11 xi) ]`.

As `xi -> 1^-`, this tends to the exact stable-branch softening threshold

`alpha_crit = 9 pi^2 A (1 + delta) / [ 8 (11 + 9 delta) ]`,

which is the same refined threshold already found in Stage 16.

So once the unique `xi_req` solving `F(xi,delta)=R_target` is known, the unique required total directional loading follows immediately from the formula above.

---

## 5. Exact required support coupling

The physical total loading still splits as

`alpha_0 = g_B^2 / varpi^2 + alpha_mix`,

with

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 )`.

Therefore the exact support/BdG loading needed to realize the selected normalization locus is

`g_B,req^2 / varpi^2 = alpha_req(xi_req,delta) - alpha_mix`

or explicitly

`g_B,req^2 / varpi^2`
`= 9 pi^2 A xi_req (xi_req + delta) / [ 8 (9 delta + 11 xi_req) ]`
`  - Chi^2 / ( Omega_U^2 Delta_0 )`.

This gives one more exact feasibility gate:

`g_B,req^2 >= 0`.

So even after the unique normalization locus is found, the completed moving-throat PDE still has to place the physical support coupling on the nonnegative side of this equation.

---

## 6. Useful asymptotics

### 6.1 Near-onset regime

For `xi << 1`, the exact D/N shape function expands as

`F(xi,delta)`
`= 1 + ( 1 + 8/(9 delta) ) xi`
`  + ( 1 + 8/(9 delta) - 28/(27 delta^2) ) xi^2 + O(xi^3)`.

So if the target is only slightly above onset,

`R_target = 1 + eps_R`,

then the unique selected branch point is approximately

`xi_req ~= eps_R / ( 1 + 8/(9 delta) )`.

The required total loading then begins as

`alpha_req ~= pi^2 A xi_req / 8`.

### 6.2 Near-softening regime

For very large target ratio `R_target`, the unique solution lies close to softening and obeys

`1 - xi_req ~= C(delta) / R_target`,

with

`C(delta) = (9 delta + 11)^4 / [ 81 (9 delta^2 + 18 delta + 11)^2 ]`.

So large normalization demand pushes the selected branch into the thin neighborhood just below softening.

---

## 7. Best current theorem gate after Stage 18

The selected-branch quadrupole normalization bottleneck has now been reduced to the smallest exact reduced form reached so far.

The completed moving-throat PDE has to determine the microscopic quantities entering

`R_target = N_Q^(target) A / ( beta_0 kappa_0^2 )`,

and the bare anisotropy ratio

`delta = DeltaK_ax / A`.

Once those are known,

1. the unique stable normalization locus `xi_req` is determined by `F(xi,delta)=R_target`,
2. the unique required total loading is `alpha_req(xi_req,delta)`,
3. and the unique required support coupling is `g_B,req^2 / varpi^2 = alpha_req - alpha_mix`.

So the open problem is no longer a vague multi-parameter search.
It is now a one-dimensional universal D/N branch-shape problem plus one exact support-coupling feasibility check.
