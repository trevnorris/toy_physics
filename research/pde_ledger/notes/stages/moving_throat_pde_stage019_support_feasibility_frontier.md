# Moving-Throat PDE — Stage 19: Dimensionless Support-Feasibility Frontier for the Selected Quadrupole Branch

## Purpose

Stage 18 reduced the selected-branch normalization problem to the unique D/N locus

`F(xi,delta) = R_target`,

and gave the exact total loading required once the unique solution `xi_req` is found.
But there was still one extra feasibility step left outside the locus itself:

the support/BdG channel has to supply a **nonnegative** additional loading,

`g_B,req^2 >= 0`.

The natural next move is therefore to isolate the exact dimensionless function that measures how much of the total directional loading can be carried by the selected branch itself.
That is what this stage does.

The result is an exact second branch function,

`G(xi,delta) = 9 xi (xi + delta) / (9 delta + 11 xi)`,

such that

`g_B,req^2 / varpi^2 = (pi^2 A / 8) [ G(xi_req,delta) - M_mix ]`,

with

`M_mix = 8 alpha_mix / (pi^2 A)`
`      = 8 Chi^2 / (pi^2 A Omega_U^2 Delta_0)`.

So the selected quadrupole branch is support-feasible iff

`M_mix <= G(xi_req,delta)`.

This turns the final support check into an exact geometric inequality in the same dimensionless branch parameter `xi` that already controlled the normalization locus.

---

## 1. Exact support-feasibility function

From Stage 18 the required total loading is

`alpha_req(xi,delta) = 9 pi^2 A xi (xi + delta) / [ 8 (9 delta + 11 xi) ]`.

Split the total loading into the mixed-sector baseline and the support contribution,

`alpha_req = alpha_mix + g_B,req^2 / varpi^2`,

with

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 )`.

Now define the dimensionless mixed baseline

`M_mix := 8 alpha_mix / (pi^2 A)`
`      = 8 Chi^2 / (pi^2 A Omega_U^2 Delta_0 )`.

Then the exact branch-dependent support-feasibility function is

`G(xi,delta) := 8 alpha_req / (pi^2 A)`
`           = 9 xi (xi + delta) / (9 delta + 11 xi)`.

The required support loading becomes

`g_B,req^2 / varpi^2 = (pi^2 A / 8) [ G(xi_req,delta) - M_mix ]`.

So the support/BdG sector is feasible exactly when

`G(xi_req,delta) >= M_mix`.

---

## 2. Exact monotonicity and endpoint values

Differentiate `G`.
One finds

`dG/dxi = 9 [ 9 delta^2 + 18 delta xi + 11 xi^2 ] / (9 delta + 11 xi)^2 > 0`.

So `G` is strictly monotone increasing on the stable branch.
Its exact endpoint values are

`G(0,delta) = 0`,

`G_max(delta) := lim_{xi -> 1^-} G(xi,delta)`
`             = 9 (1 + delta) / (9 delta + 11)`.

This is useful because it turns the support-feasibility condition into a sharp branch window.
For fixed `delta`, the selected branch can support at most

`M_mix < G_max(delta)`.

That is of course equivalent to the refined stability bound `alpha_mix < alpha_crit`, but the present stage makes the same statement directly in the dimensionless selected-branch geometry.

---

## 3. The exact admissible region in `(R_target, M_mix)` space

Stages 18 and 19 together now give two exact branch functions driven by the same parameter `xi`:

`R_target = F(xi,delta)`,

`M_crit = G(xi,delta)`.

For fixed `delta`, the stable selected quadrupole branch therefore traces an exact parametric admissibility frontier in the `(R_target, M_mix)` plane:

`xi in [0,1)  ->  ( F(xi,delta), G(xi,delta) )`.

Because both `F` and `G` are strictly monotone increasing,

- the normalization target picks out a unique `xi_req`,
- and the support feasibility test is then simply whether the actual baseline `M_mix` lies below the critical value `G(xi_req,delta)`.

So the combined theorem gate is now:

1. `R_target >= 1` to enter the unique normalization locus,
2. `M_mix <= G(xi_req,delta)` to keep the required support loading nonnegative.

That is the first exact combined reachability-plus-feasibility statement for the selected moving-throat quadrupole branch.

---

## 4. Near-onset asymptotics

For `xi << 1`, the support-feasibility function expands as

`G(xi,delta) = xi - 2 xi^2 / (9 delta) + O(xi^3)`.

Combined with the Stage-18 onset relation

`xi_req ~= (R_target - 1) / (1 + 8/(9 delta))`,

this gives the first near-onset support-feasibility estimate,

`M_crit ~= (R_target - 1) / (1 + 8/(9 delta))`

up to the first nonlinear correction.

So just above onset, the admissible mixed baseline grows linearly with the excess normalization demand.

---

## 5. Best current theorem gate after Stage 19

The selected moving-throat quadrupole problem has now split cleanly into two exact scalar branch functions:

1. the **normalization function** `F(xi,delta)`,
2. the **support-feasibility function** `G(xi,delta)`.

For any fixed `delta`, the completed moving-throat PDE must determine a physical point `(R_target, M_mix)` such that

- `R_target` lands on the unique stable normalization locus,
- and `M_mix` lies below the exact support-feasibility frontier.

So the open theorem gap is now even sharper:

not “find the right eigenvector and source map,”
not even “find the right branch stiffness,”

but simply:

> does the completed moving-throat PDE place the physical defect on the exact admissible region of the universal `(F,G)` branch geometry?
