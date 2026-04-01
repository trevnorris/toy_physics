# Moving-Throat PDE — Stage 16: Microscopic Selected-Branch Normalization Equation, Exact Stability Gate, and the First Coupling-Level Onset Criterion

## Purpose

Stage 15 removed the abstract selected-branch source-map factor on the natural D/N source branch and rewrote the invariant quadrupole target as

`beta_0 s_-^2 / (kappa_0^2 lambda_-) = 54 G c_s^5 / (5 a^5 c^5)`.

That was a real narrowing, but it still left the theorem gate partly hidden inside the shorthand spectral data `(beta_0, alpha_0, Xi_0)`.

The next step is therefore to write the selected-branch normalization problem directly in the microscopic couplings of the first explicit finite-throat kernel model and to separate three logically different questions:

1. **existence of a stable selected branch**,
2. **entry into the normalization window**,
3. **exact hit of the universal target at the physical loading**.

This stage does that.

The main result is that the selected-branch 2.5PN test is now one exact microscopic equation together with two exact necessary inequalities.
The stable natural branch must satisfy

`Delta_0 > 0`,

`A > 0`,

`alpha_0 < alpha_crit`,

and its normalization must obey

`N_-(g_B,g_U,g_W,g_R,...) = 54 G c_s^5 / (5 a^5 c^5)`

with an explicit closed formula for `N_-`.

There is also a useful exact necessary onset condition:
if the branch starts above the universal target at zero directional loading, then the natural positive-loading branch can never come back down to hit it, because the selected normalization product is strictly monotone increasing.
So the first coupling-level onset test is already exact.

---

## 1. Microscopic coupling abbreviations

For the first explicit finite-throat kernel model, define

`A = K_0 - g_U^2 / Omega_U^2`,

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`,

`Chi = Omega_U^2 g_W + g_R g_U`.

Then the Stage-15 data are

`Xi_0 = g_U^2 / Omega_U^2`,

`beta_0 = Chi^2 / Delta_0^2`,

`alpha_0 = g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`.

So the physical branch is controlled by:

- one isotropic wall shift `Xi_0`,
- one directional loading strength `alpha_0`,
- one outgoing transfer factor `beta_0`,
- and the finite-throat overlap constants `kappa_0`, `kappa_1`, `sigma` already fixed exactly.

---

## 2. Exact microscopic selected-branch normalization equation

On the natural D/N source branch,

`mhat_-^2 = s_- / kappa_0^2`,

`P_{0,-} = beta_0 s_- / lambda_-`,

so the full invariant selected-branch product is

`N_-(alpha_0) := mhat_-^2 P_{0,-}`
`              = beta_0 s_-(alpha_0)^2 / ( kappa_0^2 lambda_-(alpha_0) )`.

The exact 2.5PN target is therefore

`N_-(alpha_0) = N_Q^(target)`

with

`N_Q^(target) = 54 G c_s^5 / (5 a^5 c^5)`.

Substituting the microscopic couplings gives

`Chi^2 / Delta_0^2 * s_-(alpha_0)^2 / ( kappa_0^2 lambda_-(alpha_0) )`
`= 54 G c_s^5 / (5 a^5 c^5)`

where

`alpha_0 = g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`.

So the selected-branch normalization theorem is now a single explicit equation in the first finite-throat kernel-model couplings.

---

## 3. Exact stability gate in coupling language

Before the normalization equation can even be asked, the selected lower wall mode has to exist on a stable branch.
That imposes three exact requirements.

### 3.1 Internal passive/conservative positivity

The conservative internal block must satisfy

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma > 0`.

If `Delta_0 <= 0`, the reduced `U/W` block is already singular or overmixed before the wall problem is even formed.

### 3.2 Positive flat-branch stiffness

The flat branch must start stable,

`A = K_0 - g_U^2 / Omega_U^2 > 0`.

So the isotropic internal zero-mode shift cannot already overwhelm the bare wall stiffness.

### 3.3 Selected-branch softening bound

The directional loading must remain below the exact refined threshold

`alpha_0 < alpha_crit`,

with

`alpha_crit = A(A + Delta K_ax)`
`             / [ (A + Delta K_ax) kappa_0^2 + A kappa_1^2 ]`.

Using the exact finite-throat constants,

`kappa_0^2 = 8 / pi^2`,

`kappa_1^2 = 16 / (9 pi^2)`,

this becomes

`alpha_crit = 9 pi^2 A(A + Delta K_ax) / [ 8(11 A + 9 Delta K_ax) ]`.

So the exact coupling-level stability gate is

`g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`
`< 9 pi^2 A(A + Delta K_ax) / [ 8(11 A + 9 Delta K_ax) ]`.

This is the first fully explicit stability window for the selected moving-throat quadrupole branch in the current reduced program.

---

## 4. Exact monotonicity and the branch-onset criterion

From Stage 14, the selected normalization product is strictly monotone increasing on the stable branch.
With the Stage-15 source-map reduction, that exact statement becomes

`d N_- / d alpha_0`
`= beta_0 [ 2 s_- (d s_- / d alpha_0) lambda_- + s_-^3 ]`
`  / ( kappa_0^2 lambda_-^2 )`
`> 0`

for every stable branch point.

So the zero-loading value is now an exact **necessary onset condition** for the universal target.
At `alpha_0 = 0`,

`s_-(0) = kappa_0^2`,

`lambda_-(0) = A`,

and therefore

`N_-(0) = beta_0 kappa_0^2 / A = Chi^2 kappa_0^2 / ( A Delta_0^2 )`.

Because `N_-` only increases with positive loading, a necessary condition for the physical natural branch to hit the universal target is

`N_-(0) <= N_Q^(target)`.

Equivalently,

`Chi^2 <= N_Q^(target) A Delta_0^2 / kappa_0^2`.

This may be read as an onset stiffness bound,

`K_0 >= g_U^2 / Omega_U^2 + kappa_0^2 Chi^2 / ( N_Q^(target) Delta_0^2 )`.

If this inequality fails, then the natural positive-loading branch starts **above** the universal target and can never come back down to hit it.
That does not yet solve the normalization equation, but it is an exact necessary branch-admissibility test.

---

## 5. Weak-loading expansion of the microscopic normalization product

For a stable branch with small directional loading,

`alpha_0 << alpha_crit`,

the exact selected normalization product has the expansion

`N_-(alpha_0)`
`= beta_0 kappa_0^2 / A`
`  + alpha_0 beta_0 kappa_0^2 ( 4 A kappa_1^2 + Delta K_ax kappa_0^2 ) / ( A^2 Delta K_ax )`
`  + O(alpha_0^2)`.

Using the exact finite-throat constants gives

`N_-(alpha_0)`
`= beta_0 * 8 / (pi^2 A)`
`  + alpha_0 * 64 beta_0 (8A + 9 Delta K_ax)`
`    / ( 9 pi^4 A^2 Delta K_ax )`
`  + O(alpha_0^2)`.

Now substitute the microscopic loading strength

`alpha_0 = g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 )`.

Then the first explicit weak-loading approximation to the physical selected-branch product is

`N_-^(phys)`
`= Chi^2 * 8 / ( pi^2 A Delta_0^2 )`
`  + [ g_B^2 / varpi^2 + Chi^2 / ( Omega_U^2 Delta_0 ) ]`
`    * 64 Chi^2 (8A + 9 Delta K_ax)`
`      / ( 9 pi^4 A^2 Delta K_ax Delta_0^2 )`
`  + O(alpha_0^2)`.

This is not yet the final theorem, but it is the first concrete approximation that lets us diagnose which microscopic lane is doing the main work:

- `g_U` lowers the wall through `A`,
- `g_R` and `g_W` feed the transfer factor through `Chi` and `Delta_0`,
- `g_B` pushes the branch along the monotone loading direction,
- and `Delta K_ax` controls how fast the selected source map can grow before softening.

---

## 6. Best current theorem gate after Stage 16

The selected-branch normalization bottleneck is now sharply microscopic.
The first explicit finite-throat kernel model has reduced the theorem gap to:

1. compute the physical coupling set
   `(g_B, g_U, g_W, g_R, varpi, Omega_U, Omega_W, K_0, Delta K_ax)`
   from the completed moving-throat PDE,
2. verify the exact stability window
   `Delta_0 > 0`, `A > 0`, `alpha_0 < alpha_crit`,
3. check the exact necessary onset inequality
   `N_-(0) <= N_Q^(target)`,
4. and then test the full microscopic normalization equation
   `N_-(alpha_0) = N_Q^(target)`.

So the open problem is no longer hidden in a free source map or a generic tuning story.
It is an explicit coupling-level spectral-placement problem on the selected stable quadrupole branch.
