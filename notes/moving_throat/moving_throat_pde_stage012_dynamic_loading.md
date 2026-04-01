# Moving-Throat PDE — Stage 12: Exact Dynamic Loading from the Coupled Wall/BdG/Maxwell/Mixed Operator, and the First Selected-Mode Quadrupole Damping Coefficient

## Purpose

Stage 11 deliberately introduced a reduced loading parameter `alpha` to turn the free profile angle `theta` into a real eigenproblem. That was the right move at the time, but it also left the central spectral question unfinished:

> what does `alpha` actually come from in the coupled wall/BdG/Maxwell/mixed operator, and how does the same operator feed the outgoing `l=2` odd term into the *selected* wall quadrupole mode?

This stage answers that question for the first nontrivial reduced operator that is still microscopically faithful to the project ontology.

The main result is that the coupled operator is much cleaner than it looked before.
Using the first two N/N wall modes, the lowest D/N support/mixed half-wave, and a brane-like internal zero-mode doublet, the full reduced wall self-energy splits **exactly** into two pieces:

- an isotropic shift `Xi(omega) I_2`,
- and a rank-1 directional load `alpha(omega) v v^T`.

So Stage 11 was not just a plausible toy model. Its rank-1 loading structure is the exact Schur complement of the first coupled wall/BdG/Maxwell/mixed block.

Even better, once the mixed channel is dressed by the passive/outgoing port, the same exact decomposition survives and the first odd coefficient of the selected wall mode becomes explicit.
So the bottleneck narrows again:

- the profile-selection problem is no longer carried by a free parameter,
- the conservative loading is now an exact rational function of the reduced couplings,
- and the outgoing `i omega^5` coefficient of the selected wall mode is now a one-line projection formula.

---

## 1. Minimal coupled axial operator in the `{u_0,u_1}` wall basis

Work in the first two exact N/N wall modes,

`u_0(s) = 1 / sqrt(L)`,

`u_1(s) = sqrt(2/L) cos(pi s / L)`,

and write the wall quadrupole shape as the coefficient vector

`q = (q_0, q_1)^T`.

The D/N overlap vector from Stages 10–11 is

`v = (kappa_0, kappa_1)^T`,

with

`kappa_0 = 2 sqrt(2) / pi`,

`kappa_1 = -4 / (3 pi)`,

`sigma = v.v = kappa_0^2 + kappa_1^2 = 88 / (9 pi^2)`.

So `sigma` is exactly the Stage-10 max-coupling value `|kappa|_max^2`.

Let the bare wall operator in this basis be

`D_eta^(bare)(omega) = diag( K_0 - M_0 omega^2, K_1 - M_1 omega^2 )`,

with

`K_0 = K_eta + 6 T_Omega`,

`K_1 = K_0 + Delta K_ax`,

`Delta K_ax = T_w pi^2 / L^2`.

Now keep the smallest internal block that still matches the earlier reduced hierarchy:

- one BdG support mode `phi` on the D/N half-wave,
- one brane-like internal doublet `u = (u_0^{int}, u_1^{int})^T` in the same N/N basis,
- one mixed `A_w/F_(mu w)/J^w` mode `W` on the same D/N half-wave.

Use the reduced frequency-domain kernels

`A_phi(omega) = varpi^2 - omega^2`,

`A_U(omega)   = Omega_U^2 - omega^2`,

`A_W(omega)   = Omega_W^2 - omega^2 - Pi_out(omega)`.

The exact reduced couplings are taken to be

- wall–BdG: `lambda_B (v.q) phi`,
- wall–U:   `lambda_U q.u`,
- wall–W:   `lambda_W (v.q) W`,
- U–W:      `- lambda_R (v.u) W`.

The last sign is the one compatible with the earlier conservative self-energy convention in which the mixed-sector cross term enters with the *positive* `+ 2 G_U G_W R` numerator after elimination.

This reduced operator is the minimal dynamic refinement of Stage 11:

- `q` is the wall profile we want to select,
- `phi` is the BdG support lane already turned on at Stage 3,
- `u` is the lowest brane-like internal gauge lane,
- `W` is the first mixed `A_w/F_(mu w)/J^w` lane,
- and `Pi_out` is the passive/outgoing dressing that can carry the `i omega^5` branch.

---

## 2. Exact Schur-complement decomposition of the wall self-energy

Eliminating the internal fields `(u, W, phi)` gives an exact effective wall operator of the form

`D_wall(omega) = D_eta^(bare)(omega) - Sigma_wall(omega)`.

The first key theorem of this stage is that the wall self-energy splits **exactly** as

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`,

where

`Xi(omega) = lambda_U^2 / A_U(omega)`

and

`alpha(omega) = lambda_B^2 / A_phi(omega)`
`             + ( A_U(omega) lambda_W + lambda_R lambda_U )^2`
`               / [ A_U(omega) Delta_UW(omega) ]`.

Here

`Delta_UW(omega) = A_U(omega) A_W(omega) - lambda_R^2 sigma`.

So the full coupled wall/BdG/Maxwell/mixed operator really does reduce to exactly the Stage-11 geometry:

`D_wall(omega) = [ D_eta^(bare)(omega) - Xi(omega) I_2 ] - alpha(omega) v v^T`.

This is a stronger result than the old Stage-11 ansatz because it identifies the two physically distinct effects:

1. `Xi(omega)` is an **isotropic internal zero-mode shift** coming from the brane-like `U` doublet,
2. `alpha(omega)` is the **directional loading strength** along the D/N overlap vector `v`.

So the first loaded profile-selection problem is no longer phenomenological.
It is the exact first Schur complement of the reduced wall/BdG/Maxwell/mixed operator.

---

## 3. Conservative static loading and the refined Stage-11 profile-selection law

On the conservative static branch set

`Pi_out(0) = 0`.

Then

`Xi_0 = lambda_U^2 / Omega_U^2`,

`Delta_0 = Omega_U^2 Omega_W^2 - lambda_R^2 sigma`,

`alpha_0 = lambda_B^2 / varpi^2`
`         + ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / ( Omega_U^2 Delta_0 )`.

So the conservative loaded wall matrix is

`K_eff^(0) = [[K_0 - Xi_0, 0], [0, K_1 - Xi_0]] - alpha_0 v v^T`.

The isotropic shift moves both wall levels equally,
so the profile-angle equation keeps the same Stage-11 structure with

`K_0 -> K_0 - Xi_0`,

`K_1 -> K_1 - Xi_0`,

`alpha -> alpha_0`.

In particular, the exact conservative angle law remains

`tan(2 theta_-)`
`= 2 alpha_0 kappa_0 kappa_1 / ( Delta K_ax + alpha_0 (kappa_0^2 - kappa_1^2) )`.

So the same sign theorem survives:
for positive `alpha_0`, the selected lower eigenvector still rotates in the negative `u_1` direction, i.e. toward the Stage-10 max-coupling branch and away from the blind-angle branch.

### 3.1 Refined conservative softening threshold

Because the internal `U` doublet contributes the isotropic shift `Xi_0`, the exact softening threshold is refined to

`alpha_crit`
`= (K_0 - Xi_0)(K_1 - Xi_0)`
`  / [ (K_1 - Xi_0) kappa_0^2 + (K_0 - Xi_0) kappa_1^2 ]`.

So the internal zero-mode sector makes quadrupole softening easier in a very concrete sense:

- it does **not** change the direction-selection formula,
- but it lowers the stability margin by shifting both bare wall levels downward before the directional rank-1 loading is applied.

This is the first exact place where the BdG/Maxwell/mixed operator changes the Stage-11 theorem content rather than merely justifying it.

---

## 4. Exact outgoing dressing of the dynamic loading strength

Now restore the passive/outgoing dressing in the mixed channel,

`A_W(omega) = Omega_W^2 - omega^2 - Pi_out(omega)`.

The exact decomposition survives unchanged:

`Sigma_wall(omega) = Xi(omega) I_2 + alpha(omega) v v^T`,

with the same closed formula for `alpha(omega)` but with the dressed `Delta_UW(omega)`.

Expanding to first order in the outgoing port gives

`alpha(omega) = alpha_cons(omega) + beta(omega) Pi_out(omega) + O(Pi_out^2)`

with the exact transfer factor

`beta(omega) = [ A_U(omega) lambda_W + lambda_R lambda_U ]^2 / Delta_cons(omega)^2`,

where

`Delta_cons(omega) = A_U(omega) (Omega_W^2 - omega^2) - lambda_R^2 sigma`.

So the directional loading strength inherits the outgoing branch with the same kind of positive transfer factor already seen in the one-lane Stage-4 mixed-sector analysis.

At zero frequency,

`beta_0 = ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta_0^2 >= 0`.

This is the second key theorem of the stage:

> once the passive/outgoing port is attached to the mixed channel, the *same* wall-profile loading parameter that selects the quadrupole eigenmode also inherits the odd `i omega^5` branch with a manifestly nonnegative conservative transfer factor.

---

## 5. First selected-mode odd quadrupole coefficient

On the natural compact outgoing `l=2` branch,

`Pi_out(omega) = + i Gamma_2^(port) omega^5 + O(omega^7)`

with

`Gamma_2^(port) = a^5 / (27 c_s^5)`

from the earlier outgoing Hankel/D/N audit.

Therefore the directional loading strength has the low-frequency form

`alpha(omega) = alpha_even(omega) + i beta_5 omega^5 + O(omega^7)`

with

`beta_5 = Gamma_2^(port) beta_0`
`      = Gamma_2^(port) ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta_0^2`.

Equivalently,

`beta_5 = [ a^5 / (27 c_s^5) ]`
`        * ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / ( Omega_U^2 Omega_W^2 - lambda_R^2 sigma )^2`.

Now project this odd piece onto the *conservative* lower eigenvector `e_-` of `K_eff^(0)`.
Since the odd operator is rank-1,

`delta D_-^(odd)(omega)`
`= - i beta_5 ( v.e_- )^2 omega^5 + O(omega^7)`

in the wall-operator convention used throughout these notes.

So the selected-mode damping strength is controlled by exactly two ingredients:

1. the dynamic transfer coefficient `beta_5`,
2. the selected profile overlap factor `(v.e_-)^2`.

The second factor is just the Stage-10 overlap evaluated on the Stage-11 selected eigenmode. In angle language,

`(v.e_-)^2 = kappa(theta_-)^2`.

So the operator-selected odd coefficient is

`delta D_-^(odd)(omega)`
`= - i beta_5 kappa(theta_-)^2 omega^5 + O(omega^7)`.

This is the first exact formula in the moving-throat PDE program that feeds the passive/outgoing `l=2` odd term directly into the **selected** wall quadrupole mode rather than into a hand-chosen profile.

---

## 6. Relation to the earlier Stage-4/5 transfer factors

The conservative combination

`P_0 = Omega_U^2 lambda_W + lambda_R lambda_U`

is the same `P`-type numerator that already appeared in the earlier isotropic normalization work.
So the Stage-12 directional-loading transfer factor

`beta_0 = P_0^2 / Delta_0^2`

is not a new unrelated object.
It is the same mixed-sector transfer numerator that the Stage-4/5 bridge had already isolated, now appearing on the first real wall-profile selection operator.

That is an important unification point:

- Stage 4 first isolated the outgoing mixed-sector transfer factor on one lane,
- Stage 8 reduced the normalization problem to the radial/axial quantities `P`, `Delta`, `Q`, `X`,
- Stage 11 turned the wall profile into a real eigenproblem,
- and Stage 12 shows that the same mixed-sector numerator controls the directional loading and the selected-mode odd term.

So the different branches of the roadmap are starting to meet.

---

## 7. Best current summary after Stage 12

The bottleneck has narrowed again.

The coupled wall/BdG/Maxwell/mixed operator now tells us that the first loaded profile-selection problem is governed by two exact reduced objects:

- an isotropic shift
  `Xi_0 = lambda_U^2 / Omega_U^2`,
- and a directional load
  `alpha_0 = lambda_B^2 / varpi^2 + ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / ( Omega_U^2 Delta_0 )`.

The outgoing odd term of the selected lower quadrupole mode is then

`delta D_-^(odd)(omega)`
`= - i [ a^5 / (27 c_s^5) ]`
`    * [ ( Omega_U^2 lambda_W + lambda_R lambda_U )^2 / Delta_0^2 ]`
`    * kappa(theta_-)^2 omega^5 + O(omega^7)`.

So the Stage-11 “free parameter” bottleneck is gone.
What remains is now much sharper:

1. compute the conservative branch data `(Xi_0, alpha_0)` on the true moving-throat branch,
2. check that `alpha_0 < alpha_crit` so the selected quadrupole wall mode stays on the stable side,
3. compute the selected overlap `kappa(theta_-)^2`,
4. and then test whether the resulting operator-selected odd coefficient lands on the required normalization branch.

That is a much smaller theorem gap than we had before.
It is no longer “derive some PDE loading somehow.”
It is a specific spectral-transfer problem.
