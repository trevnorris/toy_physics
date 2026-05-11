# Moving-Throat PDE — Stage 034: Softening-Depth Normal Form, Exact Secular Reduction, and Elimination of the Selected-Mode Eigenvector Algebra

## Purpose

Stage 16 reduced the selected-branch 2.5PN normalization problem to one explicit microscopic equation,

`N_-(alpha_0) = N_Q^(target)`,

with `N_-` written in terms of the selected wall eigenvalue `lambda_-`, the selected overlap `s_-`, and the total directional loading `alpha_0`.

That was already sharp, but the formulation still carried unnecessary spectral baggage:
- the selected eigenvector `e_-`,
- the selected overlap `s_- = (v.e_-)^2`,
- and the branch stiffness `lambda_-`.

The next natural step is therefore to trade the selected eigenvalue for the more physical quantity

`x := A - lambda_-`,

namely the **softening depth** of the selected lower wall mode below the flat-branch stiffness `A`.

This removes the eigenvector algebra entirely.
In the first rank-1 selected-branch model, the full normalization problem becomes an explicit scalar problem in `x`.
The main outputs are:

1. an exact secular law `alpha_0(x)`,
2. an exact selected overlap `s_-(x)`,
3. an exact normalization product `N_-(x)`,
4. a manifestly monotone loading map `d alpha_0 / dx > 0`,
5. and an exact required support-loading formula `g_B^2/varpi^2` once `x` is known.

So Stage 17 replaces the spectral unknowns `(lambda_-, e_-)` by a single scalar deformation coordinate `x`.

---

## 1. Exact softening-depth variable

Recall the selected lower wall mode of the Stage-12/16 rank-1 reduced operator,

`D_wall = diag(A, A + DeltaK_ax) - alpha_0 v v^T`,

with

`v = (kappa_0, kappa_1)^T`.

On the stable branch,

`0 < lambda_- < A`,

so define the exact softening depth

`x := A - lambda_-`.

Then

`lambda_- = A - x`,

and stability is simply

`0 <= x < A`.

The selected branch starts at `x=0` when the directional loading vanishes and approaches softening as `x -> A^-`.

---

## 2. Exact secular equation in softening-depth form

The selected branch is the nontrivial rank-1-loaded solution of

`1 = alpha_0 [ kappa_0^2 / (A - lambda_-) + kappa_1^2 / (A + DeltaK_ax - lambda_-) ]`.

Substituting `lambda_- = A - x` gives the exact softening-depth secular equation

`1 = alpha_0 [ kappa_0^2 / x + kappa_1^2 / (x + DeltaK_ax) ]`.

Solving for the total directional loading gives

`alpha_0(x) = x (x + DeltaK_ax) / [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]`.

This is the first key simplification:

the selected-branch loading is now an explicit rational function of the softening depth.

---

## 3. Exact selected overlap in softening-depth form

The selected source/wall overlap may be written through the exact secular derivative identity,

`s_- = - d lambda_- / d alpha_0`.

Because `lambda_- = A - x`, this becomes

`s_- = d x / d alpha_0`.

Using the secular law above, the exact overlap is

`s_-(x) = [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]^2`
`         / [ kappa_0^2 (x + DeltaK_ax)^2 + kappa_1^2 x^2 ]`.

So the selected overlap is also a rational function of the same scalar variable `x`.
No explicit eigenvector is needed anymore.

---

## 4. Exact selected-branch normalization product in softening-depth form

Stage 15/16 gave the invariant selected-branch normalization quantity as

`N_-(alpha_0) = beta_0 s_-^2 / ( kappa_0^2 lambda_- )`.

Now substitute

`lambda_- = A - x`,

and the exact overlap formula from Section 3.
This gives

`N_-(x) = beta_0 [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]^4`
`         / { kappa_0^2 (A - x) [ kappa_0^2 (x + DeltaK_ax)^2 + kappa_1^2 x^2 ]^2 }`.

So the full selected-branch quadrupole normalization theorem is no longer an eigenvector problem.
It is one scalar equation in the softening depth:

`N_-(x) = N_Q^(target)`.

---

## 5. Exact monotonicity of the loading map

Differentiate the secular law.
One obtains

`d alpha_0 / dx`
`= [ kappa_0^2 (x + DeltaK_ax)^2 + kappa_1^2 x^2 ]`
`  / [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]^2`
`> 0`.

So the total directional loading grows strictly monotonically with softening depth.
This means:

- every stable selected-branch loading corresponds to exactly one softening depth,
- and vice versa.

This is useful because it turns the Stage-14 monotonicity statement into a direct scalar branch parameterization.

---

## 6. Exact required support loading once the softening depth is known

Stage 16 already split the total directional loading into

`alpha_0 = g_B^2 / varpi^2 + alpha_mix`,

with

`alpha_mix = Chi^2 / ( Omega_U^2 Delta_0 )`,

`Chi = Omega_U^2 g_W + g_R g_U`,

`Delta_0 = Omega_U^2 Omega_W^2 - g_R^2 sigma`.

So once the required softening depth `x` is known, the exact total loading is `alpha_0(x)`, and the exact required support contribution is

`g_B,req^2 / varpi^2 = alpha_0(x) - alpha_mix`

or explicitly

`g_B,req^2 / varpi^2`
`= x (x + DeltaK_ax) / [ kappa_0^2 (x + DeltaK_ax) + kappa_1^2 x ]`
`  - Chi^2 / ( Omega_U^2 Delta_0 )`.

So the support/BdG loading needed to hit the universal 2.5PN target is now an explicit function of one scalar deformation coordinate.

---

## 7. Best current theorem gate after Stage 17

The selected-branch theorem bottleneck is now even narrower than at Stage 16.
Instead of solving for a loaded eigenvector and its overlap, the reduced program only needs to determine the physical softening depth `x` of the selected lower wall mode.
Once that is known,

- the total directional loading is fixed by `alpha_0(x)`,
- the selected source map is fixed by `s_-(x)`,
- the normalization product is fixed by `N_-(x)`,
- and the required support loading follows directly.

So the remaining gap is no longer hidden in eigenvector algebra.
It is a scalar branch-placement problem.
