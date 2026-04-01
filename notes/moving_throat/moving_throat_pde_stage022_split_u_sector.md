# Moving-Throat PDE — Stage 22: First Non-Flat `U` Doublet, Exact Split Continuum Map, and the Direction-Splitting Theorem

## Purpose

Stage 21 compressed the selected quadrupole placement problem to the dimensionless continuum ratios

`eps_eta, eps_W, rho, Z_W, delta_0, Lambda`,

but it still relied on one structural simplification:

> the internal `U` sector was taken in the flat-doublet limit, so the first two N/N `U` channels were degenerate.

That was the right minimal starting point, but it is not the first genuinely nontrivial continuum operator.

The next exact question is therefore:

> what survives, and what breaks, once the first axial structure of the internal `U` sector is turned on?

This stage answers that question exactly.

The main result is that the **scalar placement map survives**, but the **directional simplicity does not**.

More concretely:

1. the direct wall softening becomes mode-dependent,
2. the bare anisotropy ratio shifts to a new exact value `delta_split`,
3. the mixed-sector blocking ratio becomes an exact split quantity `eps_W_split`,
4. the mixed loading vector rotates away from the source/support direction,
5. and the old Stage-21 one-direction picture survives **iff** the new direction-splitting invariant vanishes.

So the flat-`U` assumption was not harmless bookkeeping. It was precisely what made source, support, and mixed loading all point in the same wall-basis direction.

---

## 1. Turning on the first axial `U` structure

Keep the wall, `phi`, and `W` sectors exactly as in Stage 20, but replace the flat internal `U` block by the first nontrivial N/N continuum operator

`L_U = (mu_U/2) dot(U)^2 - (T_U/2) (U')^2 - (K_U/2) U^2.`

Then the first two N/N internal modes have exact stiffnesses

`K_(U0) = K_U,`

`K_(U1) = K_U + pi^2 T_U / L^2 = K_U (1 + delta_U),`

where the new internal axial-splitting ratio is

`delta_U := pi^2 T_U / (L^2 K_U).`

So the flat-doublet Stage-20 limit is exactly the special case

`delta_U = 0.`

The wall basis itself is unchanged:

`u_0(s) = 1/sqrt(L),`

`u_1(s) = sqrt(2/L) cos(pi s/L),`

and the same D/N overlap data remain

`kappa_0 = 2 sqrt(2)/pi,`

`kappa_1 = -4/(3 pi),`

`sigma = kappa_0^2 + kappa_1^2 = 88/(9 pi^2),`

`lambda_0 := kappa_1^2 / kappa_0^2 = 2/9.`

---

## 2. Exact split direct softening and shifted anisotropy

The wall-U coupling is still diagonal in the N/N basis, but the denominators are no longer the same. The direct softening therefore becomes

`A_0 = [K_eta^(eff) - c_(etaU)^2 / K_U] / mu_eta,`

`A_1 = [K_eta^(eff)(1 + delta_0) - c_(etaU)^2 / K_(U1)] / mu_eta,`

with

`K_eta^(eff) := K_eta + 6 T_Omega,`

`delta_0 := pi^2 T_w / (L^2 K_eta^(eff)),`

`eps_eta := c_(etaU)^2 / (K_U K_eta^(eff)).`

After exact rearrangement,

`A_0 = K_eta^(eff) (1 - eps_eta) / mu_eta,`

`A_1 = A_0 (1 + delta_split),`

with the exact shifted anisotropy ratio

`delta_split = [ delta_0 + eps_eta delta_U / (1 + delta_U) ] / (1 - eps_eta).`

So internal axial structure raises the bare selected-branch anisotropy even if the wall sector itself is unchanged.

For small `delta_U`,

`delta_split = delta_0/(1 - eps_eta) + [eps_eta/(1 - eps_eta)] delta_U + O(delta_U^2).`

---

## 3. Exact split mixed blocking ratio

The `U/W` block also feels the split directly through the overlap-weighted inverse kernel

`S_U = kappa_0^2 / K_U + kappa_1^2 / K_(U1).`

Using the Stage-21 flat-doublet ratio

`eps_W := c_(UW)^2 sigma / (K_U K_W^(eff)),`

with

`K_W^(eff) := K_W + pi^2 T_W / (4 L^2),`

one finds the exact split blocking ratio

`eps_W_split = eps_W [ 1 - (2/11) delta_U / (1 + delta_U) ].`

So the first axial `U` structure lowers the mixed blocking ratio relative to the flat-doublet value.

For small `delta_U`,

`eps_W_split = eps_W [ 1 - (2/11) delta_U ] + O(delta_U^2).`

---

## 4. Exact mixed loading vector and the direction-splitting theorem

The decisive new object is the mixed loading vector.

Define the flat-mode interference ratio

`rho_0 := c_(UW) c_(etaU) / (K_U c_(etaW)).`

Then the two wall-basis components of the mixed loading vector are

`z_0 = kappa_0 g_W (1 + rho_0),`

`z_1 = kappa_1 g_W (1 + rho_0/(1 + delta_U)),`

where

`g_W = c_(etaW) / sqrt(mu_eta mu_W).`

Equivalently,

`z_1 / z_0 = (kappa_1 / kappa_0) R_U,`

with the exact direction-ratio factor

`R_U := [ 1 + rho_0/(1 + delta_U) ] / (1 + rho_0).`

This isolates the first real failure mode of the flat-`U` simplification.

In Stage 20–21 the mixed loading direction was proportional to the source/support vector `v = (kappa_0,kappa_1)^T`.
Once `delta_U != 0`, the mixed loading vector is instead proportional to

`z = ( z_0, z_1 )^T,`

which is generally **not collinear** with `v`.

The exact direction-splitting invariant is

`D_dir := kappa_0 z_1 - kappa_1 z_0`

`      = - kappa_0 kappa_1 g_W rho_0 delta_U / (1 + delta_U).`

So the exact collinearity theorem is

`D_dir = 0  <=>  delta_U = 0  or  rho_0 = 0,`

which means:

- flat internal `U` doublet, or
- zero `U/W` interference,

are the only ways to keep the old one-direction Stage-21 geometry exactly intact.

For small `delta_U`,

`R_U = 1 - [rho_0/(1 + rho_0)] delta_U + O(delta_U^2).`

So on the natural constructive branch `rho_0 > 0`, the first axial `U` structure pushes the mixed loading direction away from the flat-doublet direction in a controlled linear way.

---

## 5. Exact split continuum placement map

Even though the directions split, the scalar placement data themselves still factorize cleanly.

Using the Stage-21 dimensionless overlap ratio

`Z_W := c_(etaW)^2 / ( K_eta^(eff) K_W^(eff) ),`

and the same radiative demand scale

`Lambda := 27 pi^2 G c_s^5 K_W^(eff) / (20 a^5 c^5 mu_W),`

the exact split placement formulas are

`M_mix^(split U)`
`= 8 Z_W (1 + rho_0)^2`
`  / [ pi^2 (1 - eps_eta) (1 - eps_W_split) ],`

`R_target^(split U)`
`= Lambda (1 - eps_eta) (1 - eps_W_split)^2`
`  / [ Z_W (1 + rho_0)^2 ].`

And the exact product law survives:

`R_target^(split U) M_mix^(split U)`
`= 8 Lambda (1 - eps_W_split) / pi^2.`

So the Stage-21 factorization survives at the scalar-placement level.
What changes is the **directional geometry** seen by the selected branch.

For small `delta_U`,

`M_mix^(split U) = M_mix^(flat) [ 1 - 2 eps_W delta_U / (11(1 - eps_W)) ] + O(delta_U^2),`

`R_target^(split U) = R_target^(flat) [ 1 + 4 eps_W delta_U / (11(1 - eps_W)) ] + O(delta_U^2).`

So positive internal axial splitting lowers the mixed baseline and raises the normalization demand ratio before any selected-branch/source-map correction is even considered.

---

## 6. Best current theorem statement after Stage 22

The first non-flat `U` doublet does **not** destroy the continuum placement map.
It does something subtler and more important.

It separates two statements that were accidentally fused in Stages 20–21:

1. **scalar placement factorization**, which survives, and
2. **source/loading collinearity**, which does not survive generically.

That means the real structural role of the flat-`U` assumption is now clear.
It was not merely making the formulas shorter. It was enforcing the coincidence of

- the source map direction,
- the support direction,
- and the mixed loading direction.

Once the first axial `U` structure is turned on, the next theorem problem is no longer “what are the scalar placement ratios?”
Those are still exact.

The next theorem problem is:

> how does the selected-branch normalization law deform when the source vector and the loading vector are no longer the same?

That is the target of Stage 23.
