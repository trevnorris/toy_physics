# Moving-Throat PDE — Stage 062: Parent-Action Projection of the Microscopic Support/Source Gain

## Purpose

Stage 061 reduced the support/source theorem gap to one microscopic gain,

`G_micro = chi_sigma Lambda_phi^2 / K_X,`

but that quantity was still written in the language of the reduced axial source/support operator.

The next honest step is therefore to derive that gain from the **parent 4D action** rather than leave it as a phenomenological effective constant.

This stage does that by projecting the linearized GNLS/confinement sector of the parent theory onto one compressional source channel and one throat-support channel.

The main results are:

1. starting from the parent matter energy

   `H_psi = int d^4X [ (hbar^2/2m)|D_i psi|^2 + V_conf rho + U(rho) ],`

   with the frozen `n=5` EOS

   `U(rho) = K rho^5 / 4,`

   `h(rho) = dU/drho = (5K/4) rho^4,`

   the exact local compressional stiffness is

   `h'(rho_*) = 5 K rho_*^3 = m c_(s*)^2 / rho_*;`
2. after linearizing about a static throat branch `rho = rho_* + delta rho`, the compressional part of the matter energy is

   `(1/2) h'(rho_*) (delta rho)^2;`
3. if the support field enters the confinement as the first linear throat-shape perturbation

   `delta V_conf(s,y) = - g_phi chi_phi(y) phi(s),`

   and the source density is truncated to one transverse compression channel

   `delta rho(s,y) = sigma(s) chi_sigma(y),`

   then the parent 4D energy reduces exactly to the Stage-060 form

   `F_red[sigma,phi]`
   `= int_0^L ds [ (Theta_sigma/2) sigma^2 - Lambda_phi sigma phi`
   `               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2,`

   with

   `Theta_sigma = h'(rho_*) N_(sigma sigma),`

   `Lambda_phi = g_phi O_(sigma phi);`
4. here the parent overlap invariants are

   `N_(sigma sigma) = int d^3y chi_sigma(y)^2,`

   `N_(phi phi)     = int d^3y chi_phi(y)^2,`

   `O_(sigma phi)   = int d^3y chi_sigma(y) chi_phi(y);`
5. the effective source susceptibility is therefore

   `chi_sigma^(eff) = 1 / Theta_sigma`
   `               = rho_* / [ m c_(s*)^2 N_(sigma sigma) ];`
6. and the microscopic gain becomes the explicit parent-action quantity

   `G_micro`
   `= chi_sigma^(eff) Lambda_phi^2 / K_X`
   `= rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ];`
7. introducing the source/support coherence factor

   `C_(sigma phi)^2 := O_(sigma phi)^2 / [ N_(sigma sigma) N_(phi phi) ],`

   one obtains the exact factorization

   `G_micro = [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2,`

   with `0 <= C_(sigma phi)^2 <= 1` by Cauchy–Schwarz;
8. finally, inserting `kappa = K_X L^2 / T_X` gives the microscopic fixed-point strength

   `Xi_micro = kappa G_micro`
   `        = rho_* g_phi^2 O_(sigma phi)^2 L^2 / [ m c_(s*)^2 T_X N_(sigma sigma) ].`

So the Stage-061 gain is no longer a free microscopic placeholder. It is the compressional susceptibility of the `n=5` GNLS medium times the square of one parent confinement/support overlap, divided by the baseline support stiffness.

---

## 1. Parent 4D matter energy and the `n=5` compressional stiffness

The parent 4D theory already fixes the matter sector to the gauged GNLS form with confinement and the frozen stiff EOS,

`H_psi = int d^4X [ (hbar^2/2m)|D_i psi|^2 + V_conf(X;a,L) rho + U(rho) ],`

`U(rho) = K rho^5 / 4,`

`h(rho) = dU/drho = (5K/4) rho^4,`

`c_s^2(rho) = (1/m) dP/drho = (5K/m) rho^4.`

Therefore

`h'(rho) = 5K rho^3 = m c_s^2(rho) / rho.`

Linearizing about a static throat branch `rho = rho_* + delta rho` and subtracting the equilibrium linear term leaves the local compressional quadratic energy

`delta H_comp = (1/2) int d^4X h'(rho_*) (delta rho)^2.`

So the parent matter sector already provides the exact local scalar stiffness needed by the Stage-060 source entropy term.

---

## 2. Parent-action reduction to one source channel and one support channel

Let the support field be the first axial throat-support amplitude `phi(s)` and let its leading linear effect on the confinement be

`delta V_conf(s,y) = - g_phi chi_phi(y) phi(s),`

where:

- `s in [0,L]` is the throat axis,
- `y` denotes the transverse/cross-sectional coordinates,
- `chi_phi(y)` is the transverse support profile,
- `g_phi` is the parent confinement-loading amplitude.

Now truncate the density perturbation to one transverse compression channel,

`delta rho(s,y) = sigma(s) chi_sigma(y),`

with source profile `chi_sigma(y)`.

Insert this into the parent linearized energy. The compressional term becomes

`(1/2) h'(rho_*) int_0^L ds sigma(s)^2 int d^3y chi_sigma(y)^2,`

while the linear support/source coupling becomes

`- int_0^L ds sigma(s) phi(s) g_phi int d^3y chi_sigma(y) chi_phi(y).`

Define the overlap invariants

`N_(sigma sigma) = int d^3y chi_sigma^2,`

`N_(phi phi)     = int d^3y chi_phi^2,`

`O_(sigma phi)   = int d^3y chi_sigma chi_phi.`

Then the parent energy reduces exactly to

`F_red[sigma,phi]`
`= int_0^L ds [ (Theta_sigma/2) sigma^2 - Lambda_phi sigma phi`
`               + (T_X/2) phi_s^2 + (K_X/2) phi^2 ] + (K_m/2) phi(0)^2,`

with

`Theta_sigma = h'(rho_*) N_(sigma sigma),`

`Lambda_phi = g_phi O_(sigma phi).`

This is precisely the Stage-060 reduced support/source free energy, now derived as a one-channel projection of the parent 4D action.

---

## 3. Exact parent formula for the effective source susceptibility

The reduced source susceptibility is

`chi_sigma^(eff) = 1 / Theta_sigma`
`               = 1 / [ h'(rho_*) N_(sigma sigma) ].`

Using the exact `n=5` identity above,

`h'(rho_*) = m c_(s*)^2 / rho_*,`

one gets

`chi_sigma^(eff) = rho_* / [ m c_(s*)^2 N_(sigma sigma) ].`

So the source compliance is no longer a free “entropy constant.” It is fixed by the local GNLS compressibility of the parent medium, dressed only by the chosen transverse source normalization.

---

## 4. Exact parent formula for the microscopic gain

Stage 061 defined the microscopic gain by

`G_micro = chi_sigma Lambda_phi^2 / K_X.`

Insert the projected parent coefficients:

`G_micro`
`= [ 1 / ( h'(rho_*) N_(sigma sigma) ) ] [ g_phi^2 O_(sigma phi)^2 ] / K_X`
`= rho_* g_phi^2 O_(sigma phi)^2 / [ m c_(s*)^2 K_X N_(sigma sigma) ].`

This is the first explicit parent-action formula for the gain that controls the operator phase diagram.

It says that the support/source drive is large when:

- the local medium is easily compressible (`rho_* / m c_(s*)^2` large),
- the confinement/support loading amplitude `g_phi` is large,
- the source and support transverse channels overlap strongly,
- and the baseline support stiffness `K_X` is small.

---

## 5. Coherence factor and the exact overlap decomposition

Define the source/support coherence factor

`C_(sigma phi)^2 := O_(sigma phi)^2 / [ N_(sigma sigma) N_(phi phi) ].`

Then

`0 <= C_(sigma phi)^2 <= 1`

by Cauchy–Schwarz, and the microscopic gain factorizes as

`G_micro`
`= [ rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X) ] C_(sigma phi)^2.`

This is useful because it separates the microscopic problem into two independent pieces:

1. a **strength scale**

   `rho_* g_phi^2 N_(phi phi) / (m c_(s*)^2 K_X),`

2. a **profile-coherence factor**

   `C_(sigma phi)^2 in [0,1].`

So the unresolved PDE-side question is no longer a single opaque parameter. It is the combination of a parent confinement-loading strength and the coherence with which the source and support channels line up on the true branch.

---

## 6. Microscopic fixed-point strength `Xi_micro`

Stage 61 used

`Xi_micro = kappa G_micro,`

with

`kappa = K_X L^2 / T_X.`

Therefore the parent projected formula becomes

`Xi_micro`
`= rho_* g_phi^2 O_(sigma phi)^2 L^2 / [ m c_(s*)^2 T_X N_(sigma sigma) ].`

So the full operator-selected branch point

`Pe = Xi_micro Delta(Pe;kappa,eta)`

is now determined by three parent-action ingredients only:

- the local compressional stiffness of the `n=5` GNLS medium,
- the confinement/support loading overlap,
- and the support tension/Robin geometry sector entering `(T_X,kappa,eta)`.

That is the cleanest microscopic restatement yet of the support/source theorem gap.
