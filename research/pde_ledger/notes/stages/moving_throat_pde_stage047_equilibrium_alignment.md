# Moving-Throat PDE — Stage 47: Parent Equilibrium Source/Support Alignment and the Exact Matched-Layer Gain

## Purpose

Stages 45–46 reduced the support/source theorem gap to two parent quantities:

- the confinement-loading amplitude `g_phi`,
- and the source/support coherence `C_(sigma phi)^2`.

But both were still being treated as independent branch data.

The next honest step is to ask whether the **parent equilibrium equations themselves** tie the source channel to the support channel.

This stage shows that they do.

Within the same local compressional linearization already used in Stages 43–46, a quasi-static wall/support displacement does not excite an arbitrary source profile. It excites a very specific one:

`chi_sigma(y) = g_phi chi_phi(y) / H(y),`

where

`H(y) := h'(rho_*(y))`

is the local compressional stiffness of the static branch.

That gives four useful results immediately:

1. the source/support overlap invariants become

   `O_(sigma phi) = g_phi I_1,`

   `N_(sigma sigma) = g_phi^2 I_2,`

   with

   `I_1 = int d^3y chi_phi(y)^2 / H(y),`

   `I_2 = int d^3y chi_phi(y)^2 / H(y)^2;`
2. the coherence factor is no longer arbitrary:

   `C_(sigma phi)^2 = I_1^2 / [ N_(phi phi) I_2 ] <= 1;`
3. the exact eliminated-source support softening is

   `Delta K_X^(eq) = g_phi^2 I_1,`

   so the corresponding exact equilibrium gain is

   `G_eq = g_phi^2 I_1 / K_X;`
4. in the thin active layer where `H(y)` is approximately constant, the branch becomes exactly matched:

   `C_(sigma phi)^2 = 1,`

   `G_eq = g_phi^2 N_(phi phi) / [ K_X H_w ],`

   reproducing the Stage-45/46 best-alignment formulas.

So the parent equilibrium problem has already removed one large ambiguity: the support-induced source channel is not a free profile choice.

---

## 1. Local static response law from the parent equilibrium branch

Stage 45 projected the parent GNLS/confinement energy onto one source channel `sigma(s)` and one support channel `phi(s)` using the local compressional quadratic term

`(1/2) h'(rho_*) (delta rho)^2`

and the linear support loading

`delta V_conf = - g_phi chi_phi(y) phi(s).`

Inside the same local static linearization, the parent equilibrium response of the density perturbation is

`H(y) delta rho(s,y) + delta V_conf(s,y) = 0,`

with

`H(y) := h'(rho_*(y)).`

Therefore the support-induced density perturbation is

`delta rho(s,y) = phi(s) chi_sigma(y),`

with the exact aligned profile

`chi_sigma(y) = g_phi chi_phi(y) / H(y).`

So the source channel is fixed pointwise by the support loading and the local compressibility.

---

## 2. Exact overlap invariants on the equilibrium-aligned branch

Using the Stage-45 overlap definitions,

`N_(phi phi) = int d^3y chi_phi^2,`

`O_(sigma phi) = int d^3y chi_sigma chi_phi,`

`N_(sigma sigma) = int d^3y chi_sigma^2,`

the equilibrium-aligned profile gives

`O_(sigma phi) = g_phi I_1,`

`N_(sigma sigma) = g_phi^2 I_2,`

with

`I_1 := int d^3y chi_phi(y)^2 / H(y),`

`I_2 := int d^3y chi_phi(y)^2 / H(y)^2.`

Therefore the coherence factor becomes

`C_(sigma phi)^2 = O_(sigma phi)^2 / [ N_(phi phi) N_(sigma sigma) ]`
`                = I_1^2 / [ N_(phi phi) I_2 ].`

This is already a strong theorem statement. The coherence is no longer an unconstrained branch datum. It is fixed by how much the compressional stiffness varies across the active support layer.

---

## 3. Exact coherence bound and the matched-layer limit

By Cauchy–Schwarz,

`I_1^2 <= N_(phi phi) I_2,`

so

`0 <= C_(sigma phi)^2 <= 1.`

The equality condition is also exact: `C_(sigma phi)^2 = 1` if and only if `1/H(y)` is constant on the support of `chi_phi(y)`.

So in the physically important thin-wall / matched-layer regime where the support lives in a narrow active layer and the local compressional stiffness is nearly constant there,

`H(y) ~ H_w,`

one gets

`I_1 = N_(phi phi) / H_w,`

`I_2 = N_(phi phi) / H_w^2,`

and therefore

`C_(sigma phi)^2 = 1.`

This is the first exact parent reason for expecting the near-matched branch rather than a strongly misaligned one.

---

## 4. Exact eliminated-source support softening

The reduced static source/support energy has the form

`F[sigma,phi] = (1/2) Theta sigma^2 - Lambda sigma phi + (1/2) K_X phi^2.`

Eliminating the static source amplitude gives

`sigma_stat = Lambda phi / Theta,`

so the effective support energy is

`F_eff(phi) = (1/2) [ K_X - Lambda^2 / Theta ] phi^2.`

Therefore the support softening is exactly

`Delta K_X = Lambda^2 / Theta.`

On the equilibrium-aligned branch, the direct parent elimination gives

`Delta K_X^(eq) = g_phi^2 I_1.`

So the exact equilibrium gain is

`G_eq = Delta K_X^(eq) / K_X = g_phi^2 I_1 / K_X.`

This is slightly stronger than the Stage-45/46 formula because it no longer requires the source/support branch data to be treated as independent objects.

---

## 5. Constant-compressibility reduction and contact with Stages 45–46

In the matched-layer limit `H(y) ~ H_w`,

`I_1 = N_(phi phi) / H_w,`

so

`G_eq = g_phi^2 N_(phi phi) / [ K_X H_w ].`

Using the Stage-45 `n=5` identity

`H_w = h'(rho_w) = m c_(s,w)^2 / rho_w,`

this becomes

`G_eq = rho_w g_phi^2 N_(phi phi) / [ m c_(s,w)^2 K_X ].`

That is exactly the Stage-45/46 best-alignment gain with

`C_(sigma phi)^2 = 1.`

So the earlier best-case branch is not arbitrary. It is the natural thin-layer limit of the parent equilibrium-aligned source/support branch.

---

## 6. What Stage 47 changes

Before this stage, the unresolved microscopic support/source data were:

- one parent loading amplitude `g_phi`,
- and one independent coherence factor `C_(sigma phi)^2`.

After this stage, the situation is sharper.

The parent equilibrium equations already tie the source profile to the support profile. The remaining branch data are now:

1. the confinement-loading amplitude `g_phi`,
2. the support profile `chi_phi(y)`,
3. the local compressional stiffness profile `H(y)` across the active layer.

The coherence factor is no longer free. It is a derived quantity,

`C_(sigma phi)^2 = I_1^2 / [ N_(phi phi) I_2 ],`

and it is automatically near 1 when the active wall layer is thin enough that `H(y)` is nearly constant there.

That is the point where it becomes worthwhile to stop speaking abstractly about `g_phi` and evaluate it on a concrete parent confinement branch.
