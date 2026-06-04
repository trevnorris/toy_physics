# Moving-Throat PDE — Stage 070: Explicit GNLS Wall-Shell Reduction for the First Support Branch

## Purpose

The reduced support/source program ended with the wall figure of merit

`W_wall = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell),`

but the actual branch data `(T_X, J_1, kappa, eta)` were still being treated as external inputs.

The next honest step is to derive those quantities from the parent GNLS action on the first explicit wall-support branch.

This stage does that for a thin active shell around the throat wall.

The main result is that, once the support field is taken to be a linearized density mode living on the wall shell,

- the axial support tension is fixed by the GNLS gradient term,
- the axial support stiffness is fixed by compressibility plus transverse curvature,
- the reduced wall figure of merit simplifies enormously,
- and on the matched thin-wall branch it is exactly the same object as the Stage-058/059 fixed-point coupling.

So the explicit moving-throat branch is now much more concrete.

---

## 1. Parent quadratic wall-shell support energy

Start from the parent GNLS matter energy density near a static wall background `rho_w`:

`H_psi = (hbar^2 / 2m) |grad sqrt(rho)|^2 + U(rho) + V_conf rho.`

Linearize around the wall layer with a real density perturbation `delta rho`. To quadratic order in the compressive support sector, the standard GNLS expansion gives

`E^(2)[delta rho]`
`= 1/2 int ds d^3y [ (hbar^2 / (4 m rho_w)) |grad delta rho|^2 + H_w (delta rho)^2 ],`

where

`H_w := h'(rho_w) = m c_(s,w)^2 / rho_w.`

Now project the perturbation onto a separated wall-shell mode,

`delta rho(s,y) = q(s) chi_phi(y),`

with axial coordinate `s in [0,L]` and transverse wall-shell coordinates `y`.

Then the quadratic support energy becomes

`E^(2)[q] = 1/2 int_0^L ds [ T_X q'(s)^2 + K_X q(s)^2 ],`

with

`T_X = (hbar^2 / (4 m rho_w)) N_(phi phi),`

`K_X = H_w N_(phi phi) + (hbar^2 / (4 m rho_w)) G_(phi phi),`

where

`N_(phi phi) := int d^3y chi_phi(y)^2,`

`G_(phi phi) := int d^3y |grad_y chi_phi(y)|^2.`

So the support tension and stiffness are no longer abstract branch data; they are explicit shell integrals of the wall-support profile.

---

## 2. Thin-shell wall profile

Take the explicit wall family already introduced in Stage 065,

`V_conf(r;a) = V0 f((r-a)/ell),`

with shell coordinate

`xi := (r-a)/ell.`

The induced support profile is naturally

`chi_phi(y) = f'(xi),`

and the loading amplitude is

`g_phi = V0 / ell.`

In the thin-shell approximation, the three-dimensional shell measure is

`d^3y = 4 pi r^2 dr ~ 4 pi a^2 ell dxi.`

Define the dimensionless wall-shape moments

`I_f := int dxi f'(xi)^2,`

`I_g := int dxi f''(xi)^2.`

Then

`N_(phi phi) = 4 pi a^2 ell I_f,`

`G_(phi phi) = 4 pi a^2 I_g / ell.`

Substituting these into the parent quadratic support coefficients gives the exact thin-shell branch formulas

`T_X = pi a^2 ell I_f hbar^2 / (m rho_w),`

`K_X = 4 pi a^2 ell I_f H_w + pi a^2 I_g hbar^2 / (m rho_w ell).`

Equivalently,

`kappa := K_X L^2 / T_X`
`      = 4 (m c_(s,w) L / hbar)^2 + (I_g / I_f) (L / ell)^2.`

This is the first explicit parent formula for the geometry/support parameter `kappa`.

---

## 3. Wall overlap moment and exact collapse of the wall figure of merit

Stage 065 already showed that, for an almost constant-compressibility active layer,

`J_1 = I_f / H_w.`

Insert this together with the explicit `T_X` above into the wall figure of merit:

`W_wall = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

The shell geometry and wall-shape factors cancel exactly, leaving

`W_wall = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).`

So the first explicit parent branch has a very strong simplification:

> on the matched thin-wall branch, `W_wall` is independent of the throat radius `a` and independent of the detailed wall-shape moment `I_f`.

It depends only on the wall-layer density, sound speed, amplitude, axial length, and shell width.

---

## 4. Identification with the Stage-058/059 fixed-point coupling

Stage 058 introduced the coupled support/source fixed-point strength `Xi`.

On the matched thin-wall branch, the exact Stage-064 gain is

`G_eq = g_phi^2 I_1 / K_X,`

with

`I_1 = int d^3y chi_phi(y)^2 / H(y).`

Therefore the fixed-point coupling is

`Xi = kappa G_eq = g_phi^2 I_1 L^2 / T_X.`

Using `g_phi = V0/ell`, `I_1 = N_(phi phi)/H_w`, and the explicit wall-shell `T_X` above gives

`Xi = 4 rho_w^2 V0^2 L^2 / (hbar^2 c_(s,w)^2 ell^2).`

Therefore

`Xi = W_wall`

exactly on the explicit matched thin-wall branch.

So the explicit branch solve has now tied together two quantities that had previously appeared in different parts of the program:

- the Stage-066 wall figure of merit,
- and the Stage-058/059 support/source fixed-point coupling.

They are the same object on this branch.

---

## 5. What Stage 070 changes

Before this stage, the explicit branch still looked as if it depended on a large set of independent support quantities.

After this stage, the first wall-shell branch fixes

- `T_X`,
- `K_X`,
- `kappa`,
- `J_1`,
- and `W_wall = Xi`

directly from the parent GNLS shell reduction.

That means the next phase is no longer to invent more branch coefficients.
It is to choose one concrete wall profile and one concrete mouth closure, and then evaluate the branch point in terms of a very small set of parent dimensionless ratios.
