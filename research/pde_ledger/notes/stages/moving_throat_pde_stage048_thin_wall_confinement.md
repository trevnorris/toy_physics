# Moving-Throat PDE — Stage 48: Explicit Thin-Wall Confinement Branch and Parent Thresholds for the Wall Amplitude

## Purpose

Stage 47 showed that the parent equilibrium branch already aligns the source and support channels and turns the exact gain into

`G_eq = g_phi^2 I_1 / K_X,`

with

`I_1 = int d^3y chi_phi(y)^2 / H(y).`

The next honest step is to evaluate that gain on the first explicit parent wall family instead of leaving `g_phi` and `I_1` abstract.

This stage does that for the natural moving-wall confinement form

`V_conf(r;a) = V0 f((r-a)/ell),`

where:

- `a` is the throat radius,
- `ell` is the wall thickness / active-layer width,
- `V0` is the wall amplitude,
- and `f` is a fixed dimensionless wall-shape profile.

The main results are:

1. the support loading amplitude is exactly

   `g_phi = V0 / ell;`
2. the exact shell integral entering the equilibrium gain is

   `I_1 = 4 pi ell [ a^2 J_1 + 2 a ell J_2 + ell^2 J_3 ],`

   where

   `J_n := int dxi xi^n f'(xi)^2 / H(xi);`
3. for a centered symmetric wall layer, `J_2 = 0`, so

   `I_1 = 4 pi ell [ a^2 J_1 + ell^2 J_3 ];`
4. the exact equilibrium gain becomes

   `G_eq = 4 pi V0^2 [ a^2 J_1 / ell + 2 a J_2 + ell J_3 ] / K_X;`
5. in the thin-wall limit `ell << a`, the leading gain is

   `G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell);`
6. comparing this with the Stage-44 fail/succeed surfaces gives the explicit wall-amplitude thresholds

   `V0_fail^2 = K_X ell G_fail / (4 pi a^2 J_1),`

   `V0_suff^2 = K_X ell G_suff / (4 pi a^2 J_1);`
7. after inserting

   `G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

   `G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

   `kappa = K_X L^2 / T_X,`

   the explicit `K_X` dependence cancels:

   `V0_fail^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_inf ],`

   `V0_suff^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_0 ];`
8. if the active layer also has nearly constant compressibility,

   `H(xi) ~ H_w,`

   `J_1 = I_f / H_w,`

   with

   `I_f := int dxi f'(xi)^2,`

   so the thresholds reduce to

   `V0_fail^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_inf ],`

   `V0_suff^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_0 ].`

So the Stage-46 parent-overlap theorem has now been converted into an explicit wall-amplitude test on the first concrete parent confinement branch.

---

## 1. The explicit wall family

Take the parent confinement near the throat wall to be

`V_conf(r;a) = V0 f((r-a)/ell),`

with dimensionless wall coordinate

`xi := (r-a)/ell.`

A support displacement `a -> a + phi(s)` gives

`delta V_conf = - (partial_a V_conf) phi(s)`
`            = + (V0/ell) f'(xi) phi(s).`

Comparing with the Stage-45/47 loading form

`delta V_conf = - g_phi chi_phi(y) phi(s),`

one identifies the branch data

`g_phi = V0/ell,`

`chi_phi(y) = f'(xi)`

up to the overall sign convention of `chi_phi`.

So the wall amplitude and wall thickness directly determine the parent loading strength.

---

## 2. Exact shell integral for the equilibrium gain

Stage 47 showed that the exact equilibrium gain is

`G_eq = g_phi^2 I_1 / K_X,`

with

`I_1 = int d^3y chi_phi(y)^2 / H(y).`

For a spherical shell-like wall in the three transverse coordinates of the parent reduction, the volume element is

`d^3y = 4 pi r^2 dr = 4 pi ell (a + ell xi)^2 dxi.`

Using `chi_phi = f'(xi)` gives

`I_1 = 4 pi ell int dxi (a + ell xi)^2 f'(xi)^2 / H(xi)`
`    = 4 pi ell [ a^2 J_1 + 2 a ell J_2 + ell^2 J_3 ],`

with the wall-profile moments

`J_1 := int dxi f'(xi)^2 / H(xi),`

`J_2 := int dxi xi f'(xi)^2 / H(xi),`

`J_3 := int dxi xi^2 f'(xi)^2 / H(xi).`

If the active layer is centered symmetrically around the nominal wall location, the odd moment vanishes,

`J_2 = 0,`

and then

`I_1 = 4 pi ell [ a^2 J_1 + ell^2 J_3 ].`

---

## 3. Exact and thin-wall gains

Insert `g_phi = V0/ell` into the exact equilibrium gain:

`G_eq = g_phi^2 I_1 / K_X`
`     = 4 pi V0^2 [ a^2 J_1 / ell + 2 a J_2 + ell J_3 ] / K_X.`

For the centered branch `J_2 = 0`,

`G_eq = 4 pi V0^2 [ a^2 J_1 / ell + ell J_3 ] / K_X.`

So in the thin-wall regime `ell << a`, the leading term is

`G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell).`

This is the first explicit parent scaling law for the support/source gain. It says the branch is helped by:

- larger wall amplitude `V0`,
- larger wall area `~ a^2`,
- thinner active wall width `ell`,
- larger weighted wall-profile moment `J_1`,
- and smaller baseline support stiffness `K_X`.

---

## 4. Exact wall-amplitude fail/succeed thresholds

The Stage-44 phase diagram is

`G_eq <= G_fail`  -> fail,

`G_eq >= G_suff`  -> succeed.

Using the thin-wall leading gain gives the direct wall-amplitude thresholds

`V0_fail^2 = K_X ell G_fail / (4 pi a^2 J_1),`

`V0_suff^2 = K_X ell G_suff / (4 pi a^2 J_1).`

So the first explicit parent wall family no longer speaks in terms of abstract gain. It speaks directly in terms of the physical wall amplitude `V0`.

---

## 5. Cancellation of `K_X` after inserting the operator geometry law

Now insert the Stage-44 operator formulas

`G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

with

`kappa = K_X L^2 / T_X.`

Then the `K_X` in the prefactor cancels exactly, leaving

`V0_fail^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_inf(kappa,eta) ],`

`V0_suff^2 = T_X ell Pe_req / [ 4 pi a^2 L^2 J_1 Delta_0(kappa,eta) ].`

This is structurally important. Once the support geometry functions are inserted, the thin-wall wall-amplitude thresholds are controlled by:

- the support tension scale `T_X`,
- the wall width `ell`,
- the branch geometry `a` and `L`,
- the demanded transport Peclet `Pe_req`,
- and the wall-profile overlap moment `J_1`.

The baseline support stiffness no longer appears explicitly in the prefactor.

---

## 6. Constant-compressibility wall layer

If the active wall layer is also nearly constant in compressional stiffness,

`H(xi) ~ H_w,`

then

`J_1 = I_f / H_w,`

with

`I_f := int dxi f'(xi)^2.`

So the wall-amplitude thresholds become

`V0_fail^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_inf ],`

`V0_suff^2 = H_w T_X ell Pe_req / [ 4 pi a^2 L^2 I_f Delta_0 ].`

Using the Stage-45 identity

`H_w = h'(rho_w) = m c_(s,w)^2 / rho_w,`

this may also be written as

`V0_fail^2 = m c_(s,w)^2 T_X ell Pe_req / [ 4 pi a^2 rho_w L^2 I_f Delta_inf ],`

`V0_suff^2 = m c_(s,w)^2 T_X ell Pe_req / [ 4 pi a^2 rho_w L^2 I_f Delta_0 ].`

So the explicit wall family now ties the theorem gap to a very concrete question:

> is the actual parent wall amplitude `V0` on the moving-throat branch large enough, for the actual wall width `ell` and active-layer compressibility, to clear these thresholds?

---

## 7. What Stage 48 changes

Before this stage, the theorem gap was still phrased in terms of a parent loading amplitude `g_phi` and an overlap/coherence structure.

After this stage, the first explicit parent branch has collapsed that gap to a physical wall-amplitude test.

The remaining branch data are now:

- wall amplitude `V0`,
- wall width `ell`,
- throat radius `a`,
- throat length `L`,
- wall-profile moments such as `I_f` or `J_1`,
- and the axial support geometry functions `Delta_0(kappa,eta)`, `Delta_inf(kappa,eta)`.

That is a much narrower place to be.

The next honest step is therefore no longer another abstract overlap theorem. It is to evaluate the actual moving-throat branch values of

- `V0`,
- `ell`,
- `a`, `L`,
- and the axial support functions,

and test them against the explicit wall-amplitude fail/succeed surfaces derived here.
