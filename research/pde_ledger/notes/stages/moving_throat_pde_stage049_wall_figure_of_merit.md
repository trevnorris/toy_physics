# Moving-Throat PDE — Stage 49: Dimensionless Wall Figure of Merit for the First Explicit Parent Branch

## Purpose

Stages 47–48 turned the parent support/source theorem gap into an explicit wall-amplitude test on the first thin-wall confinement branch.

The next useful simplification is to compress that explicit branch into the smallest possible control parameter.

This stage does that.

For the first explicit thin-wall wall family,

`V_conf(r;a) = V0 f((r-a)/ell),`

with the equilibrium-aligned source/support branch from Stage 47, define the dimensionless wall figure of merit

`W_wall := 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

Then the entire support/source success problem collapses to the exact comparison

`W_wall <= Pe_req / Delta_inf`  -> fail,

`W_wall >= Pe_req / Delta_0`    -> succeed,

and only the narrow intermediate band still needs the full fixed-point solve.

So the first explicit parent branch is no longer controlled by a diffuse set of parameters. It is controlled by one dimensionless number.

---

## 1. Definition of the wall figure of merit

Stage 48 showed that the thin-wall gain of the first explicit parent confinement family is

`G_eq^(tw) = 4 pi a^2 V0^2 J_1 / (K_X ell),`

with

- `a` the throat radius,
- `L` the throat length,
- `V0` the wall amplitude,
- `ell` the wall width,
- `J_1 = int dxi f'(xi)^2 / H(xi)` the weighted wall-profile moment,
- and `K_X` the axial support stiffness.

Using the Stage-44 geometry parameter

`kappa = K_X L^2 / T_X,`

a natural dimensionless wall control variable is

`W_wall := kappa G_eq^(tw)`
`        = 4 pi a^2 L^2 J_1 V0^2 / (T_X ell).`

This is the smallest parent quantity that still carries all of the support/source information relevant for the thin-wall matched branch.

---

## 2. Exact fail/succeed thresholds in wall form

The Stage-44 operator theorem gave

`G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ].`

Multiplying by `kappa` immediately yields the wall-threshold pair

`W_fail = Pe_req / Delta_inf(kappa,eta),`

`W_suff = Pe_req / Delta_0(kappa,eta).`

Therefore the first explicit parent branch satisfies the exact theorem:

- if `W_wall <= W_fail`, it fails,
- if `W_wall >= W_suff`, it succeeds,
- only if `W_fail < W_wall < W_suff` does one still need the full fixed-point solve.

So the wall-amplitude thresholds from Stage 48 are equivalent to one dimensionless wall figure-of-merit comparison.

---

## 3. Monotonicity of the wall control variable

The figure of merit is

`W_wall = 4 pi J_1 L^2 a^2 V0^2 / (T_X ell).`

So it is strictly monotone in the physically relevant directions:

- larger wall amplitude `V0` increases `W_wall`,
- larger throat radius `a` increases `W_wall`,
- larger throat length `L` increases `W_wall`,
- thinner wall width `ell` increases `W_wall`,
- larger support tension `T_X` suppresses `W_wall`,
- and larger weighted wall-profile moment `J_1` increases `W_wall`.

This makes the branch geometry immediately interpretable.

---

## 4. Constant-compressibility wall form

If the active wall layer is also nearly constant in compressional stiffness,

`H(xi) ~ H_w,`

then Stage 48 gave

`J_1 = I_f / H_w,`

with

`I_f := int dxi f'(xi)^2.`

The wall figure of merit becomes

`W_H = 4 pi a^2 L^2 I_f V0^2 / (H_w T_X ell).`

The same exact theorem holds:

`W_H <= Pe_req / Delta_inf`  -> fail,

`W_H >= Pe_req / Delta_0`    -> succeed.

So in the matched, constant-compressibility wall layer, the entire parent branch comparison is a direct contest between one support/loading number `W_H` and one demanded transport window set by the axial functions `Delta_inf`, `Delta_0`.

---

## 5. What Stage 49 changes

Before this stage, the explicit branch comparison still looked like a threshold on the wall amplitude `V0`.

After this stage, the first explicit parent support/source family has a one-number theorem gate:

`W_wall` (or `W_H` in the constant-compressibility limit).

That is the cleanest current status of the reduced theorem program.

The remaining work is no longer to invent more algebra. It is to compute the actual moving-throat branch values of

- `a`, `L`,
- the active wall width `ell`,
- the support tension scale `T_X`,
- the wall moment `J_1` (or `I_f/H_w`),
- and the axial support functions `Delta_0`, `Delta_inf`,

and then evaluate whether the real branch lands below, within, or above the exact wall-control window.
