# Moving-Throat PDE — Stage 051: Exact Lowest-Twin Sufficiency Criterion on the Physical Tracking Branch

## Purpose

Stage 050 reduced the explicit D/N support comparison to one sharp question:

> does the physical moving-throat branch lie in the regime `zeta_req <= 1`, so that the symmetric lowest twin support lane is already sufficient?

That statement is still useful, but it is written in terms of the reduced support ratio `zeta_req`.
The next honest step is to eliminate `zeta_req` and translate the test into the continuum variables already frozen in Stages 047–050.

This stage does that exactly.

The main result is that the lowest symmetric twin lane is sufficient **iff**
the actual selected-branch point satisfies

`Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1-eps) / pi^2,`

where

`Pi_tr := F_tr G_tr`

is the exact tracking-branch product function. Equivalently, the full twin-sufficiency question collapses to a single radiative-demand inequality.

So the support problem is now no longer “do we have enough support in some vague sense?”
It is a one-line theorem test on the physical branch product.

---

## 1. Exact tracking-branch product

From the earlier stages,

`R_target = F_tr(xi_req,delta;R_tr),`

`M_req    = G_tr(xi_req,delta;R_tr),`

with

`G_tr(xi,delta;R)`
`= 9 xi (xi + delta) / [ 9 delta + (9 + 2 R^2) xi ],`

`F_tr(xi,delta;R)`
`= [ 9 delta + (9 + 2 R^2) xi ]^2 [ 9 delta + (9 + 2 R) xi ]^2`
`  / [ 81 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

Multiplying them gives the exact branch product

`Pi_tr(xi,delta;R)`
`:= F_tr(xi,delta;R) G_tr(xi,delta;R)`

`= xi (xi + delta) [ 9 delta + (9 + 2 R) xi ]^2 [ 9 delta + (9 + 2 R^2) xi ]`
`  / [ 9 (1 - xi) ( 9 delta^2 + 18 delta xi + (9 + 2 R^2) xi^2 )^2 ].`

This is the unique product that matters for the lowest-twin question.

Two endpoint facts are exact:

`Pi_tr(0,delta;R) = 0,`

`lim_(xi -> 1^-) Pi_tr(xi,delta;R) = +infinity.`

So for every finite microscopic support scale there is at least one stable-side depth where the lowest-twin criterion is saturated.

---

## 2. Elimination of `zeta_req`

Stage 048 gave the exact support-demand map

`S_req := M_req / M_mix,`

`zeta_req = (S_req - 1) / [ 1 + eps (S_req - 2) ],`

with the mixed-only product law from Stage 047

`R_target M_mix = 8 Lambda (1 - eps) / pi^2.`

Since

`R_target = F_tr(xi_req,delta;R_tr),`

`M_req    = G_tr(xi_req,delta;R_tr),`

it follows that

`S_req = M_req / M_mix`
`     = [ R_target M_req ] / [ R_target M_mix ]`
`     = Pi_tr(xi_req,delta;R_tr) / [ 8 Lambda (1 - eps) / pi^2 ].`

Now impose the exact Stage-050 twin threshold

`zeta_req <= 1 <=> S_req <= 2.`

Then the lowest symmetric twin lane is sufficient **iff**

`Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1 - eps) / pi^2.`

That is the sharpest version yet of the support theorem.

---

## 3. Exact threshold scales

The same inequality can be written as exact thresholds for several microscopic variables.

### Radiative-demand threshold

Define the exact lowest-twin required radiative scale

`Lambda_twin,req`
`:= (pi^2 / [16 (1 - eps)]) Pi_tr(xi_req,delta;R_tr).`

Then

`Lambda >= Lambda_twin,req`
`<=>`
`zeta_req <= 1.`

### Mixed-baseline threshold

Since `S_req <= 2` is equivalent to `M_req <= 2 M_mix`, the exact mixed-only baseline threshold is simply

`M_mix >= M_mix^(twin,req)`
`:= (1/2) G_tr(xi_req,delta;R_tr).`

So the lowest twin succeeds exactly when the mixed branch already supplies at least half of the total load demanded by the physical selected point.

### Wall-to-mixed overlap threshold

Using the Stage-047 coherent map

`M_mix = 8 Z_W (1 + chi_0)^2 / [ pi^2 (1 - eps_eta) (1 - eps) ],`

the exact wall/mixed overlap threshold is

`Z_W >= Z_W^(twin,req)`

with

`Z_W^(twin,req)`
`:= pi^2 (1 - eps_eta) (1 - eps) G_tr(xi_req,delta;R_tr)`
`   / [ 16 (1 + chi_0)^2 ].`

So even before the full PDE is solved, the lowest-twin support test is already equivalent to an explicit lower bound on the wall-to-mixed overlap strength.

---

## 4. Exact twin-saturation depth at fixed mixed baseline

The twin-sufficiency boundary corresponds to

`M_req = 2 M_mix,`

i.e.

`G_tr(xi,delta;R_tr) = 2 M_mix.`

Because `G_tr` is strictly increasing on the stable branch, this equation has a unique stable-side root. Solving the resulting quadratic gives

`xi_(2x)(M_mix,delta;R)`
`= [ 2 M_mix (9 + 2 R^2) - 9 delta`
`    + sqrt( [2 M_mix (9 + 2 R^2) - 9 delta]^2 + 648 M_mix delta ) ] / 18.`

This is the exact depth at which the lowest twin lane is just saturated by the mixed baseline.

It is useful because it shows the twin threshold is not a loose scaling estimate.
It lives on an exact algebraic branch of the same selected-mode family.

---

## 5. Best current theorem statement after Stage 051

The lowest-twin question is no longer phrased in terms of an abstract support ratio.

### What is exact now

- the physical selected-branch product is

  `Pi_tr = F_tr G_tr`,
- the mixed-only product scale is

  `8 Lambda (1 - eps) / pi^2`,
- the lowest symmetric twin lane succeeds **iff**

  `Pi_tr(xi_req,delta;R_tr) <= 16 Lambda (1 - eps) / pi^2`,
- the exact required radiative scale is

  `Lambda_twin,req = pi^2 Pi_tr / [16(1-eps)]`,
- the exact required mixed baseline is

  `M_mix^(twin,req) = G_tr/2`,
- and the exact twin-saturation depth at fixed mixed baseline is the closed quadratic root `xi_(2x)` above.

### What this means physically

The support problem has collapsed again.

The first explicit coherent twin lane is sufficient or insufficient according to a **single branch product inequality**.
So the next honest question is no longer “how do we parameterize support?”
It is:

> **does the completed moving-throat operator produce enough mixed baseline / radiative scale to satisfy the exact twin-sufficiency product test, or must the physical lowest support lane become non-twin before the normalization target can be met?**
