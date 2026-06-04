# Moving-Throat PDE — Stage 057: Physical `(Pe, kappa, eta)` Placement Map for the Lowest Support Lane

## Purpose

Stage 056 converted the old abstract overlap-bias variable `alpha` into the physical transport Peclet number

`Pe = v_sigma L / D_sigma.`

The last remaining Stage-055 unknown is therefore the support-compliance variable previously written as `x`.

This stage converts the whole explicit non-twin family into a directly physical parameter map.

The main results are:

1. the old abstract support parameter is exactly

   `x = pi^2 / (kappa + pi^2/4),`

   where

   `kappa := K_X L^2 / T_X`

   is the baseline support stiffness ratio;
2. if the Robin mouth coefficient is written as

   `eta := hL = K_m L / T_X`,

   with `y(eta)` the unique root of `y tan y = eta`, then the support softening factor becomes

   `A_K(eta;kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2);`
3. the physical explicit lowest-lane ratio is therefore

   `zeta_0^(Pe+R)(Pe,eta;kappa)`
   `= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2);`
4. on the constructive branch `Pe >= 0`, this map is monotone

   - increasing in `Pe`,
   - decreasing in `eta`,
   - decreasing in `kappa`;
5. its exact closure ceiling is

   `zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa;`
6. therefore the Stage-052 demand is reachable on this first physical family **iff**

   `zeta_req <= zeta_max(kappa),`

   equivalently, when `zeta_req > pi^2/4`,

   `kappa <= pi^4 / [4(4 zeta_req - pi^2)].`

So the Stage-055 reachability problem is now fully rewritten in terms of three physical throat-operator ratios:

- `Pe` — axial support-source transport bias,
- `eta` — mouth Robin compliance,
- `kappa` — baseline support stiffness ratio.

---

## 1. Exact physical support parameters

Carry forward the finite-throat support operator with base stiffness `K_X`, axial tension `T_X`, and Robin mouth coefficient `h`.

Define the two physical dimensionless support ratios

`kappa := K_X L^2 / T_X,`

`eta := hL = K_m L / T_X`

(the last equality if the Robin boundary arises from a mouth spring `K_m`).

The mixed D/N lane has

`K_W^(eff) = K_X + pi^2 T_X/(4L^2)`
`          = (T_X/L^2) (kappa + pi^2/4),`

so the Stage-054 support parameter becomes

`x = pi^2 T_X / (L^2 K_W^(eff))`
`  = pi^2 / (kappa + pi^2/4).`

Thus the old abstract `x` is exactly the inverse stiffness ratio of the physical support operator.

---

## 2. Exact Robin softening in physical form

Let `y(eta)` be the unique root of

`y tan y = eta,`

with `0 < y < pi/2`.

Then the lowest Robin support lane has

`K_(phi,0)^(eff) = K_X + T_X y^2 / L^2`
`                = (T_X/L^2) (kappa + y^2).`

So the exact softening factor simplifies to

`A_K(eta;kappa)`
`= K_W^(eff)/K_(phi,0)^(eff)`
`= (kappa + pi^2/4)/(kappa + y(eta)^2).`

This is algebraically cleaner than the Stage-054 `x`-form and makes the physical support meaning explicit.

---

## 3. Exact physical explicit lowest-lane family

Combining Stage 056 with the Robin support branch gives the first fully physical explicit family:

`zeta_0^(Pe+R)(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

where

`Omega_Pe`
`= pi Pe (2 Pe exp(Pe) + pi)`
`  / [ (4 Pe^2 + pi^2)(exp(Pe)-1) ]`

and `y tan y = eta`.

This is the direct physical replacement for the Stage-055 abstract family

`zeta_0^(exp+R)(alpha,eta)`
`= Omega_exp(alpha)^2 / [1 - x/4 + x y(eta)^2/pi^2].`

The old variables are now identified as

`alpha = Pe,`

`x = pi^2/(kappa + pi^2/4).`

---

## 4. Exact monotone placement map

The physical family has a clean monotone structure.

### Increasing in transport bias `Pe`

From Stage 056,

`dOmega_Pe/dPe > 0`

on the constructive branch `Pe >= 0`, so

`partial_Pe zeta_0^(Pe+R) > 0.`

### Decreasing in compliance parameter `eta`

Because `y tan y = eta` has

`dy/deta = 1 / (tan y + y sec^2 y) > 0,`

one has

`partial_eta zeta_0^(Pe+R) < 0.`

So weaker mouth pinning (smaller `eta`) always helps the support lane.

### Decreasing in stiffness ratio `kappa`

Differentiating directly gives

`partial_kappa zeta_0^(Pe+R)`
`= Omega_Pe^2 ( y^2 - pi^2/4 ) / (kappa + y^2)^2 < 0,`

since `0 < y < pi/2`.

So a softer baseline support branch (smaller `kappa`) always helps.

This means the physical branch placement is completely ordered:

- more axial source drift helps,
- more mouth compliance helps,
- less baseline support stiffness helps.

---

## 5. Exact closure window in physical form

On the constructive branch `Pe >= 0`, the exact baseline point is

`Pe = 0,`

`eta = +infinity,`

so

`zeta_0^(Pe+R) = 1.`

The closure ceiling is reached at

`Pe -> +infinity,`

`eta -> 0^+,`

which gives

`zeta_max(kappa)`
`= (pi^2/4) (kappa + pi^2/4) / kappa.`

Therefore the exact constructive-branch window is

`1 <= zeta_0^(Pe+R)(Pe,eta;kappa) <= zeta_max(kappa)`

in closure.

This is the first explicit reachability ceiling written entirely in physical operator variables.

---

## 6. Exact reachability criterion and stiffness ceiling

The Stage-052 support demand is reachable on this physical family exactly when

`zeta_req <= zeta_max(kappa)`

that is,

`zeta_req <= (pi^2/4)(kappa + pi^2/4)/kappa.`

When `zeta_req > pi^2/4`, this is equivalent to the exact stiffness ceiling

`kappa <= kappa_max(zeta_req)`
`:= pi^4 / [4(4 zeta_req - pi^2)].`

So above the pure-overlap ceiling, the question is no longer “maybe enough compliance exists.”
It is the exact physical inequality above.

---

## 7. Exact parameter thresholds inside the physical family

Because the map is monotone in each physical parameter, the threshold surfaces are exact.

### Required overlap/transport threshold at fixed `(kappa,eta)`

The support demand requires

`Omega_Pe^2 >= Omega_req^2(kappa,eta;zeta_req)`

with

`Omega_req^2 := zeta_req (kappa + y(eta)^2)/(kappa + pi^2/4).`

Since the constructive `Omega_Pe` branch runs continuously from `1` to `pi/2`, a physical transport solution exists whenever

`1 <= Omega_req^2 <= pi^2/4`,

and the least constructive transport bias is the least nonnegative root of

`Omega_Pe^2 = Omega_req^2.`

### Required compliance threshold at fixed `(Pe,kappa)`

Equivalently,

`y(eta)^2 <= y_req^2(Pe,kappa;zeta_req)`

with

`y_req^2 := (Omega_Pe^2/zeta_req)(kappa + pi^2/4) - kappa.`

So if `0 <= y_req < pi/2`, the exact Robin threshold is

`eta <= eta_req := y_req tan(y_req).`

### Required stiffness threshold at fixed `(Pe,eta)`

If overlap alone is not already enough, i.e. if `Omega_Pe^2 < zeta_req`, then the exact stiffness ceiling is

`kappa <= kappa_req(Pe,eta;zeta_req)`

with

`kappa_req`
`= [ Omega_Pe^2 pi^2/4 - zeta_req y(eta)^2 ] / [ zeta_req - Omega_Pe^2 ],`

provided the numerator is nonnegative.

So all three physical threshold surfaces are now explicit.

---

## 8. Best current theorem statement after Stage 057

### What is exact now

- the old abstract overlap parameter is the physical transport Peclet number
  `Pe = v_sigma L / D_sigma`,
- the old abstract support parameter is
  `x = pi^2/(kappa + pi^2/4)`
  with `kappa = K_X L^2/T_X`,
- the Robin mouth parameter is
  `eta = hL = K_m L / T_X`,
- the physical explicit lowest-lane family is

  `zeta_0^(Pe+R)(Pe,eta;kappa)`
  `= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)`,

- this map is monotone increasing in `Pe` and monotone decreasing in both `eta` and `kappa`,
- and the exact constructive-branch ceiling is

  `zeta_max(kappa) = (pi^2/4)(kappa + pi^2/4)/kappa`.

### What this means physically

The lowest-lane support problem is no longer formulated in abstract deformation parameters at all.

It has collapsed to three concrete moving-throat operator ratios:

- axial source Peclet number `Pe`,
- mouth compliance `eta`,
- baseline support stiffness ratio `kappa`.

So the remaining gap is now as narrow as it can be without the full PDE:

> compute the physical branch point `(Pe, eta, kappa)` from the completed moving-throat operator and check whether it satisfies the exact threshold `zeta_req <= zeta_max(kappa)` and the corresponding monotone placement inequalities above.
