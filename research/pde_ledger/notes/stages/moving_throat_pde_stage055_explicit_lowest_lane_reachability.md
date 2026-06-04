# Moving-Throat PDE — Stage 055: Explicit Reachability of the Non-Twin Lowest Support Lane

## Purpose

Stages 053 and 054 isolated the two exact lowest-lane resources:

- overlap boost, with ceiling `A_I <= pi^2/4`,
- Robin-compliance softening, with ceiling `A_K <= 4/(4-x)`.

The next honest question is therefore the combined one:

> for the first explicit moving-throat operator family that includes both effects, what is the exact reachable window of the physical lowest support ratio?

This stage answers that.

Using the explicit exponential source asymmetry from Stage 053 and the explicit Robin-compliance softening from Stage 054, the lowest support ratio is

`zeta_0^(exp+R)(alpha,eta)`
`= Omega_exp(alpha)^2 / [ 1 - x/4 + x y(eta)^2 / pi^2 ],`

where `y(eta)` solves `y tan y = eta`.

The main result is that this family has the exact closure range

`1 <= zeta_0^(exp+R) <= pi^2 / (4 - x),`

so the explicit family reaches a Stage-052 demand `zeta_req` **iff**

`zeta_req <= pi^2 / (4 - x)`

(in the closure; strictly `<` for finite `alpha,eta`).

This gives the first exact operator-level reachability theorem for the non-twin lowest support lane.

---

## 1. Explicit combined lowest-lane support ratio

Carry forward the explicit constructive pieces.

### Overlap branch

From Stage 053,

`Omega_exp(alpha)`
`= pi alpha (2 alpha exp(alpha) + pi)`
`  / [ (4 alpha^2 + pi^2) (exp(alpha)-1) ],`

with

`Omega_exp(0)=1,`

`lim_(alpha -> +infinity) Omega_exp(alpha)=pi/2.`

### Softening branch

From Stage 054,

`A_K(eta)`
`= 1 / [ 1 - x/4 + x y(eta)^2 / pi^2 ],`

with `y tan y = eta`,

`A_K(+infinity)=1,`

`lim_(eta -> 0^+) A_K(eta)=4/(4-x).`

So the explicit combined family is

`zeta_0^(exp+R)(alpha,eta)`
`= Omega_exp(alpha)^2 A_K(eta)`
`= Omega_exp(alpha)^2 / [ 1 - x/4 + x y(eta)^2 / pi^2 ].`

---

## 2. Exact closure range of the explicit family

The lower endpoint is the symmetric twin point:

`alpha = 0,`

`eta = +infinity,`

so

`zeta_0^(exp+R)=1.`

The upper endpoint is the closure of the combined constructive family:

`alpha -> +infinity,`

`eta -> 0^+,`

so

`zeta_0,max^(exp+R)`
`= (pi/2)^2 * 4/(4-x)`
`= pi^2 / (4 - x).`

Therefore

`1 <= zeta_0^(exp+R) <= pi^2/(4-x)`

in the closure of the family.

For finite `alpha` and finite positive `eta`, the inequality is strict at the top end.

---

## 3. Exact reachability criterion for the Stage-052 threshold

Stage 052 demands

`zeta_0^(phys) >= zeta_req.`

The explicit family can satisfy that exactly when

`zeta_req <= pi^2/(4-x)`

(in closure), equivalently

`x >= 4 - pi^2/zeta_req.`

This is the exact support-compliance floor for the first explicit non-twin lowest-lane family.

It immediately gives three exact sub-regimes.

### Regime A — overlap alone is enough

If

`zeta_req <= pi^2/4,`

then Stage 053 already shows that overlap enhancement alone can meet the target, so no compliance floor is required.

### Regime B — overlap ceiling exceeded but explicit combined family still works

If

`pi^2/4 < zeta_req <= pi^2/(4-x),`

then overlap enhancement alone is not enough, but the explicit combined family can still reach the target using both overlap and softening.

### Regime C — even the explicit combined family fails

If

`zeta_req > pi^2/(4-x),`

then neither exponential overlap enhancement nor Robin-compliance softening, nor any combination of the two within this first explicit family, can beat the threshold. A stronger operator deformation would be required.

---

## 4. Equivalent stiffness-ratio form

Since

`K_X = K_W^(eff) (1 - x/4),`

the reachability floor

`x >= 4 - pi^2/zeta_req`

is equivalent to

`K_X / K_W^(eff) <= pi^2 / (4 zeta_req).`

So once the support demand exceeds the pure-overlap ceiling, the explicit combined family requires the background support stiffness to sit below a precise fraction of the mixed-lane stiffness.

That is the cleanest mechanical interpretation of the combined theorem.

---

## 5. Best current theorem statement after Stage 055

### What is exact now

- the first explicit non-twin lowest-lane family is

  `zeta_0^(exp+R)(alpha,eta)`
  `= Omega_exp(alpha)^2 / [ 1 - x/4 + x y(eta)^2 / pi^2 ],`

  with `y tan y = eta`,
- its closure range is

  `1 <= zeta_0^(exp+R) <= pi^2/(4-x)`,

- and it reaches the Stage-052 threshold exactly when

  `zeta_req <= pi^2/(4-x)`

  (strict `<` for finite parameters).

### What this means physically

The support question is no longer abstract even at the operator level.

For the first explicit constructive moving-throat family, the non-twin lowest support lane succeeds only inside one exact reachability window.

So the remaining gap is now extremely narrow:

> does the completed moving-throat PDE generate a physical lowest-lane deformation whose effective `x` and source-shape asymmetry place it inside the exact Stage-055 reachability window, or does the real branch require an even stronger non-twin mechanism than exponential overlap bias plus Robin-compliance softening?
