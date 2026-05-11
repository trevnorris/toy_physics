# Moving-Throat PDE — Stage 054: Robin-Compliance Softening of the Lowest Support Lane

## Purpose

Stage 36 put an exact ceiling on the overlap resource:

`A_I <= pi^2/4.`

So the next honest step is to make the second Stage-35 resource explicit:

> can the lowest support lane become softer than the mixed D/N lane by a concrete finite-throat mechanism, and by how much?

This stage answers that with the simplest explicit support-lane deformation:

- keep the same bulk operator and the same Neumann bottom,
- but replace the Dirichlet mouth by a finite Robin compliance.

The main results are:

1. the lowest support wavenumber is no longer fixed at `pi/(2L)` but is the unique root `y/L` of

   `y tan y = eta,`

   with `eta := h L` the dimensionless mouth-compliance parameter;
2. the support-lane stiffness becomes

   `K_(phi,0)^(eff) = K_X + T_X y^2 / L^2,`

   which is strictly smaller than the mixed-lane D/N value whenever `0 < eta < +infinity`;
3. the exact softening factor is

   `A_K(eta) := K_W^(eff) / K_(phi,0)^(eff)`
   `= 1 / [ 1 - x/4 + x y^2 / pi^2 ],`

   with `x := pi^2 T_X / (L^2 K_W^(eff))` and `0 < x < 4`;
4. it ranges from

   `A_K = 1` at the symmetric D/N point to the exact maximum

   `A_K,max = 4/(4-x)`

   in the soft-mouth limit `eta -> 0^+`;
5. therefore **pure support softening alone** can rescue the Stage-35 threshold only if

   `zeta_req <= 4/(4-x)`.

So this stage makes the support-softening side of the non-twin budget just as explicit as the overlap side.

---

## 1. Explicit compliant lowest support lane

Take the same finite interval `s in [0,L]` and the same bulk support operator, but impose

- Robin condition at the mouth `s=0`,
- Neumann condition at the bottom `s=L`.

Write the axial support mode as `psi(s)` satisfying

`psi'' + k^2 psi = 0,`

`psi'(0) = h psi(0),`

`psi'(L) = 0,`

with `h > 0` the mouth-compliance coefficient.

Solving the boundary-value problem gives the exact characteristic equation

`k tan(kL) = h.`

Defining the dimensionless variables

`y := kL,`

`eta := hL,`

the lowest support branch is the unique root

`y tan y = eta,`

with

`0 < y < pi/2`.

The symmetric D/N limit is recovered as `eta -> +infinity`, where `y -> pi/2`.
The fully soft-mouth limit is `eta -> 0^+`, where `y -> 0`.

---

## 2. Exact lowest-lane support stiffness and softening factor

Let the undeformed same-operator family have base stiffness `K_X` and axial tension `T_X`.

Then the mixed D/N lane has the frozen stiffness

`K_W^(eff) = K_X + pi^2 T_X / (4L^2),`

while the Robin-deformed lowest support lane has

`K_(phi,0)^(eff) = K_X + T_X y^2 / L^2.`

So the exact support-softening factor is

`A_K(eta) := K_W^(eff) / K_(phi,0)^(eff).`

It is convenient to introduce the dimensionless stiffness/tension ratio

`x := pi^2 T_X / (L^2 K_W^(eff)).`

Since `K_X > 0`, one has

`K_X = K_W^(eff) (1 - x/4),`

hence

`0 < x < 4.`

Substituting gives the exact reduced form

`A_K(eta)`
`= 1 / [ 1 - x/4 + x y(eta)^2 / pi^2 ].`

---

## 3. Exact softening window

The endpoint values are immediate.

### Symmetric D/N point

At `eta -> +infinity`, `y -> pi/2`, so

`A_K = 1.`

### Soft-mouth limit

At `eta -> 0^+`, `y -> 0`, so

`A_K,max = 1 / (1 - x/4) = 4 / (4 - x).`

Thus the exact softening window is

`1 <= A_K < 4/(4-x)`

for finite `eta > 0`, with the upper endpoint reached only in the closure `eta -> 0^+`.

Because the map `y -> eta = y tan y` is strictly increasing on `(0,pi/2)`, and `A_K` is strictly decreasing in `y`, the softening factor is strictly decreasing in `eta`.
So weaker mouth pinning means stronger support softening.

---

## 4. Exact pure-softening rescue criterion

Stage 35 showed that, at equal overlap,

`A_K >= zeta_req`

is the exact rescue condition.

Stage 37 now bounds the left-hand side by

`A_K <= 4/(4-x)`

(with closure equality as `eta -> 0^+`).

Therefore pure support softening alone is possible only if

`zeta_req <= 4/(4-x).`

Equivalently, for `zeta_req > 1`, the support-tension parameter must satisfy

`x >= 4 - 4/zeta_req.`

So once the demanded support ratio exceeds the symmetric twin value, the first explicit softening branch supplies an exact compliance floor.

---

## 5. Exact Robin threshold at fixed `zeta_req`

From

`A_K = 1 / [ 1 - x/4 + x y^2 / pi^2 ] >= zeta_req,`

one obtains the exact eigenvalue bound

`y^2 <= y_req^2`

with

`y_req^2 := (pi^2/x) ( 1/zeta_req - 1 + x/4 ).`

Therefore the exact Robin-compliance threshold is

`eta <= eta_req := y_req tan(y_req)`

whenever the right-hand side is nonnegative.

If

`1/zeta_req - 1 + x/4 < 0,`

then `y_req^2 < 0` and pure softening is impossible on this branch, regardless of the Robin parameter.

So the support-softening question is no longer qualitative either. It is an exact root-placement problem on the finite throat.

---

## 6. Best current theorem statement after Stage 37

### What is exact now

- the Robin-deformed lowest support branch is determined by

  `y tan y = eta`,

- the exact support-softening factor is

  `A_K = 1 / [ 1 - x/4 + x y^2 / pi^2 ]`,

- it obeys the exact window

  `1 <= A_K <= 4/(4-x)`

  (closure at `eta -> 0^+`),
- and pure softening rescue is possible only if

  `zeta_req <= 4/(4-x)`.

### What this means physically

The second Stage-35 resource is now explicit.

The lowest support lane can indeed become softer than the mixed D/N lane, but only within a finite operator-controlled window set by `x`.

So the remaining question is already narrower:

> if the physical branch needs more than the Stage-36 overlap ceiling `pi^2/4`, does the moving-throat operator also supply enough Robin-type support softening to bridge the remaining gap?
