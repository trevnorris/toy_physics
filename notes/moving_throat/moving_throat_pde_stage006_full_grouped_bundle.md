# Moving-Throat PDE — Stage 6: Full Grouped Real `P2` Bundle, Exact Projectors, and the Isotropic Normalization Test

## Purpose

Stage 5 isolated the exact algebraic object that still has to hit the universal quadrupole target:

`mhat_0^2 P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

But Stage 5 still spoke mostly in one-lane language. It showed how a single grouped lane talks to the outgoing `l=2` port, and it translated the result into the grouped-`P2` normalization stack. What it did **not** yet do was write the full coupled `20/21/22` bundle in one place and show exactly how the microscopic wall/BdG/Maxwell/mixed data feed

- the conservative grouped coefficients `D_{A0}, D_{A2}, D_{A4}`,
- the outgoing-transfer coefficients `N_{A0}, N_{A2}, N_{A4}`,
- the grouped isotropy defects,
- and the final isotropic normalization product.

That is the whole job of this stage.

The main new outputs are:

1. an exact weighted-projector calculus for the grouped real `P2` bundle,
2. explicit full-bundle formulas for the coupled low-frequency coefficients,
3. a clean separation of the microscopic sources of anisotropy,
4. the exact isotropic-branch normalization test in terms of the full coupled bundle,
5. and first-order anisotropy transport laws showing how conservative or transfer anisotropy leaks into the grouped response and outgoing prefactor.

So this stage is not yet the full nonlinear PDE theorem. It is the first exact bookkeeping layer in which the **entire** reduced `20/21/22` bundle is carried at once.

---

## 1. Exact weighted projector calculus for the grouped real `P2` bundle

The grouped real `P2` bundle is not an ordinary Euclidean triple. The natural grouped metric is

`Ggrp = diag(1,2,2)`.

This is the same weighting that already appeared in the grouped trace/anomaly formulas:

`xbar = (x_20 + 2 x_21 + 2 x_22)/5`,

`a_x = (2 x_20 - x_21 - x_22)/10`,

`b_x = (x_21 - x_22)/2`.

A particularly useful `Ggrp`-orthogonal basis is

`e_bar = (1,1,1)`,

`e_a   = (4,-1,-1)`,

`e_b   = (0,1,-1)`.

These vectors satisfy

`e_bar^T Ggrp e_a = 0`,

`e_bar^T Ggrp e_b = 0`,

`e_a^T   Ggrp e_b = 0`,

with norms

`e_bar^T Ggrp e_bar = 5`,

`e_a^T   Ggrp e_a   = 20`,

`e_b^T   Ggrp e_b   = 4`.

So for any grouped coefficient vector

`x = (x_20, x_21, x_22)^T`,

we have the exact decomposition

`x = xbar e_bar + a_x e_a + b_x e_b`.

Equivalently, the three exact projectors are

`P_bar = (1/5) [[1,2,2],[1,2,2],[1,2,2]]`,

`P_a   = (1/20) [[16,-8,-8],[-4,2,2],[-4,2,2]]`,

`P_b   = (1/4) [[0,0,0],[0,2,-2],[0,-2,2]]`,

so that

`P_bar + P_a + P_b = I3`,

`P_i P_j = delta_ij P_i`.

This is the cleanest exact way to organize grouped-lane isotropy.

- `P_bar` extracts the isotropic grouped trace,
- `P_a` and `P_b` extract the two independent anisotropy defects.

So at any stage of the coupled derivation, the grouped-lane isotropy condition is simply

`P_a x = P_b x = 0`.

---

## 2. Full coupled lane-by-lane reduced bundle

At linear order on an isotropic reference throat, the grouped real `20/21/22` channels do not mix directly with one another. Each lane is coupled internally, but the three lanes remain separate copies of the same reduced mechanism unless anisotropy is introduced.

So for each grouped lane

`A in {20,21,22}`

we keep:

- a wall/worldtube amplitude `q_A`,
- stable BdG support modes `X_(A,alpha)` with frequencies `varpi_(A,alpha)`,
- one or more localized brane-like gauge coordinates `U_(A,r)` with frequencies `Omega_(U,A,r)`,
- one or more mixed `A_w/F_{mu w}/J^w` coordinates `W_(A,r)` with frequencies `Omega_(W,A,r)`,
- internal gauge-sector mixing `R_(A,r)`.

The reduced quadratic Lagrangian in lane `A` is therefore

`L_A`
`= (1/2) M_A qdot_A^2 - (1/2) K_A q_A^2`
`  + sum_alpha [ (1/2) Xdot_(A,alpha)^2 - (1/2) varpi_(A,alpha)^2 X_(A,alpha)^2 + c_(A,alpha) q_A X_(A,alpha) ]`
`  + sum_r [ (1/2) Udot_(A,r)^2 - (1/2) Omega_(U,A,r)^2 U_(A,r)^2`
`          + (1/2) Wdot_(A,r)^2 - (1/2) Omega_(W,A,r)^2 W_(A,r)^2`
`          + R_(A,r) U_(A,r) W_(A,r)`
`          + g_(U,A,r) q_A U_(A,r) + g_(W,A,r) q_A W_(A,r) ]`.

This is the full reduced grouped-lane bundle we need for the next theorem gate.

It combines, in one place,

- the geometry wall backbone,
- the BdG support self-energy,
- the conservative localized-Maxwell/mixed self-energy,
- and the outgoing transfer channel.

---

## 3. Exact full-bundle low-frequency coefficients

### 3.1 BdG support contributions

For each lane `A`, define the exact BdG moments

`B_(A,0) = sum_alpha c_(A,alpha)^2 / varpi_(A,alpha)^2`,

`B_(A,2) = sum_alpha c_(A,alpha)^2 / varpi_(A,alpha)^4`,

`B_(A,4) = sum_alpha c_(A,alpha)^2 / varpi_(A,alpha)^6`.

These are just the Stage-3 conservative support moments collected lane by lane.

### 3.2 Conservative Maxwell/mixed contributions

For each port-active internal pair `(U_(A,r),W_(A,r))`, define

`Delta_(A,r) = Omega_(U,A,r)^2 Omega_(W,A,r)^2 - R_(A,r)^2`,

`S_(A,r) = Omega_(U,A,r)^2 + Omega_(W,A,r)^2`,

`Q_(A,r) = g_(U,A,r)^2 Omega_(W,A,r)^2 + 2 g_(U,A,r) g_(W,A,r) R_(A,r) + g_(W,A,r)^2 Omega_(U,A,r)^2`,

`G_(A,r) = g_(U,A,r)^2 + g_(W,A,r)^2`.

Then the conservative low-frequency coefficients contributed by that port are

`Z_(A,0)^(r) = Q_(A,r) / Delta_(A,r)`,

`Z_(A,2)^(r) = [ Q_(A,r) S_(A,r) - G_(A,r) Delta_(A,r) ] / Delta_(A,r)^2`,

`Z_(A,4)^(r) = [ Q_(A,r)(S_(A,r)^2 - Delta_(A,r)) - S_(A,r) G_(A,r) Delta_(A,r) ] / Delta_(A,r)^3`.

Summing over all internal Maxwell/mixed branches gives

`Z_(A,n) = sum_r Z_(A,n)^(r)`

for `n=0,2,4`.

### 3.3 Outgoing-transfer coefficients

For each port-active internal pair, the exact outgoing transfer factor begins with

`N_(A,0)^(r) = [ Omega_(U,A,r)^2 g_(W,A,r) + R_(A,r) g_(U,A,r) ]^2 / Delta_(A,r)^2`.

At the next two even orders,

`N_(A,2)^(r)`
`= 2 P_(A,r) [ P_(A,r) S_(A,r) - Delta_(A,r) g_(W,A,r) ] / Delta_(A,r)^3`,

`N_(A,4)^(r)`
`= [ Delta_(A,r)^2 g_(W,A,r)^2 - 2 Delta_(A,r) P_(A,r)^2 - 4 Delta_(A,r) P_(A,r) S_(A,r) g_(W,A,r) + 3 P_(A,r)^2 S_(A,r)^2 ] / Delta_(A,r)^4`,

where

`P_(A,r) = Omega_(U,A,r)^2 g_(W,A,r) + R_(A,r) g_(U,A,r)`.

Summing over outgoing ports gives

`N_(A,n) = sum_r N_(A,n)^(r)`

for `n=0,2,4`.

### 3.4 Total conservative wall operator coefficients

Putting everything together, the full coupled conservative grouped-lane operator is

`D_A^(cons)(omega) = D_(A,0) + D_(A,2) omega^2 + D_(A,4) omega^4 + O(omega^6)`

with

`D_(A,0) = K_A - B_(A,0) - Z_(A,0)`,

`D_(A,2) = -[ M_A + B_(A,2) + Z_(A,2) ]`,

`D_(A,4) = -[ B_(A,4) + Z_(A,4) ]`.

These are the exact full-bundle low-frequency coefficients the previous stage was asking for.

Nothing symbolic is left hidden here:

- `B_(A,n)` is the BdG support part,
- `Z_(A,n)` is the conservative localized-Maxwell/mixed part,
- `N_(A,n)` is the outgoing-transfer part.

So each grouped lane now has fully explicit microscopic reduced coefficients.

---

## 4. Exact grouped trace/anomaly bookkeeping for the microscopic coefficients

For any lane family `X_(A,n)` — for example `K_A`, `M_A`, `B_(A,0)`, `Z_(A,2)`, or `N_(A,0)` — define

`Xbar_n = (X_(20,n) + 2 X_(21,n) + 2 X_(22,n))/5`,

`a_(X,n) = (2 X_(20,n) - X_(21,n) - X_(22,n))/10`,

`b_(X,n) = (X_(21,n) - X_(22,n))/2`.

Then the full coupled coefficients decompose exactly as

`Dbar_0 = Kbar - Bbar_0 - Zbar_0`,

`a_(D,0) = a_K - a_(B,0) - a_(Z,0)`,

`b_(D,0) = b_K - b_(B,0) - b_(Z,0)`,

`Dbar_2 = -[ Mbar + Bbar_2 + Zbar_2 ]`,

`a_(D,2) = -[ a_M + a_(B,2) + a_(Z,2) ]`,

`b_(D,2) = -[ b_M + b_(B,2) + b_(Z,2) ]`,

`Dbar_4 = -[ Bbar_4 + Zbar_4 ]`,

`a_(D,4) = -[ a_(B,4) + a_(Z,4) ]`,

`b_(D,4) = -[ b_(B,4) + b_(Z,4) ]`.

Similarly, for the outgoing-transfer bundle,

`Nbar_n = (N_(20,n) + 2 N_(21,n) + 2 N_(22,n))/5`,

`a_(N,n) = (2 N_(20,n) - N_(21,n) - N_(22,n))/10`,

`b_(N,n) = (N_(21,n) - N_(22,n))/2`.

This is the first exact microscopic answer to the question:

> where can grouped-lane anisotropy enter?

Answer: only through anisotropy in the wall baseline (`K_A,M_A`), the BdG support data (`B_(A,n)`), or the Maxwell/mixed conservative and outgoing-transfer bundles (`Z_(A,n), N_(A,n)`).

That is the complete reduced anisotropy ledger.

---

## 5. Isotropic branch and the exact full-bundle normalization test

The natural isotropic grouped-lane branch is the one on which all grouped anisotropy defects vanish:

`a_(D,0)=b_(D,0)=a_(D,2)=b_(D,2)=a_(D,4)=b_(D,4)=0`,

`a_(N,0)=b_(N,0)=a_(N,2)=b_(N,2)=a_(N,4)=b_(N,4)=0`.

Then the three grouped lanes collapse to common coefficients

`D_(20,n)=D_(21,n)=D_(22,n)=D_n`,

`N_(20,n)=N_(21,n)=N_(22,n)=N_n`.

On that branch, the exact Stage-5 formulas become

`u_2 = - D_2 / D_0`,

`u_4 = (D_2^2 - D_0 D_4) / D_0^2`,

`P_0 = N_0 / D_0`,

`P_2 = (D_0 N_2 - 2 D_2 N_0) / D_0^2`,

`P_4 = [ D_0^2 N_4 - 2 D_0 ( D_2 N_2 + D_4 N_0 ) + 3 D_2^2 N_0 ] / D_0^3`.

The universal normalization condition is therefore

`mhat_0^2 P_0 = 54 G c_s^5 / (5 a^5 c^5)`.

Using the explicit full-bundle coefficient definitions above, this becomes

`mhat_0^2 * N_0 / ( K - B_0 - Z_0 ) = 54 G c_s^5 / (5 a^5 c^5)`.

This is the exact isotropic normalization test for the full coupled reduced bundle.

So the remaining theorem gap is now completely sharp in reduced microscopic language:

> does the true moving-throat branch produce coupled bundle data `(K, M, B_n, Z_n, N_n)` whose isotropic static ratio `N_0/(K-B_0-Z_0)` lands on the universal target?

That is the whole question.

---

## 6. Constant-prefactor branch conditions

The 2.5PN normalization package repeatedly singles out the especially simple branch on which the outgoing prefactor is constant at the orders we retain.

In the present full-bundle notation, that means

`P_2 = 0`,

`P_4 = 0`.

The exact algebraic conditions are therefore

`N_2 = 2 D_2 N_0 / D_0`,

`N_4 = [ 2 D_0 ( D_2 N_2 + D_4 N_0 ) - 3 D_2^2 N_0 ] / D_0^2`.

So the constant-prefactor outgoing branch does **not** require the higher coefficients to vanish separately. It requires the outgoing-transfer moments `N_2,N_4` to be correlated with the conservative bundle moments `D_0,D_2,D_4` in exactly the way above.

That is another exact reduced theorem target for the completed PDE.

---

## 7. First-order anisotropy transport laws

One of the most useful things we can do at this stage is linearize around the isotropic branch.

Write

`D_(A,n) = D_n + delta D_(A,n)`,

`N_(A,n) = N_n + delta N_(A,n)`

with `delta D` and `delta N` small grouped-lane anisotropy corrections.

Then the first two objects that matter most are:

### 7.1 Response-moment anisotropy

Since

`u_2^(A) = - D_(A,2) / D_(A,0)`,

the first-order variation is

`delta u_2^(A) = -[ delta D_(A,2) + u_2 delta D_(A,0) ] / D_0`.

So the grouped anisotropy defects of the conservative response obey

`a_(u,2) = -[ a_(D,2) + u_2 a_(D,0) ] / D_0`,

`b_(u,2) = -[ b_(D,2) + u_2 b_(D,0) ] / D_0`.

This shows that even if the `omega^2` coefficient itself were isotropic, anisotropy in the static coefficient `D_0` would still leak into the normalized grouped response.

### 7.2 Outgoing-prefactor anisotropy

Since

`P_0^(A) = N_(A,0) / D_(A,0)`,

the first-order variation is

`delta P_0^(A) = [ delta N_(A,0) - P_0 delta D_(A,0) ] / D_0`.

Therefore

`a_(P,0) = [ a_(N,0) - P_0 a_(D,0) ] / D_0`,

`b_(P,0) = [ b_(N,0) - P_0 b_(D,0) ] / D_0`.

This is the most useful anisotropy diagnostic of the stage.

It says the outgoing-prefactor isotropy can fail in only two ways:

1. direct anisotropy in the outgoing transfer bundle `N_(A,0)`,
2. conservative anisotropy in `D_(A,0)` that is reweighted by the isotropic prefactor `P_0`.

So when the next PDE step is computed, we will know immediately whether any anisotropy is entering through the support side, the Maxwell/mixed transfer side, or both.

---

## 8. Expert reading of the normalization formula

The exact isotropic normalization product

`P_0 = N_0 / D_0`

has a simple physical reading once the branch is stable (`D_0 > 0`).

- Increasing the static outgoing-transfer weight `N_0` raises `P_0` directly.
- Increasing the conservative support renormalizations `B_0` or `Z_0` lowers the denominator `D_0 = K-B_0-Z_0` and therefore also raises `P_0`.

So if the natural branch undershoots the universal target, there are only two broad ways to move toward it:

1. increase the outgoing transfer strength,
2. or soften the static wall operator through conservative support dressing while keeping the branch stable.

This does **not** solve the theorem, but it tells us how the microscopic bundle parameters have to move if the completed PDE is to land on the right normalization.

---

## 9. What this stage achieves physically

This stage closes the first full-bundle bookkeeping gap.

### 9.1 It gives the exact coupled grouped-lane coefficients

The bundle coefficients `D_(A,0), D_(A,2), D_(A,4)` and `N_(A,0), N_(A,2), N_(A,4)` are now written explicitly in terms of

- wall data,
- BdG support data,
- conservative localized-Maxwell/mixed data,
- and outgoing-transfer data.

So the next theorem gate is not missing algebra anymore.

### 9.2 It localizes the sources of anisotropy exactly

The projector calculus makes it impossible to hide where grouped-lane anisotropy comes from.
Each microscopic sector carries its own exact `(bar,a,b)` defects, and the full bundle inherits them linearly.

### 9.3 It turns the final normalization problem into one ratio

On the natural isotropic branch, everything collapses to the exact ratio

`P_0 = N_0 / D_0`.

That is the reduced microscopic version of the remaining 2.5PN/4PN normalization bottleneck.

---

## 10. Immediate next step

The right next derivation is now even narrower than before.

Do **not** reopen the whole theory at once.
Instead:

1. compute the actual overlap integrals that determine the full-bundle data `B_(A,n), Z_(A,n), N_(A,n)` on the true moving-throat branch,
2. project them with `(P_bar, P_a, P_b)`,
3. check whether the grouped anisotropy defects vanish on the natural branch,
4. and then test the single ratio

`mhat_0^2 N_0 / (K-B_0-Z_0)`

against the universal target.

That is now the smallest honest next theorem gate.
