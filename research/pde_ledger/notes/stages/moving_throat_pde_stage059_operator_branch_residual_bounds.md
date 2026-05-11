# Moving-Throat PDE — Stage 059: Exact Residual Bounds on the Operator-Selected Branch

## Purpose

Stage 41 turned the transport bias `Pe` into the root of a coupled support/source fixed-point equation,

`Pe = Xi Delta(Pe;kappa,eta),`

and trapped every constructive branch point inside the exact interval

`Xi Delta_0(kappa,eta) <= Pe_* <= Xi Delta_inf(kappa,eta).`

That means the Stage-40 physical placement map can now be evaluated on a **true operator-selected branch** rather than on an arbitrary external `Pe`.

This stage uses that bracket to derive exact success/no-go tests for the normalization residual without solving the full root.

The main results are:

1. the operator-selected physical lowest-lane ratio is

   `zeta_phys(Xi,eta;kappa)`
   `= zeta_0^(Pe_*+R)`
   `= Omega_(Pe_*)^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

   with `y tan y = eta`;
2. because `zeta_0^(Pe+R)` is monotone increasing in `Pe`, the exact Stage-41 branch bracket immediately gives

   `zeta_-(Xi,eta;kappa) <= zeta_phys(Xi,eta;kappa) <= zeta_+(Xi,eta;kappa),`

   where

   `zeta_- = Omega_(Xi Delta_0)^2 (kappa + pi^2/4)/(kappa + y^2),`

   `zeta_+ = Omega_(Xi Delta_inf)^2 (kappa + pi^2/4)/(kappa + y^2);`
3. therefore the exact residual

   `R_phys := zeta_req - zeta_phys`

   obeys

   `R_- <= R_phys <= R_+,`

   with

   `R_- = zeta_req - zeta_+,`

   `R_+ = zeta_req - zeta_-;`
4. this yields two exact theorem gates **before** solving the full branch root:

   - guaranteed success if `zeta_- >= zeta_req`,
   - guaranteed failure inside this operator family if `zeta_+ < zeta_req`;
5. inside the reachable Stage-40 window, define the unique constructive branch point `Pe_req` by

   `Omega_(Pe_req)^2 = zeta_req (kappa + y^2)/(kappa + pi^2/4)`

   (`Pe_req` exists exactly when the Stage-40 ceiling holds);
6. then the exact operator-coupling thresholds are

   `Xi_fail = Pe_req / Delta_inf(kappa,eta),`

   `Xi_suff = Pe_req / Delta_0(kappa,eta),`

   and they satisfy

   `Xi_fail <= Xi_suff`;
7. so the coupled operator now has a precise three-zone structure:

   - `Xi <= Xi_fail` : impossible,
   - `Xi >= Xi_suff` : guaranteed,
   - `Xi_fail < Xi < Xi_suff` : full root solve needed, but only inside this narrow interval.

This is the first exact residual-bracketing theorem on the operator-selected branch.

---

## 1. Operator-selected physical lowest-lane ratio

Carry forward the Stage-40 physical family

`zeta_0^(Pe+R)(Pe,eta;kappa)`
`= Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2),`

with `y tan y = eta` and the constructive-branch monotonicity

`dOmega_Pe/dPe > 0.`

Now evaluate it on the Stage-41 branch point `Pe_*`, defined implicitly by

`Pe_* = Xi Delta(Pe_*;kappa,eta).`

The operator-selected physical support ratio is therefore

`zeta_phys(Xi,eta;kappa)`
`= Omega_(Pe_*)^2 (kappa + pi^2/4)/(kappa + y(eta)^2).`

This is the first branch-selected version of the Stage-40 placement map.

---

## 2. Exact branch brackets for the physical ratio

Because Stage 41 proved

`Xi Delta_0 <= Pe_* <= Xi Delta_inf`

and Stage 40 already proved that `zeta_0^(Pe+R)` is strictly increasing in `Pe`, one gets the exact operator-side support bracket

`zeta_-(Xi,eta;kappa) <= zeta_phys(Xi,eta;kappa) <= zeta_+(Xi,eta;kappa),`

where

`zeta_-(Xi,eta;kappa)`
`= Omega_(Xi Delta_0)^2 (kappa + pi^2/4)/(kappa + y^2),`

`zeta_+(Xi,eta;kappa)`
`= Omega_(Xi Delta_inf)^2 (kappa + pi^2/4)/(kappa + y^2).`

So the full root solve is not needed to get a rigorous placement interval for the physical branch.

---

## 3. Exact residual bracket

Define the physical support residual

`R_phys(Xi,eta,kappa;zeta_req)`
`:= zeta_req - zeta_phys(Xi,eta;kappa).`

Then the branch bracket gives

`R_- <= R_phys <= R_+,`

with

`R_-(Xi,eta,kappa;zeta_req)`
`:= zeta_req - zeta_+(Xi,eta;kappa),`

`R_+(Xi,eta,kappa;zeta_req)`
`:= zeta_req - zeta_-(Xi,eta;kappa).`

This is the exact operator-selected residual envelope.

It has two immediate consequences.

### Guaranteed success

If

`zeta_-(Xi,eta;kappa) >= zeta_req,`

then the lower branch bracket already clears the target, so

`R_phys <= 0`

for every constructive branch root.

### Guaranteed failure in this operator family

If

`zeta_+(Xi,eta;kappa) < zeta_req,`

then even the upper branch bracket cannot reach the target, so

`R_phys > 0`

for every constructive branch root in this coupled operator family.

So the full fixed-point solve is only needed in the intermediate window.

---

## 4. Exact coupling thresholds `Xi_fail` and `Xi_suff`

Inside the Stage-40 reachability window, define `Pe_req` as the unique constructive solution of

`Omega_(Pe_req)^2 = zeta_req (kappa + y^2)/(kappa + pi^2/4).`

Because `Omega_Pe` is strictly increasing from `1` to `pi/2`, this `Pe_req` exists iff the Stage-40 ceiling is satisfied.

Then the exact coupling thresholds are

`Xi_fail = Pe_req / Delta_inf(kappa,eta),`

`Xi_suff = Pe_req / Delta_0(kappa,eta).`

Their meaning is exact.

### No-go threshold

If `Xi <= Xi_fail`, then

`Xi Delta_inf <= Pe_req,`

hence

`zeta_+(Xi,eta;kappa) <= zeta_req`,

so the coupled operator family cannot reach the target.

### Guaranteed-success threshold

If `Xi >= Xi_suff`, then

`Xi Delta_0 >= Pe_req,`

hence

`zeta_-(Xi,eta;kappa) >= zeta_req`,

so every constructive branch root reaches the target.

Since `Delta_inf >= Delta_0 > 0`, one has

`Xi_fail <= Xi_suff.`

So the physical theorem gap is now reduced to a bounded coupling window, not a wide-open parameter hunt.

---

## 5. Weak-coupling branch law for the residual

From Stage 39,

`Omega_Pe = 1 + ((4-pi)/(2pi)) Pe + O(Pe^2),`

so

`Omega_Pe^2 = 1 + ((4-pi)/pi) Pe + O(Pe^2).`

Combining this with the Stage-41 weak-coupling branch law

`Pe_* = Xi Delta_0 + O(Xi^2)`

gives

`zeta_phys(Xi,eta;kappa)`
`= A_K(eta;kappa)`
`  [ 1 + ((4-pi)/pi) Xi Delta_0(kappa,eta) + O(Xi^2) ],`

where

`A_K(eta;kappa) = (kappa + pi^2/4)/(kappa + y(eta)^2).`

So the operator-selected physical branch departs from the symmetric support point with a completely explicit linear slope in the coupling strength `Xi`.

---

## 6. Best current theorem statement after Stage 42

### What is exact now

- the physical support branch is no longer just the Stage-40 family; it is that family evaluated on a genuine coupled-operator branch point `Pe_*`;
- every constructive branch root obeys the exact interval

  `Xi Delta_0 <= Pe_* <= Xi Delta_inf`;
- the physical support ratio and the exact residual therefore obey rigorous lower and upper bounds;
- and the full support question has collapsed to one narrow coupling window between the exact thresholds

  `Xi_fail = Pe_req / Delta_inf,`

  `Xi_suff = Pe_req / Delta_0.`

### What this means physically

The remaining theorem gap is no longer “derive everything about the lowest support lane.”
It is now:

> compute the real moving-throat coupling strength `Xi`, compare it to the exact interval `[Xi_fail, Xi_suff]`, and then solve the fixed-point root only if the real branch lands inside that already narrow intermediate window.
