# Moving-Throat PDE — Stage 061: Microscopic Gain Thresholds and the Exact Operator Phase Diagram

## Purpose

Stage 060 replaced the phenomenological support/source strength by the explicit microscopic gain

`Xi_micro = chi_sigma Lambda_phi^2 L^2 / T_X.`

The next honest step is to compare that gain directly with the exact Stage-059 thresholds

`Xi_fail = Pe_req / Delta_inf(kappa,eta),`

`Xi_suff = Pe_req / Delta_0(kappa,eta).`

This stage does that and pushes the result one step further by isolating a cleaner dimensionless microscopic control parameter.

The main results are:

1. using

   `kappa = K_X L^2 / T_X,`

   the microscopic coupling can be rewritten as

   `Xi_micro = kappa G_micro,`

   where the exact dimensionless support/source gain is

   `G_micro := chi_sigma Lambda_phi^2 / K_X;`
2. the exact no-go / sufficiency thresholds therefore become

   `G_fail = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

   `G_suff = Pe_req / [ kappa Delta_0(kappa,eta) ],`

   with the exact phase diagram

   - fail if `G_micro <= G_fail`,
   - succeed if `G_micro >= G_suff`,
   - only the interval `G_fail < G_micro < G_suff` requires the full root solve;
3. equivalently, the original microscopic parameters satisfy the exact threshold surfaces

   `chi_sigma <= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_inf ]`  -> fail,

   `chi_sigma >= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_0 ]`    -> succeed,

   and, at fixed `chi_sigma`,

   `Lambda_phi^2 <= T_X Pe_req / [ chi_sigma L^2 Delta_inf ]`  -> fail,

   `Lambda_phi^2 >= T_X Pe_req / [ chi_sigma L^2 Delta_0 ]`    -> succeed;
4. in the soft-support limit `kappa -> 0^+`, the exact endpoint values are

   `Delta_0 -> 1/2,`

   `Delta_inf -> 1,`

   so the microscopic gain thresholds diverge as

   `G_fail ~ Pe_req / kappa,`

   `G_suff ~ 2 Pe_req / kappa;`
5. in the highly compliant-mouth limit `eta -> +infty`, the endpoint data simplify exactly to

   `Delta_0^(inf) = (1-sech(sqrt(kappa)))/kappa,`

   `Delta_inf^(inf) = tanh(sqrt(kappa))/sqrt(kappa),`

   and therefore

   `G_fail^(inf) = Pe_req / [ sqrt(kappa) tanh(sqrt(kappa)) ],`

   `G_suff^(inf) = Pe_req / [ 1 - sech(sqrt(kappa)) ];`
6. in the combined stiff-support / compliant-mouth regime `kappa >> 1`, `eta >> 1`, these reduce to

   `G_fail ~ Pe_req / sqrt(kappa),`

   `G_suff ~ Pe_req.`

So the explicit microscopic closure now has a very sharp interpretation:

- too-soft support (`kappa << 1`) requires unrealistically large gain,
- a sufficiently compliant mouth makes the sufficiency threshold saturate at an `O(Pe_req)` microscopic gain,
- and only the bounded intermediate band between `G_fail` and `G_suff` still requires the full fixed-point solve.

---

## 1. Exact microscopic gain `G_micro`

Stage 060 gave

`Xi_micro = chi_sigma Lambda_phi^2 L^2 / T_X.`

But Stage 057 already uses the support stiffness ratio

`kappa = K_X L^2 / T_X.`

Therefore

`Xi_micro = kappa ( chi_sigma Lambda_phi^2 / K_X ).`

Define the exact dimensionless support/source gain

`G_micro := chi_sigma Lambda_phi^2 / K_X.`

Then

`Xi_micro = kappa G_micro.`

This is the cleanest microscopic control parameter so far because it removes the explicit length and tension scales from the branch-strength comparison.

---

## 2. Exact operator phase diagram in microscopic variables

Stage 059 already proved that the operator-selected branch succeeds or fails according to the comparison of `Xi` with the exact thresholds

`Xi_fail = Pe_req / Delta_inf(kappa,eta),`

`Xi_suff = Pe_req / Delta_0(kappa,eta).`

Substitute `Xi_micro = kappa G_micro`. Then the exact microscopic thresholds are

`G_fail(kappa,eta) = Pe_req / [ kappa Delta_inf(kappa,eta) ],`

`G_suff(kappa,eta) = Pe_req / [ kappa Delta_0(kappa,eta) ].`

So the exact phase diagram is now

- `G_micro <= G_fail`  -> impossible inside this operator family,
- `G_micro >= G_suff`  -> guaranteed success,
- `G_fail < G_micro < G_suff` -> only then is the full root solve needed.

This is the first operator phase diagram written directly in microscopic support/source variables.

---

## 3. Threshold surfaces for `chi_sigma` and `Lambda_phi`

Undoing the definition of `G_micro` gives exact threshold surfaces in the original microscopic variables.

At fixed `Lambda_phi`, success and failure are controlled by

`chi_sigma <= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_inf(kappa,eta) ]`  -> fail,

`chi_sigma >= T_X Pe_req / [ Lambda_phi^2 L^2 Delta_0(kappa,eta) ]`    -> succeed.

At fixed `chi_sigma`, they are controlled by

`Lambda_phi^2 <= T_X Pe_req / [ chi_sigma L^2 Delta_inf(kappa,eta) ]`  -> fail,

`Lambda_phi^2 >= T_X Pe_req / [ chi_sigma L^2 Delta_0(kappa,eta) ]`    -> succeed.

So the microscopic theorem gap is no longer “somehow get a large enough `Xi`.” It is now a concrete competition between source susceptibility and support/source coupling on one side, and the exact support geometry functions `Delta_0`, `Delta_inf` on the other.

---

## 4. Exact soft-support limit `kappa -> 0^+`

Using the exact endpoint formulas from Stage 058,

`Delta_0(kappa,eta)`
`= eta (cosh(sqrt(kappa))-1)`
`  / [ kappa ( sqrt(kappa) sinh(sqrt(kappa)) + eta cosh(sqrt(kappa)) ) ],`

`Delta_inf(kappa,eta)`
`= [ cosh(sqrt(kappa)) + (eta/sqrt(kappa)) sinh(sqrt(kappa)) - 1 ]`
`  / [ sqrt(kappa) sinh(sqrt(kappa)) + eta cosh(sqrt(kappa)) ],`

the exact soft-support limits are

`lim_(kappa->0+) Delta_0 = 1/2,`

`lim_(kappa->0+) Delta_inf = 1.`

Therefore

`G_fail ~ Pe_req / kappa,`

`G_suff ~ 2 Pe_req / kappa`

as `kappa -> 0^+`.

So a very soft baseline support channel is strongly disfavored: the microscopic gain required for success diverges like `1/kappa`.

---

## 5. Exact highly compliant-mouth limit `eta -> +infty`

At fixed `kappa`, the mouth-compliant limit is exact.

From the same endpoint formulas,

`lim_(eta->+infty) Delta_0(kappa,eta) = (1-sech(sqrt(kappa)))/kappa,`

`lim_(eta->+infty) Delta_inf(kappa,eta) = tanh(sqrt(kappa))/sqrt(kappa).`

So the exact microscopic gain thresholds become

`G_fail^(inf)(kappa) = Pe_req / [ sqrt(kappa) tanh(sqrt(kappa)) ],`

`G_suff^(inf)(kappa) = Pe_req / [ 1 - sech(sqrt(kappa)) ].`

This is useful because it isolates the best-case mouth-compliance regime of the same operator family.

Two consequences are immediate.

### Small `kappa` inside the compliant-mouth branch

If `kappa << 1`, then

`tanh(sqrt(kappa)) ~ sqrt(kappa),`

`1 - sech(sqrt(kappa)) ~ kappa/2,`

so the exact soft-support divergence is recovered:

`G_fail^(inf) ~ Pe_req / kappa,`

`G_suff^(inf) ~ 2 Pe_req / kappa.`

### Stiff-support inside the compliant-mouth branch

If `kappa >> 1`, then

`tanh(sqrt(kappa)) -> 1,`

`sech(sqrt(kappa)) -> 0,`

so

`G_fail^(inf) ~ Pe_req / sqrt(kappa),`

`G_suff^(inf) ~ Pe_req.`

That means sufficiently strong baseline support stiffness is not itself the main obstacle once the mouth is compliant enough; the sufficiency threshold saturates at an `O(Pe_req)` microscopic gain.

---

## 6. What Stage 061 changes

The theorem gap has narrowed again.

Before Stage 061, the support question was:

> what is the physical `Xi`, and is it above or below the exact interval `[Xi_fail, Xi_suff]`?

After Stage 061, it is:

> what is the physical dimensionless gain `G_micro = chi_sigma Lambda_phi^2/K_X`, and does it lie above or below the exact geometry-controlled thresholds `G_fail(kappa,eta)` and `G_suff(kappa,eta)`?

That is a stronger result because it separates the operator problem into two transparent pieces:

- the support geometry sector `(kappa,eta)`, which sets the exact threshold functions,
- and the microscopic source/support gain `G_micro`, which the completed moving-throat operator must supply.

The sharpest practical lesson is now clear:

- very soft support (`kappa << 1`) is strongly unfavorable,
- a highly compliant mouth (`eta >> 1`) helps substantially,
- and on the best-case compliant branch, success is controlled by whether `G_micro` can reach an order-`Pe_req` value.

So the next honest step is to compute `G_micro = chi_sigma Lambda_phi^2/K_X` from a more explicit moving-throat branch and compare it directly to the exact threshold surfaces derived here.
