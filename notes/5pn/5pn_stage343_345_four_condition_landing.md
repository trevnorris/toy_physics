# 5PN continuation — Stages 343–345: landing on the four reduced finish-line conditions

This session takes the reduced finish line from Stages 340–342 and asks the next obvious question:

> can the current reduced moving-throat hierarchy actually land on the four surviving conditions simultaneously, or is there still a hidden algebraic incompatibility?

The answer is better than expected.

Within the present reduced hierarchy there is **no algebraic contradiction** among the four surviving conditions. The whole problem splits cleanly into:

1. a three-equation **orbit-lock landing surface** for
   - `delta ln R_tr = 0`,
   - `delta ln R_target = 0`,
   - `delta ln epsilon_eta = 0`,
2. and a one-equation **outgoing landing surface** for
   - `N_Q = 1`, equivalently `chi_Q = 1` on the natural source-map branch.

So the remaining obstruction is now purely PDE-side branch realization, not another reduced-sector inconsistency.

---

## Stage 343 — exact orbit-lock landing surface

On the coherent branch,

a) the direct tracking observable is

a) `R_tr = (1 + chi_0/(1+delta_U)) / (1 + chi_0)`,

b) the target observable is

b) `R_target = Lambda (1-epsilon_eta)(1-epsilon)^2 / [ Z_W (1+chi_0)^2 ]`,

c) and the third direct observable is simply

c) `epsilon_eta`.

The exact first-order landing conditions are therefore

1. `delta ln R_tr = 0`,
2. `delta ln R_target = 0`,
3. `delta ln epsilon_eta = 0`.

The script shows these are exactly equivalent to the solved co-drift system

`dln delta_U = - (1 + delta_U)/(1 + chi_0) * dln chi_0`,

`dln Lambda = dln Z_W + 2 chi_0 dln chi_0/(1 + chi_0) + 2 epsilon dln epsilon/(1 - epsilon)`,

`dln epsilon_eta = 0`.

So the first three finish-line conditions live on one explicit orbit-lock surface in the coherent branch variables.

A useful structural fact is that this whole orbit packet is support-blind at the reduced level:

`d_zeta ln R_tr = d_zeta ln R_target = d_zeta ln epsilon_eta = 0`.

So support enhancement still matters for the isotropic normalization baseline, but it does **not** move the branch on or off the weak-axisymmetric orbit-lock surface.

---

## Stage 344 — exact outgoing landing surface

On the natural source-map branch,

`N_Q = 1 / chi_Q`.

So the fourth finish-line condition is exactly

`N_Q = 1  <=>  chi_Q = 1`.

That already shows one exact outgoing landing surface: the canonical outgoing branch with `chi_Q = 1`.

The session also rewrote the first-order Family-1 outgoing defect in the form

`Delta_Q = - [sigma_* / (1 - sigma_*)] Xi_slip deltaPi_tan`,

so that

`N_Q - 1 = 1/(1 + Delta_Q) - 1`.

Therefore, on the exact lower parent compensation family, the first-order outgoing defect vanishes automatically whenever

`Xi_slip = 0`.

Then

`Delta_Q = 0`,

`N_Q - 1 = 0`.

So there are two exact reduced ways to land the fourth condition:

1. the canonical outgoing branch `chi_Q = 1`,
2. or, at first order, the lower parent compensation family with `Xi_slip = 0`.

---

## Stage 345 — the combined four-condition landing theorem

Putting the two pieces together gives the current cleanest reduced theorem:

### Orbit-lock surface

`dln delta_U = - (1 + delta_U)/(1 + chi_0) * dln chi_0`,

`dln Lambda = dln Z_W + 2 chi_0 dln chi_0/(1 + chi_0) + 2 epsilon dln epsilon/(1 - epsilon)`,

`dln epsilon_eta = 0`.

### Outgoing surface

Either

`chi_Q = 1`,

or, on the first-order lower compensation family,

`Xi_slip = 0`.

### Combined reduced finish line

If the branch lies on the orbit-lock surface and also lies on either outgoing surface above, then all four reduced finish-line conditions vanish simultaneously:

`delta ln R_tr = 0`,

`delta ln R_target = 0`,

`delta ln epsilon_eta = 0`,

`N_Q - 1 = 0`.

So the four-condition finish line is now shown to be **algebraically reachable** inside the current reduced hierarchy.

---

## Best current verdict

This session removes the last serious worry that there might still be some hidden reduced-sector contradiction among the four finish-line conditions.

There is not.

What remains open is exactly what the compact moving-throat master has been saying:

- whether the completed PDE actually realizes the orbit-lock surface,
- and whether it actually realizes the canonical / lower-compensation outgoing surface.

So the open problem is now purely one of **branch realization**, not another reduced-algebra obstruction.
