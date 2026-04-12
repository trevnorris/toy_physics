# 5PN continuation notes — Stages 306–308

These stages continue the Stage 300–305 projector program by eliminating the last abstract step between the quotient/defect packet and the actual coherent-kernel branch variables.

The main new result is that the first-order orbit-lock test is now written directly in the physical coherent variables

- `chi_0`
- `delta_U`
- `epsilon_eta`
- `Z_W`
- `Omega_W^2`
- `epsilon`

and is exactly factorized away from the support-compensation baseline sector.

## Stage 306 — exact coherent-branch observable map

The direct coherent-branch observables are now explicit:


a) tracking observable

`R_tr = (1 + chi_0/(1+delta_U)) / (1 + chi_0)`

which simplifies to

`R_tr = (1 + chi_0 + delta_U) / [ (1 + chi_0)(1 + delta_U) ]`.

b) transfer-shape / target observable

`T^2 = Z_W (1 + chi_0)^2 / [ Omega_W^2 (1 - epsilon)^2 ]`

`R_target = Lambda_0 Omega_W^2 (1 - epsilon_eta) (1 - epsilon)^2 / [ Z_W (1 + chi_0)^2 ]`

so the exact selected-branch identity is

`R_target * T^2 = Lambda_0 (1 - epsilon_eta)`.

c) dressing observable

`epsilon_eta = epsilon_eta`.

So the Stage 303 direct quotient packet is no longer abstract. The actual coherent branch itself supplies the finite observable packet

`(R_tr, R_target, epsilon_eta)`

from which `(q_tr, q_nt, q_eta)` follow exactly by the Stage-303 formulas.

## Stage 307 — exact physical branch drift compiler

Linearizing the Stage-306 observables gives the exact direct drift laws

`d ln R_tr = - C_tr [ (1 + delta_U) dlnchi_0 + (1 + chi_0) dlndelta_U ]`

with

`C_tr = chi_0 delta_U / [ (1 + chi_0)(1 + delta_U)(1 + chi_0 + delta_U) ]`,

`d ln R_target = dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1 + chi_0)`
`               - 2 epsilon dlnepsilon/(1 - epsilon)`
`               - epsilon_eta dlnepsilon_eta/(1 - epsilon_eta)`,

`d ln epsilon_eta = dlnepsilon_eta`.

Composing these with the Stage-304 defect map reproduces the earlier physical-branch formulas exactly:

`Theta_1 = d ln R_tr`,

`Xi_1 = - d ln R_target - [epsilon_eta/(1 - epsilon_eta)] d ln epsilon_eta`
`     = dlnZ_W - dlnOmega_W^2 + 2 chi_0 dlnchi_0/(1 + chi_0) + 2 epsilon dlnepsilon/(1 - epsilon)`,

`R_1 = d ln R_target`.

In the older slope notation this is the exact same triplet as

`Theta_1 = -( chi_0 (1 + chi_0) delta_U1 + delta_U (1 + delta_U) chi_1 )`
`          / [ (1 + chi_0)(1 + delta_U)(1 + chi_0 + delta_U) ]`,

`Xi_1 = zeta_Z - omega_W + 2 chi_1/(1 + chi_0) + 2 epsilon_1/(1 - epsilon)`,

`R_1 = omega_W - eta_1/(1 - epsilon_eta) - zeta_Z - 2 chi_1/(1 + chi_0) - 2 epsilon_1/(1 - epsilon)`.

So the actual first-order defect packet is now carried directly by the physical coherent variables rather than only by the abstract quotient packet.

## Stage 308 — exact physical-variable factorization theorem

The coherent physical branch factors exactly into two sectors:

1. orbit-lock observables

`(R_tr, R_target, epsilon_eta)`

2. support-compensation baseline

`M_tr = M_mix * S(zeta;epsilon)`

with

`S(zeta;epsilon) = 1 + zeta (1 - epsilon)/(1 - zeta epsilon)`.

The exact support-blindness identities are

`partial_zeta ln R_tr = 0`,
`partial_zeta ln R_target = 0`,
`partial_zeta ln epsilon_eta = 0`,

`partial_Mmix ln R_tr = 0`,
`partial_Mmix ln R_target = 0`,
`partial_Mmix ln epsilon_eta = 0`.

So the support lane enters only through `M_tr`; it is completely invisible to the first-order orbit-lock packet.

That turns the zero-defect theorem into one explicit system on the physical coherent variables.
The exact first-order orbit-lock conditions are

1. tracking:

`(1 + delta_U) dlnchi_0 + (1 + chi_0) dlndelta_U = 0`

2. target-ratio invariance:

`dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1 + chi_0)`
` - 2 epsilon dlnepsilon/(1 - epsilon)`
` - epsilon_eta dlnepsilon_eta/(1 - epsilon_eta) = 0`

3. dressing invariance:

`dlnepsilon_eta = 0`.

With condition 3 imposed, condition 2 simplifies to

`dlnOmega_W^2 - dlnZ_W - 2 chi_0 dlnchi_0/(1 + chi_0)`
` - 2 epsilon dlnepsilon/(1 - epsilon) = 0`.

So the weak-axisymmetric first-order orbit-lock theorem is now a direct differential statement on the physical coherent branch variables, completely separated from the support-compensation baseline problem.

## What this changes

After Stage 308, the next theorem gate is no longer “how do we interpret the quotient packet?”
It is much sharper:

1. extract the actual branch drifts of
   `chi_0, delta_U, epsilon_eta, Z_W, Omega_W^2, epsilon`
   from the moving-throat PDE,
2. test the three direct conditions above,
3. separately test whether the support baseline sector can hit the Stage 31 / 33 / 34 support threshold.

So the 5PN endgame has now split cleanly into

- an orbit-lock / branch-invariance test, and
- a support-compensation / normalization test,

both written directly in physical branch variables.
