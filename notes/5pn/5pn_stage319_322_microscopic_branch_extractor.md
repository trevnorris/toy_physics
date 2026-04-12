# 5PN continuation notes — Stages 319–322

These stages take the Stage-315–318 reduced continuum orbit tester and complete the missing bridge back to the actual microscopic coherent-kernel variables.

The net effect is that the weak-axisymmetric orbit-lock problem is no longer formulated in an abstract five-drift packet. It is now directly readable from the microscopic grouped drifts of

- `lambda_W`
- `c_etaU`
- `gamma`
- `K_U`
- `K_eta^(eff)`
- `K_W^(eff)`
- `mu_W`
- `T_U`

while the explicit support variables `lambda_phi` and `K_phi^(eff)` are shown to live only on the separate isotropic support-compensation branch.

## Stage 319 — exact microscopic -> reduced drift extractor

Files:
- `fivepn_stage319_322_common.py`
- `5pn_stage319_microscopic_reduced_drift_extractor.py`
- `5pn_stage319_microscopic_reduced_drift_extractor_output.txt`

The exact reduced Stage-318 tester inputs are

- `dln chi0`
- `dln deltaU`
- `dln ZhatW`
- `dln epsilonW`
- `dln epsilon_eta`

with

`chi0 = gamma c_etaU / K_U`,

`deltaU = pi^2 T_U / (L^2 K_U)`,

`ZhatW = lambda_W^2 mu_W / (K_eta^(eff) (K_W^(eff))^2)`,

`epsilonW = gamma^2 lambda_W^2 sigma / (K_U K_W^(eff))`,

`epsilon_eta = c_etaU^2 / (K_U K_eta^(eff))`.

The exact microscopic grouped-drift extractor is therefore

`dln chi0     = gamma1 + c1 - kappaU`,

`dln deltaU   = tau1 - kappaU`,

`dln ZhatW    = 2 lambda1 + mu1 - kappaEta - 2 kappaW`,

`dln epsilonW = 2 gamma1 + 2 lambda1 - kappaU - kappaW`,

`dln epsilon_eta = 2 c1 - kappaU - kappaEta`.

So the reduced five-drift packet is an exact linear image of the eight microscopic grouped drifts. No extra reduced closure assumptions are needed to pass from the microscopic branch to the Stage-318 tester.

## Stage 320 — exact microscopic actual-branch tester

Files:
- `5pn_stage320_microscopic_actual_branch_tester.py`
- `5pn_stage320_microscopic_actual_branch_tester_output.txt`

Composing the Stage-319 extractor with the Stage-317 reduced tester gives the direct microscopic monomial-drift packet

`Sigma_tr`
`= (1+deltaU_*) (gamma1 + c1 - kappaU) + (1+chi0_*) (tau1 - kappaU)`,

`Sigma_nt`
`= (2 lambda1 + mu1 - kappaEta - 2 kappaW)`
`  + E_* (2 gamma1 + 2 lambda1 - kappaU - kappaW)`
`  - F_* (tau1 - kappaU)`,

`Sigma_eta = 2 c1 - kappaU - kappaEta`.

So the physical defect packet is still the same triangular normal form,

`Theta_1 = - C_tr,* Sigma_tr`,

`Xi_1    =   A_tr,* Sigma_tr + Sigma_nt`,

`R_1     = - Xi_1 - epsilon_eta,* Sigma_eta/(1-epsilon_eta,*)`.

The direct microscopic tester has exact rank `3` inside the eight-drift kernel space, so its zero-defect kernel is `5`-dimensional.

A convenient exact solve is

`tau1 = kappaU - (1+deltaU_*) (gamma1 + c1 - kappaU)/(1+chi0_*)`,

`kappaEta = 2 c1 - kappaU`,

`mu1 = 2 c1 - kappaU + 2 kappaW - 2 lambda1`
`      - E_* (2 gamma1 + 2 lambda1 - kappaU - kappaW)`
`      - F_* (1+deltaU_*) (gamma1 + c1 - kappaU)/(1+chi0_*)`.

So the codimension-3 similarity-orbit closure from Stages 312–314 now appears directly inside the Stage-318 actual-branch tester language.

## Stage 321 — finite microscopic branch packet

Files:
- `5pn_stage321_finite_microscopic_branch_packet.py`
- `5pn_stage321_finite_microscopic_branch_packet_output.txt`

The finite Stage-318 quotient packet is now evaluated directly on the microscopic coherent-kernel state

`(lambda_W, c_etaU, gamma, K_U, K_eta^(eff), K_W^(eff), mu_W, T_U)`.

The three direct microscopic monomials are

`C_tr,*`
`= (gamma c_etaU / K_U)^(1+deltaU_*)`
`  (pi^2 T_U / (L^2 K_U))^(1+chi0_*)`,

`C_nt,*`
`= [lambda_W^2 mu_W / (K_eta^(eff) (K_W^(eff))^2)]`
`  [gamma^2 lambda_W^2 sigma / (K_U K_W^(eff))]^(E_*)`
`  [pi^2 T_U / (L^2 K_U)]^(-F_*)`,

`epsilon_eta = c_etaU^2 / (K_U K_eta^(eff))`.

And the finite quotient packet is exactly

`q_tr  = ln(C_tr,* / C_tr,*,ref)`,

`q_nt  = ln(C_nt,* / C_nt,*,ref)`,

`q_eta = ln(epsilon_eta / epsilon_eta,ref)`.

The direct physical branch observables are also explicit in microscopic variables:

`R_tr = (1 + chi0/(1+deltaU)) / (1 + chi0)`,

`R_target = Lambda_0 (1-epsilon_eta) (1-epsilon)^2 / [ ZhatW (1+chi0)^2 ]`,

with

`epsilon = epsilonW [ 1 - (2/11) deltaU/(1+deltaU) ]`.

So the actual moving-throat branch no longer needs a separate reduced-state postulate before the finite quotient packet can be compiled: it is already a direct log-ratio packet of three microscopic monomials and three direct microscopic branch observables.

## Stage 322 — exact support-blind orbit lock vs support-sensitive isotropic split

Files:
- `5pn_stage322_support_blind_orbit_split.py`
- `5pn_stage322_support_blind_orbit_split_output.txt`

The final bridge of this batch is the exact separation theorem.

Introduce the explicit support variables

- `lambda_phi`
- `K_phi^(eff)`

through the coherent support ratio

`zeta = lambda_phi^2 K_W^(eff) / (lambda_W^2 K_phi^(eff))`.

Then:

1. the orbit-lock observables

   `R_tr`, `R_target`, `epsilon_eta`

   are exactly independent of `lambda_phi` and `K_phi^(eff)`;

2. the finite quotient packet `(q_tr, q_nt, q_eta)` is also exactly independent of them;

3. but the isotropic support lane depends on them through the single enhancement ratio

   `zeta = lambda_phi^2 K_W^(eff) / (lambda_W^2 K_phi^(eff))`,

   with exact support-enhancement factor

   `S(zeta;epsilon) = 1 + zeta (1-epsilon)/(1-zeta epsilon)`.

The exact logarithmic support sensitivities are

`d ln zeta / d ln lambda_phi = 2`,

`d ln zeta / d ln K_phi^(eff) = -1`,

and therefore

`d M_tr / d ln lambda_phi > 0`,

`d M_tr / d ln K_phi^(eff) < 0`

on the physical branch.

So the weak-axisymmetric orbit-lock test and the isotropic support-compensation / normalization test are rigorously separated at the microscopic level.

## Bottom line after Stage 322

The continuation point has tightened again.

1. The weak-axisymmetric orbit-lock problem is now a direct tester on the microscopic drift vector of
   `(lambda_W, c_etaU, gamma, K_U, K_eta^(eff), K_W^(eff), mu_W, T_U)`.
2. Its exact zero-defect kernel is still the codimension-3 similarity-orbit closure already isolated in the notes.
3. The explicit support variables `(lambda_phi, K_phi^(eff))` do **not** enter that orbit-lock packet at all.
4. They enter only through the separate isotropic support-enhancement lane `zeta`.

So the next honest theorem gate is now even smaller:

> extract the actual microscopic branch drifts and finite kernel states from the completed moving-throat operator, feed them into the Stage-320/321 orbit packet, and test the support-enhancement lane separately through `zeta`.
