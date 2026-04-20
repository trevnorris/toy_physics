# Moving-Throat PDE — Phase 0 Theorem-Target Sheet

## 1. Parent equations already frozen by the current stack

### Matter / GNLS sector
Use the exact bulk matter equation

i hbar D_t psi = [ -(hbar^2 / 2m) D_i D_i + V_conf(X; a, L) + h(rho) ] psi,

with
- rho = |psi|^2,
- P(rho) = K_EOS rho^5,
- h(rho) = dU/drho,
- exact bulk continuity
  partial_t rho + partial_i j^i = 0.

### Gauge sector
Use the exact localized 4+1 Maxwell equation

partial_M [ Z(w) F^{MN} ] + (1/xi) partial^N (partial · A) = mu_0 (J_psi^N + J_ext^N).

### Geometry sector in the current parent files
At present the geometry is only carried through effective throat DOFs
- a(t),
- L(t),

with baseline evolution
- M_a a_ddot + Gamma_a a_dot = - partial H_tot / partial a,
- M_L L_ddot + Gamma_L L_dot = - partial H_tot / partial L.

This is enough for the conservative reduced hierarchy, but not yet a true moving-throat PDE.

## 2. What must not be applied too early
Do not impose the far-field brane Maxwell reduction while building the PDE.
In particular, do not set to zero at the PDE level:
- A_w,
- partial_w A_mu,
- J^w,
- F_{mu w}.

Those are controlled brane-effective suppressions, not part of the microscopic ontology.

Do not replace the throat by the worldtube/particle limit.
That comes later, after the response operator has been derived.

## 3. Proposed geometry lift needed for a genuine PDE
This is the first new ingredient that is not yet frozen by the current files.

Introduce a distributed throat geometry variable, for example:
- a level-set Sigma(X,t)=0,
- or an equivalent moving wall / shape field.

Then promote the confinement potential from
- V_conf(X; a, L)

to something like
- V_conf(X; Sigma).

Design requirement:
- the old collective variables a(t), L(t) must be recoverable as the lowest moments of Sigma,
- grouped P2 and geometry channels must appear as explicit deformations of Sigma,
- the stationary reference throat must reduce to the already preferred finite geometry branch.

## 4. First linearization target
Choose a stationary finite-throat background
- psi_0,
- A_{M0},
- Sigma_0,

and linearize:
- psi = psi_0 + delta psi,
- A_M = A_{M0} + delta A_M,
- Sigma = Sigma_0 + delta Sigma.

The first theorem target is not the full nonlinear PDE. It is the linearized finite-throat response problem on the natural compact/passive/outgoing branch.

## 5. First observables the PDE must produce
The PDE must be asked to output the dynamic response data that the reduced hierarchy still leaves open:
- Omega_{1 perp},
- Omega_10,
- Omega_0,
- Omega_20,
- Omega_21,
- Omega_22,
- Omega_g.

It must also output the grouped low-frequency conservative data needed downstream:
- u2^(20),
- u2^(21),
- u2^(22),
- and, on the isotropic quadrupole branch, the canonical pair (Kbar_0, Kbar_2).

## 6. First drive/measure protocol
Use an operational mouth/worldtube response setup.

Effort variables:
- scalar / enthalpy perturbation,
- grouped real P2 drive,
- geometry drive.

Measured variables:
- normal mouth/worldtube flux,
- optional leakage flux in w,
- grouped response amplitudes.

Define the linear response operator
- j_i(omega) = sum_j Z_eff,ij(omega) u_j(omega).

## 7. Acceptance tests
The first acceptable PDE branch must pass all of the following:

### 2PN gate
It must reproduce the fact that the zero-frequency throat operator is fixed and that what remains open is the dynamic pole structure.

### 3PN gate
It must provide a microscopic origin for:
- the richer grouped real P2 constitutive law,
- the unique geometry completion.

### 2.5PN gate
It must stay on the compact/passive/outgoing branch and deliver the passive/outgoing STF quadrupole normalization.

### 4PN gate
Its quadrupole normalization must be compatible with the hereditary bridge target.

## 8. Working theorem statement for the next stage
Derive the linearized finite-throat Dirichlet-to-Neumann / mouth response operator on the natural compact passive outgoing branch, extract the grouped low-frequency data and dynamic pole scales, and test whether the resulting branch closes the grouped P2/geometry conservative lanes and the quadrupole normalization bridge.
