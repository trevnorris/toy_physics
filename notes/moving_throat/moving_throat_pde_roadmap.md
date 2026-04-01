# Moving-Throat PDE Roadmap

## Bottom-line goal
Do **not** jump straight to a fully nonlinear moving-throat PDE theorem.
The correct next target is a **linearized finite-throat passive/outgoing response problem** whose low-frequency data reproduce the already frozen conservative and dissipative ledgers.

## Why this is the right target
The current stack already says:
- the parent 4D action is exact at the bulk level,
- the open problem after 2PN is the **dynamic throat-response completion**,
- the open problem after 3PN is the microscopic derivation of the **grouped real P2 constitutive law** and the **geometry completion**,
- the open problem after 2.5PN and 4PN is the **passive/outgoing STF quadrupole normalization**.

So the missing PDE must be asked to do one very specific job:
produce the dynamic pole data, the grouped P2/geometry response, and the quadrupole normalization on the natural branch.

## Phase 0 — Freeze the theorem target
Target theorem:

Given a stationary finite throat background, derive the **linearized moving-throat Dirichlet-to-Neumann / response operator** on the natural compact, passive, outgoing branch, and extract:
- scalar/geometry response,
- grouped real P2 response,
- dynamic pole scales,
- and the quadrupole normalization needed for 2.5PN and 4PN.

Acceptance outputs:
- dynamic observables
  - Omega_{1\perp}, Omega_{10}, Omega_0, Omega_{20}, Omega_{21}, Omega_{22}, Omega_g,
- grouped low-frequency coefficients
  - u2^(20), u2^(21), u2^(22),
- isotropic quadrupole data
  - (Kbar_0, Kbar_2),
- odd quadrupole normalization
  - gamma_quad^eff.

## Phase 1 — Promote geometry from collective coordinates to a distributed throat variable
Current parent files only carry geometry through effective DOFs a(t), L(t).
To define a genuine moving-throat PDE, replace that finite-dimensional closure by a distributed geometry representation, e.g. one of:
- a level-set field Sigma(X,t)=0,
- a moving wall/shape field R(theta,phi,z,t),
- or an equivalent boundary embedding.

Design rule:
- a(t) and L(t) must reappear as the lowest collective moments of the new geometry field,
- the grouped P2 and geometry lanes must appear as explicit shape/response channels rather than as symbolic residuals.

## Phase 2 — Write the finite-throat linearized bulk problem
Start from the exact parent bulk system:
- GNLS matter sector,
- localized 4+1 Maxwell sector,
- moving geometry/confinement sector.

Then linearize about a stationary finite throat background:
- psi = psi_0 + delta psi,
- A_M = A_{M0} + delta A_M,
- Sigma = Sigma_0 + delta Sigma.

Keep mixed channels alive:
- A_w,
- J^w,
- F_{\mu w},
- E_w,
- C_a.

Do **not** impose the far-field zero-mode suppression at this stage.
That limit is for the brane-effective reduction, not for the PDE that is supposed to generate the response data.

## Phase 3 — Define the mouth/worldtube response operator
Choose a drive/measure protocol at the mouth/worldtube and compute the linear response matrix

j_i(omega) = sum_j Z_eff,ij(omega) u_j(omega).

Recommended first port basis:
- scalar / breathing port,
- grouped real P2 ports,
- geometry port.

The low-frequency target is a stable expansion away from resonances:
Z_eff(omega) ~ Z^(0) + i omega Z^(1) + omega^2 Z^(2) + ...

The immediate output of this phase is the dynamic pole data and grouped low-frequency coefficients.

## Phase 4 — Match back to the frozen PN ledger
Use the extracted low-frequency data to check three already-frozen gates:
1. 2PN: the dynamic pole structure and zero-frequency packaging,
2. 3PN: the richer grouped real P2 closure and the geometry completion,
3. 2.5PN / 4PN: the passive/outgoing quadrupole normalization.

This is the decisive falsification stage.
If the extracted data miss any of these gates, the PDE/branch choice is wrong.

## Phase 5 — Close scalar and dipole side conditions from the same PDE
Do not restart the old generic scalar/dipole doom stories.
Instead, use the new PDE to answer the narrow remaining questions:
- whether projection-locking or an equivalent suppression mechanism holds for the scalar soft modes,
- whether the dipole outgoing branch is indeed cubic and finite-size demoted,
- whether the completed PDE stays on the compact passive/outgoing branch already isolated by the 2.5PN audit.

## Phase 6 — Close the quadrupole normalization
Once the same PDE gives the natural isotropic passive/outgoing quadrupole branch, extract the canonical pair:
- Kbar_0,
- Kbar_2.

Then compute gamma_quad^eff and test the target
- gamma_quad^eff = 2G / (5 c^5).

If this closes, then:
- the normalized 2.5PN odd quadrupole coefficient closes,
- the normalized 4PN hereditary coefficient closes,
- and the present conditional 4PN theorem becomes unconditional inside the same hierarchy.

## The first step to do now
Start with a **Ground-Truth / Theorem-Target sheet** containing:
1. the exact parent equations we are allowed to use,
2. the precise branch assumptions,
3. the list of observables the PDE must output,
4. the first drive/measure port basis,
5. the acceptance tests against 2PN / 3PN / 2.5PN / 4PN.

That document is the minimum scaffolding needed before any trustworthy PDE work begins.
