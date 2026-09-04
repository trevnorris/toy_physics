---
title: Moving-throat magnetism calculation
type: source
status: current
sources:
  - software/em_charge_attribute/magnetism_moving_throat_result.md
  - software/em_charge_attribute/magnetism_moving_throat_check.py
  - software/em_charge_attribute/magnetism_moving_throat_check.wl
last_updated: 2026-09-03
---

# Moving-throat magnetism calculation

## Purpose and result

This target-blind Tier-A calculation asks whether slowly translating signed
throats source a magnetic-like transverse interaction. It adds a scoped brane
transverse-displacement row, derives the convection-like source
$\mathbf J_T=q_Ts\eta_a\mathbf V$, and compares two routes to the
velocity-dependent pair force.

Both routes give the same Darwin tensor, $1/R$ interaction potential,
$1/R^2$ force, and quadratic velocity order. For side-by-side motion the
velocity term reverses between parallel and antiparallel velocities. But the
coefficient comparison remains conditional because the nonlinear throat
normalization $q_T$ and the electric coefficient $A_E$ are both open.

The production landing is therefore `R1_REQUIRED(electric_bc_selection)`, with
additional open blockers for the direct moving-throat source, magnitude, and
route consistency. On a separately selected like-repelling electric branch,
the direct parallel-current term would be attractive; this package does not
select that upstream branch.

Sources:

- `software/em_charge_attribute/magnetism_moving_throat_result.md` — heading `Concise result (body; ≤2 pages)`, especially `Q-CURRENT`, `Q-BOOST`, `Q-DIRECT`, `Q-COMPARE / Q-MAG`, and `SEALED §4 landing`

## Two calculation routes

- **Route A** starts from a Lorentz-completed conditional electric anchor and
  derives the general two-velocity Darwin correction through
  $O(v^2/c_\gamma^2)$.
- **Route B** starts from the translated throat source and transverse Euler
  equation, reconstructing the transverse Green tensor without consuming the
  Route-A result.
- Their coefficient ratio is
  $r_{BA}=q_T^2/(\rho_{\rm br}A_E)$, while cone agreement requires
  $c_E^2\rho_{\rm br}/\mu_R=1$. Neither equality is derived here.

Sources:

- `software/em_charge_attribute/magnetism_moving_throat_result.md` — headings `Independent derivations and dimensions` and `Parallel, antiparallel, and zero controls`
- `software/em_charge_attribute/magnetism_moving_throat_check.py` — direct-source, route-independence, force, and comparison checks

## Limitations and characterized departure

The finite throat profile, active-drain time arrow, and positive transverse
kinetic/stiffness coefficients are declared ingredients. The proposed magnetic
field $\nabla\times\mathbf u_T$ is time-reversal even in this construction,
whereas Maxwell magnetic field is time-reversal odd; the report records this
as a characterized departure. Higher velocity orders, possible contamination
by other operators, cone locking, and the non-conservative active-flux force
remain open, so full-force integrability and emergent Lorentz invariance are
not established.

Sources:

- `software/em_charge_attribute/magnetism_moving_throat_result.md` — appendix heading `Hooks` and final `DECIDED`/`R1` paragraph

## Related memory

- `memory/sources/software/em-electric-sign.md` owns the unresolved electric
  boundary anchor used here.
