---
title: Puncture-deflection electric-sign calculation
type: source
status: current
sources:
  - software/em_charge_attribute/puncture_deflection_electric_sign_result.md
  - software/em_charge_attribute/puncture_deflection_electric_sign_independent_verdict.md
  - software/em_charge_attribute/puncture_deflection_electric_sign_check.py
  - software/em_charge_attribute/puncture_deflection_electric_sign_check.wl
  - software/em_charge_attribute/puncture_deflection_electric_sign_independent_recompute.py
last_updated: 2026-09-03
---

# Puncture-deflection electric-sign calculation

## Purpose and result

This is a target-blind Tier-A calculation of the far-field interaction between
two signed throat deflections at $R\gg r_e$. It identifies the geometric
field as the normal embedding displacement $\xi_w=\ell h$, reduces the live
parent coupling to the dimensionless $h$ channel, derives the coupled static
kernel, and evaluates the interaction under several possible mouth boundary
ensembles.

The calculation does **not** earn a unique electric sign. The existing model
does not select the boundary class of a signed mouth, and the admissible
ensembles disagree:

- held value (`V`) gives a positive pair coefficient;
- fixed source (`J`) gives a negative coefficient;
- fixed monopole (`M`) is sign-indefinite;
- the mixed ensemble can be negative, zero, or positive over its allowed
  parameters.

All have a $1/R$ potential and $1/R^2$ force at this order. Because the
boundary-class ambiguity comes first, the report's production landing is
`R1_REQUIRED(bc_selection)`, not a sign success or a model inconsistency.

Sources:

- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` — heading `Concise result (body; ≤2 pages)`, especially `Q-BC`, `Q-FORCE / Q-COMBINE`, and `§4 landing`
- `software/em_charge_attribute/puncture_deflection_electric_sign_independent_verdict.md` — complete independent recompute summary

## What is decided

Within the declared quadratic reduction, the package fixes the field identity,
the positive coupled kernel for $D=B_{\rm eff}K_h-C_{hu}^2>0$, the ensemble
Legendre transforms, the conditional `REPLACE` and `ADD` coefficients, and the
$1/R^2$ force law. An independent Python recomputation confirms the central
kernel, determinant, ensemble signs, mixed range, and falloff.

The Python and Wolfram production paths also enumerate the sealed decision
table and exercise mutations intended to catch wrong functionals, injected
landings, dimension errors, and broken controls. A successful run establishes
agreement about the implemented Tier-A calculation, not physical selection of
one ensemble.

Sources:

- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` — appendix headings `Field, action, kernel, and ledger`, `Full ensemble definitions`, and `Production checks`
- `software/em_charge_attribute/puncture_deflection_electric_sign_check.py` — ensemble and sealed-adjudicator implementation
- `software/em_charge_attribute/puncture_deflection_electric_sign_independent_recompute.py` — independent kernel and sign recomputation

## Open inputs and limitations

The missing nonlinear parent-throat/core functional must determine whether the
mouth holds a value, source, conormal flux, or mixed relation. `REPLACE` versus
`ADD`, mixed parameters, the deflection normalization, and the existence of a
signed radial monopole also remain unresolved. The magnitude contains free
core factors, the calculation is far-field only, and it supplies no local
density prediction or close-range result.

Sources:

- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` — `Q-BC`, `Q-MAG / hooks`, and the final `DECIDED`/`R1-deferred` paragraph

## Related memory

- `memory/sources/software/em-magnetism.md` uses this electric coefficient as
  a deliberately conditional anchor.
