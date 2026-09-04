---
title: Emergent electromagnetic construction
type: source
status: current
sources:
  - software/em_charge_attribute/reports/emergent_em_construction.md
  - software/em_charge_attribute/run_emergent_em_construction.sh
  - software/em_charge_attribute/emergent_em_sympy.py
  - software/em_charge_attribute/emergent_em_dual.wl
  - software/em_charge_attribute/compare_engines.py
last_updated: 2026-09-03
---

# Emergent electromagnetic construction

## Purpose and present standing

This package asks whether an added constrained internal degree of freedom can
support a compact U(1) sector with the desired electric/magnetic interaction
signs. It instantiates that degree of freedom with the link spins of a pinned
quantum-spin-ice model and identifies a physical `+w` or `-w` throat with a
flux-string endpoint of charge $Q=+1$ or $Q=-1$.

The strongest honest conclusion is a **consistency and identification result**,
not a derivation from the existing continuum throat equations. The report's
verification note says the composite throat operator is described but never
constructed by the code, deconfinement is imported from the cited spin-ice
literature, and several advertised controls are Boolean or tautological guards.

Sources:

- `software/em_charge_attribute/reports/emergent_em_construction.md` — opening verification note and heading `Honest scope and remaining obligations`
- `software/em_charge_attribute/emergent_em_sympy.py` — functions `microscopic_mapping`, `throat_embedding`, and `main`

## What the construction establishes conditionally

- Rewriting the microscopic ice constraint as an oriented divergence gives
  integer vertex charges $Q\in\{-2,-1,0,1,2\}$. A link flip creates an
  opposite-charge pair, a path product leaves charge only at its endpoints,
  and a closed ring move has zero divergence.
- Compact link flux and ring exchange supply the usual compact gauge
  redundancy. In the assumed Coulomb phase, the transverse projector has rank
  two and the quadratic IR Hamiltonian has two photon-like modes.
- The additive electric charge comes from the divergence sum and flux-string
  dressing, not from relabeling the geometric `+w/-w` sign as an additive
  charge.
- Eliminating one Maxwell field gives a repulsive like-density channel and an
  attractive parallel-current contribution with the same $1/k^2$ kernel.
  The deliberately substituted scalar mediator instead gives an attractive
  density channel and no transverse-current channel.

Sources:

- `software/em_charge_attribute/reports/emergent_em_construction.md` — headings `Constraint, defect charge, and compact gauge field`, `The ±w throat is an endpoint, not a renamed Z2 charge`, and `One Maxwell field and the dual sign`

## Evidence and limits

The SymPy and Wolfram programs independently check the finite incidence
algebra, hopping continuity identities, transverse projector, static kernels,
dimensions, and selected negative controls. `compare_engines.py` compares
their structured results. These calculations validate the encoded algebra;
they do not establish that the original throat model contains the proposed
spin-ice link variables or dynamically binds a throat to a dressed endpoint.

The make-or-break open question is therefore whether a microscopic parent of
the toy model actually realizes this independent constrained internal sector.
The package also leaves gravity/superfluid coexistence, 4+1D embedding, charge
normalization, and realistic spectra outside scope.

Sources:

- `software/em_charge_attribute/reports/emergent_em_construction.md` — headings `Dual-engine log` and `Honest scope and remaining obligations`
- `software/em_charge_attribute/emergent_em_dual.wl` — emitted checks `phase_existence = CITED_NOT_COMPUTED` and `single_kernel`

## Related memory

- `memory/sources/software/em-native-p.md` records why the existing native
  polar sector does not generically supply a Gauss constraint.
- `memory/sources/software/em-electric-sign.md` and
  `memory/sources/software/em-magnetism.md` test the far-field force sectors.
