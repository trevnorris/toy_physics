---
title: Charge and electromagnetism
type: topic
status: current
sources:
  - research/4d/paper/4d.tex
  - research/4d_em_fields/paper/4d_em_fields.tex
  - research/em_fields/paper/em_fields.tex
  - research/4d_plasma/paper/4d_plasma.tex
  - research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md
  - software/em_charge_attribute/reports/native_p_constraint_gate.md
  - software/em_charge_attribute/reports/emergent_em_construction.md
  - software/em_charge_attribute/puncture_deflection_electric_sign_result.md
  - software/em_charge_attribute/magnetism_moving_throat_result.md
  - software/em_charge_attribute/reports/u1_body_dynamics.md
last_updated: 2026-09-03
---

# Charge and electromagnetism

## Current bottom line

The model currently **assumes a localized Maxwell sector and adopts a
consistent charge dictionary**. It has not derived electromagnetism from the
native throat or superfluid variables. The native polar-field route fails its
generic quadratic Gauss-constraint test; an added spin-ice degree of freedom
offers a possible realization but is only identified with the throat, not
dynamically derived. Electric boundary selection, moving-throat source
normalization, full body mechanics, and radiation remain open.

## Current charge dictionary

For current reading, use

$$
\eta_Q=\pm1,\qquad q_\star=\eta_Q e_\star,\qquad
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}}.
$$

- $\eta_Q$ is the oriented puncture/throat branch sign.
- $q_\star$ is the fixed charge assigned to the microscopic branch.
- $q_{\rm eff}$ is the brane coupling after Maxwell zero-mode normalization.
- circulation and fluxoid winding belong to the magnetic/vortical sector;
  they do not define electric charge.
- throat breathing changes geometry and support energy, not $\eta_Q$ or
  $q_\star$.

The ontology is internally consistent but still adopted: the model does not
derive the magnitude $e_\star$ or an additive microscopic charge carrier from
the throat equations.

Sources:

- `research/4d/paper/4d.tex` — labels `sec:fields`, `eq:qeff-canonical`, `eq:vorticity-gauge`, and `eq:fluxoid-quantization`
- `research/4d_em_fields/paper/4d_em_fields.tex` — labels `sec:canonical_charge` and `eq:qeff_canonical`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — headings `Electric charge variables`, `Circulation and magnetism`, and `Minimal claim boundary`

## Scoped correction to the older EM paper

The older `em_fields.tex` used circulation, throat area/geometry, and breathing
to motivate several effective definitions of electric charge. Those
definitions are superseded for the **charge ontology** by the
$\eta_Q,q_\star,q_{\rm eff}$ dictionary above. This is not a whole-paper
retraction: its potential identities, exterior $1/r$ solution, and conditional
Magnus/Lorentz analogy remain useful within their assumptions.

Similarly, the 1PN/2PN response coefficients concern gravitational/internal
response. The historical gravity-side symbol $q$ is now $\kappa_\rho$, and the
2PN residual coefficients $q_i$ are ADM bookkeeping, not electric charge.

Sources:

- `research/em_fields/paper/em_fields.tex` — labels `eq:q-def`, `subsec:gauss-breathing`, `eq:q-from-volume`, and `sec:discussion`
- `research/4d_em_fields/paper/4d_em_fields.tex` — label `sec:intro`
- `research/4d_1pn_full/paper/4d_1pn_full.tex` — label `sec:intro`
- `research/4d_2pn/paper/4d_2pn.tex` — labels `sec:intro-motivation` and `sec:carry-forward`

## What the localized Maxwell sector provides

The parent action yields localized Maxwell equations and current-consistency
conditions. Under the controlled zero-mode ansatz it produces a brane Maxwell
equation with thickness-normalized coupling. Gaussian localization supplies a
massless mode plus a gapped transverse tower, giving Coulomb behavior with
Yukawa corrections. Mixed $A_w,F_{\mu w},J^w$ channels survive in the parent
theory even when suppressed by the brane far-field reduction.

This establishes a workable effective Maxwell carrier because it was placed
in the parent action. It does not explain why that carrier emerges from the
superfluid or close the exterior circulation/intake plumbing of an actual
throat.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — labels `eq:integrated_exact`, `eq:brane_maxwell_exact_weighted`, `eq:zero_mode_ansatz`, `eq:gaussian-spectrum`, and `eq:coulomb-yukawa`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — headings `What EM structure is present` and `What is reduced versus still open`

## Candidate microscopic carriers

### Native polar field

The quadratic Hamiltonian/Dirac search finds no first-class Gauss constraint
for either native $P^a$ theory on the regular kinetic stratum. Tuned
rank-drop/common-null loci do not yield a nonzero gradient Gauss descendant in
the tested analysis. Because the non-common tuned loci were sampled rather
than exhaustively classified, the result is a generic quadratic no-Gauss
statement, not a universal no-go theorem.

Source:

- `software/em_charge_attribute/reports/native_p_constraint_gate.md` — headings `Completeness and stratum scope`, `THEORY-A: computed H₂, Dirac chain, and G search`, `THEORY-C: computed H₂, Dirac chain, and G search`, and `Decision-table result`

### Added spin-ice degree of freedom

The spin-ice construction shows how an **additional** constrained link degree
of freedom can carry integer divergence charge, a compact U(1) field, two
transverse modes, and the repulsive-density/attractive-parallel-current sign
pair. It identifies `+w/-w` with flux-dressed $Q=\pm1$ endpoints. The code does
not construct or dynamically bind the composite throat operator, and the
deconfined phase is cited rather than computed.

Source:

- `software/em_charge_attribute/reports/emergent_em_construction.md` — opening verification note and headings `Constraint, defect charge, and compact gauge field`, `The ±w throat is an endpoint, not a renamed Z2 charge`, and `Honest scope and remaining obligations`

## Force and dynamics status

- **Electric sign:** the coupled kernel and $1/R^2$ falloff are fixed, but the
  mouth boundary ensemble is not. Held value is positive, fixed source is
  negative, fixed monopole is indefinite, and mixed conditions span outcomes.
  Current landing: `R1_REQUIRED(bc_selection)`.
- **Magnetism:** two calculations agree on the conditional Darwin tensor,
  falloff, and velocity order. The coefficient still depends on the unresolved
  electric branch and throat-current normalization. The candidate magnetic
  field is time-reversal even, a characterized departure from Maxwell.
- **U1 body dynamics:** Phase A is conditionally well posed, but Phase B1 body
  mechanics and Phase B2 return/intake/radiative quantities remain unresolved;
  Phase C is not run.

Sources:

- `software/em_charge_attribute/puncture_deflection_electric_sign_result.md` — body items `Q-BC`, `Q-FORCE / Q-COMBINE`, and `§4 landing`
- `software/em_charge_attribute/magnetism_moving_throat_result.md` — body items `Q-DIRECT`, `Q-COMPARE / Q-MAG`, and `SEALED §4 landing`, plus appendix `Hooks`
- `software/em_charge_attribute/reports/u1_body_dynamics.md` — headings `Verdict and Phase-A halt`, `Phase B1 — indexed mechanics remediation 3`, and `Phase B2 — Intake response and radiative residues`

## Current answer to “has EM emerged?”

No. The repository has an assumed localized Maxwell action, a corrected charge
ontology, several controlled reductions, and useful candidate-mechanism tests.
It does not yet have a single realized throat branch that derives the charge
carrier, selects the electric boundary class, fixes the moving source, and
closes conservative plus radiative dynamics.

Related pages:

- `memory/topics/projection-and-reduction.md`
- `memory/scripts/em-charge-attribute.md`
- `memory/conflicts.md`
