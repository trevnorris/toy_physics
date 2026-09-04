---
title: Foundations
type: topic
status: current
sources:
  - research/4d/paper/4d.tex
  - research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
  - research/pde/paper/pde.tex
last_updated: 2026-09-03
---

# Foundations

## Working picture

The project is a weak-field, slow-motion toy universe built in four spatial
dimensions plus time. Its declared parent theory contains a confined gauged
nonlinear Schrödinger matter field and a localized Maxwell sector. A
three-dimensional brane observer sees projected fields near $w=0$, while a
defect is modeled as a finite throat joining the brane mouth to a bulk
interior.

This is an organizing analogy, not a claim about our universe and not yet a
microscopic material theory. The papers mix exact consequences of declared
actions with controlled reductions and phenomenological closures, so the
status of each statement matters as much as its formula.

Sources:

- `research/4d/paper/4d.tex` — labels `sec:action` and `sec:discussion-limitations`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — labels `subsec:intro_ontology` and `subsec:intro_scope`
- `research/pde/paper/pde.tex` — labels `sec:introduction` and `sec:status-parent`

## Geometry and objects

The older variables $a(t)$ and $L(t)$ describe mouth radius and throat depth.
They remain useful collective coordinates, but the moving-throat framework
promotes the fundamental geometry to the distributed level set
$\Sigma(\mathbf X,t)=r-R(\Omega,w,t)$. Real $\ell=2$ wall/support harmonics
then become physical response channels rather than merely names for PN
coefficient lanes.

Keep the following objects separate:

- the **mouth**, where brane observables and boundary data are defined;
- the **finite interior/support**, which supplies axial modes and endpoint
  conditions;
- the **material and gauge fields**, which respond to the geometry;
- collective moments such as $a,L$;
- reduced readouts such as $D_0$, $P_0$, and $N_n$, which summarize a chosen
  response calculation and are not themselves pieces of the throat.

Sources:

- `research/pde/paper/pde.tex` — labels `eq:geometry-lift-levelset`, `eq:geometry-lift-a-moment`, `eq:geometry-lift-L-moment`, and `eq:geometry-lift-harmonic-decomp`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — labels `subsec:throat_modes` and `eq:app_mode_spectrum`

## Four layers of reasoning

1. **Parent theory.** Vary the declared bulk action and use its exact
   identities: GNLS/Maxwell equations, current continuity, and mixed
   gauge-invariant observables.
2. **Observation.** Project bulk quantities with an explicitly normalized
   observer kernel. This can leave exchange or leakage with the bulk.
3. **Controlled reduction.** Impose a zero-mode, quasi-static, small-body,
   stable-mode, or PN regime and derive the corresponding effective operator.
4. **Closure and realization.** Choose constitutive/boundary data, solve an
   actual branch, and only then compare its exported observables with a target.

Skipping from layer 1 to layer 4 is the most common source of overclaiming in
the older material.

Sources:

- `research/4d/paper/4d.tex` — labels `sec:braneobs`, `sec:emreduction`, and `sec:poisson-regime`
- `research/pde/paper/pde.tex` — labels `sec:linearized-pde`, `sec:reduced-system`, `tab:claim-status`, and `sec:discussion-gap`

## Current foundation boundary

The action-level equations and specified projection identities are the firmest
parts of the program. Geometry dynamics, material coefficients, boundary
selection, mode truncation, and physical branch realization remain adopted or
open at different stages. The memory therefore treats successful algebraic
assembly as evidence about a declared toy model, not as evidence that nature
or even a completed nonlinear throat realizes it.

Related pages:

- `memory/topics/projection-and-reduction.md`
- `memory/topics/moving-throat-dynamics.md`
- `memory/topics/status-and-reading-rules.md`
