---
title: Research memory index
type: index
status: current
sources: []
last_updated: 2026-09-03
---

# Research memory

This is the entry point for the toy-model research. Start here, open the one
topic page closest to the question, and follow its citations only when the
technical details matter.

## Current picture

The repository has a declared matter-plus-localized-Maxwell parent action,
controlled projection/reduction machinery, and conditional conservative PN
matches through local 4PN. It also isolates an outgoing quadrupole route for
2.5PN reaction and the 4PN hereditary term. These are not yet a complete
physical derivation: the nonlinear moving-throat/material branch, EM carrier
and electric boundary selection, and outgoing normalization remain open.
Tested frozen candidate branches miss their target packets without invalidating
the exact identities or every possible completion.

## Where to look

| Question | Start here |
|---|---|
| What is the parent model, and what is fundamental versus a readout? | [Foundations](topics/foundations.md) |
| What is projected, reduced, localized, or leaked through the brane? | [Projection and reduction](topics/projection-and-reduction.md) |
| What is the current charge dictionary? Has EM actually emerged? | [Charge and electromagnetism](topics/charge-and-electromagnetism.md) |
| Which PN orders are covered, and under what assumptions? | [Post-Newtonian ladder](topics/post-newtonian-ladder.md) |
| What remains open in 2.5PN radiation and the 4PN tail? | [Quadrupole normalization](topics/quadrupole-normalization.md) |
| What is the moving-throat action/branch status? | [Moving-throat dynamics](topics/moving-throat-dynamics.md) |
| What do the weak-axisymmetric tests and orbit-closure criterion establish? | [Weak-axisymmetric and orbit closure](topics/weak-axisymmetric-and-orbit-closure.md) |
| How strong is a claim, and what does a script pass mean? | [Status and reading rules](topics/status-and-reading-rules.md) |
| Which older claims were corrected, and which choices remain unresolved? | [Conflicts](conflicts.md) |
| What is the current future research program and its decision sequence? | [Lean research proposal](sources/research-proposal.md) |

## Script catalogs

- [PDE audit](scripts/pde-audit.md) — symbolic and numerical audit entry points,
  frozen-branch protocol, and claim matrix.
- [EM and charge checks](scripts/em-charge-attribute.md) — emergent-carrier,
  native-$P$, electric-sign, magnetism, and U1 body-dynamics checks.
- [Force visualizer](scripts/force-visualizer.md) — renders the configured
  comparison laws and produces a verification report; it does not derive them.
- [Stage-1 solver](scripts/stage1-solver.md) — effective-closure branch,
  convergence, freeze, and verdict workflow.

## High-priority open gates

- Select and derive a nonlinear moving-throat/material action with consistent
  source, return, boundary, and stability laws.
- Derive an EM carrier from the throat variables—or clearly retain localized
  Maxwell as an independent sector—and select the electric mouth ensemble.
- Compute the source-to-outgoing quadrupole normalization on one physical,
  frozen branch; 2.5PN and the 4PN tail share this gate.
- Export conservative, weak-axisymmetric, orbit-closure, intake/return, and
  radiative observables from that same branch without post-target refitting.
- Treat the cataloged PN ladder as ending at 4PN; there is no permitted 5PN
  result in this memory.

## Source capsules

Individual summaries live in:

- [paper capsules](sources/papers/)
- [software and result capsules](sources/software/)
- [maintained PDE-audit capsule](sources/pde-audit.md)
- [lean research-proposal capsule](sources/research-proposal.md)

The original files under `research/` and `software/` are authoritative. Root
`notes/`, root `docs/`, and all `research/pde_ledger*` trees are deliberately
outside this memory. Source units are curated reading bundles, not complete
transitive build-dependency graphs.

To check freshness, run `python3 memory/update.py status`. To understand the
small maintenance workflow, read [memory/README.md](README.md).
