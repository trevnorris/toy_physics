# Atlas and graph migration checklist

Keep only the ideas below, and re-anchor them to original `research/` or
`software/` sources. Everything else in `atlas/` and `graph/` may be discarded.

## Charge and electromagnetism

- charge is not circulation;
- effective charge/thickness is not a breathing-mode identification;
- `kappa_rho` is not electric charge;
- the Maxwell gauge-localization patch is required where the source says so;
- the fluxoid identity and electric/magnetic ontology remain distinct.

Target: `memory/topics/charge-and-electromagnetism.md`.

## Projection and reduction

- projection onto the brane is not the same operation as dimensional
  reduction;
- a zero-mode reduction does not erase mixed bulk/brane fields from the parent
  ontology.

Target: `memory/topics/projection-and-reduction.md`.

## Moving-throat dynamics

- the parent confinement term does not by itself supply a strict wall PDE;
- wall coefficients may be branch data rather than universal constants;
- readout packets are not the throat itself;
- distinguish mouth, interior, support, exporter, material closure, and branch
  realization;
- retain the open actual-branch exporter and executable-branch-solver gaps.

Target: `memory/topics/moving-throat-dynamics.md`.

## Post-Newtonian and quadrupole results

- the local 4PN result is not automatically the complete hereditary/tail
  result;
- the 2.5PN identification remains conditional where stated;
- angular normalization does not settle radial or outgoing normalization;
- similarity-orbit closure is not a complete 5PN derivation.

Targets: `memory/topics/post-newtonian-ladder.md`,
`memory/topics/quadrupole-normalization.md`, and
`memory/topics/weak-axisymmetric-and-orbit-closure.md`.

## Reader pages

Retain concise coverage of foundations, projection/open-brane behavior,
charge/electromagnetism, the PN ladder, quadrupole radiation, moving-throat
dynamics, response packets, and status-reading rules.

Discard the Atlas-only future-paper queue and obsolete atom, lepton, anomaly,
and 5PN-planning material. Do not migrate content whose only support is root
`notes/`, root `docs/`, or a PDE ledger.

## Cutover

- [x] Every retained item above appears in a current topic page with original
      source citations.
- [x] `memory/index.md` points to the resulting pages.
- [x] A fresh agent can answer representative lookup questions using only the
      memory index and linked pages, without a source-tree grep.
- [ ] Delete `atlas/` and `graph/` in a separate, reviewable commit.
