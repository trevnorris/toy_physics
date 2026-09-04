---
title: Post-Newtonian ladder
type: topic
status: current
sources:
  - research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
  - research/4d_1pn_full/paper/4d_1pn_full.tex
  - research/4d_2pn/paper/4d_2pn.tex
  - research/4d_2_5pn/paper/4d_2_5pn.tex
  - research/4d_3pn/paper/4d_3pn.tex
  - research/4d_4pn/paper/4d_4pn.tex
last_updated: 2026-09-03
---

# Post-Newtonian ladder

## How to read the ladder

The PN papers build a coherent algebraic ladder after a hierarchy of
projection, small-body, response, constitutive, and chart choices is fixed.
They are not successive numerical solutions of the parent moving-throat PDE.
At every rung, distinguish exact coefficient/Legendre-transform algebra from
the closure assumptions that supplied the coefficients.

There is no cataloged 5PN rung in this memory.

## 1PN: conservative EIH assembly

The bridge paper derives individual topology, optical, thermodynamic, and
kinematic contributions. EIH cross-tensor matching initially leaves a family;
the declared $n=5$ thermodynamic relation selects
$\alpha^2=3/4$, $a_H=0$, and $K_{\rm vec}=2/\pi^2$. Its
$\kappa_{\rm PV}=3/2$ breathing response is likewise a one-degree-of-freedom
adiabatic closure rather than a general throat response.

The full 1PN paper is the best single source for the assembled ledger. Once
the mass dressing, added mass, breathing closure, wake tensor, and static
nonlinear ansatz are accepted, the result is algebraically identical to the
conservative EIH Lagrangian and gives the usual $6\pi$ test-orbit precession.
The exact assembly does not make the input closures action-derived.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — labels `eq:family_alphaK`, `sec:1pn-interaction-thermo`, `eq:alphaK_final`, and `sec:kappaPV:adiabatic1dof`
- `research/4d_1pn_full/paper/4d_1pn_full.tex` — labels `sec:taxonomy`, `sec:kappaPV-closure`, `sec:wake-unique`, `eq:derived-eih`, and `eq:eih-residual-zero`

## 2PN: conservative generic-frame match

The 2PN construction adds a protocol-level low-frequency denominator repair,
a self/static seed, and a unique comparable-mass residual solve. With those
choices and the lower-order ledger frozen, its perturbative Legendre transform
matches the standard generic-frame ADM Hamiltonian through $c^{-4}$ with no
remaining dimensionless coefficient in that chosen conservative basis.

Dynamic throat poles, absolute normalization, damping, non-adiabatic response,
finite-size structure outside the strict gate, and radiation are not supplied
by that algebraic match.

Sources:

- `research/4d_2pn/paper/4d_2pn.tex` — labels `sec:one-body-ansatz`, `clm:adm-solve`, `eq:full-assembly-exact-match`, and `sec:status-open`

## 2.5PN: first dissipative gate

A real, even-frequency conservative kernel cannot generate local 2.5PN
radiation reaction; an odd retarded sector is required. The paper reduces the
Burke–Thorne force to the Iyer–Will member $(\alpha,\beta)=(4,5)$ and shows that
the real grouped $\ell=2$ ports have the necessary STF representation.

Under a compact, passive, outgoing, isotropic, strict-small-body closure, the
remaining scalar/dipole/internal-quadrupole channels are demoted and the
quadrupole response has the right tensor, parity, and frequency structure.
The physical normalization of that outgoing branch remains open.

Sources:

- `research/4d_2_5pn/paper/4d_2_5pn.tex` — labels `sec:benchmark-framework-conservative-nogo`, `sec:benchmark-framework-bt`, `app:quad-ledger-representation`, `app:quad-ledger-outgoing`, `sec:conclusion-proven`, and `sec:conditional-theorem-gap-normalization`

## 3PN: conservative three-lane assembly

The 3PN paper works in a fixed ADM chart. It derives a cubic Legendre compiler,
shows that Hamiltonian-first reduction is required, enlarges the grouped
$P_2$ middle block after a smaller family fails, and adds a geometry-side
static completion. Once these protocol-closure ingredients are inserted, the
three-lane result matches the target Hamiltonian algebraically.

Lower-order viability data may not be retuned to repair a 3PN mismatch. The
actual branch must still realize conservative isotropy and the separate
outgoing normalization; the conservative construction does not close the
2.5PN gate.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — labels `sec:compmass-hamiltonian-first`, `sec:grouped-richer-family`, `sec:final-main`, `sec:dict-carry`, and `sec:interface-open`

## 4PN: local reconstruction plus hereditary bridge

The 4PN paper derives a quartic chart compiler and emphasizes the same
Hamiltonian-first firewall: extra ordinary-chart powers of the symmetric mass
ratio can be translation effects rather than new physical carriers. Within
its closure hierarchy it reports a complete local instantaneous
reconstruction.

The hereditary tail coefficient is tied to the **same** effective quadrupole
normalization used at 2.5PN,

$$
C_{\rm tail}=\frac{GM}{2c^3}\gamma_{\rm quad}^{\rm eff}.
$$

Thus the canonical bridge introduces no second quadrupole-normalization datum,
but it does not solve the original normalization or transport gate.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — labels `eq:4pn-legendre-h4`, `eq:4pn-local-hamiltonian-first-scheme`, `eq:4pn-local-final-theorem`, `eq:4pn-tail-exact-bridge`, `sec:tail-no-new-datum`, and `sec:fixed-open-still-open`

## Current status

- Conservative 1PN, 2PN, 3PN, and the local 4PN sector are assembly results
  inside increasingly rich declared closures.
- 2.5PN isolates the necessary retarded STF quadrupole route but leaves its
  physical outgoing normalization open.
- The 4PN hereditary sector inherits that same open normalization and an
  associated transport condition.
- None of these algebraic successes replaces an actual branch solve or permits
  lower-order coefficients to be refit after seeing a higher-order target.

Related pages:

- `memory/topics/quadrupole-normalization.md`
- `memory/topics/moving-throat-dynamics.md`
- `memory/topics/status-and-reading-rules.md`
