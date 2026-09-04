---
title: Quadrupole normalization
type: topic
status: current
sources:
  - research/4d_2_5pn/paper/4d_2_5pn.tex
  - research/4d_3pn/paper/4d_3pn.tex
  - research/4d_4pn/paper/4d_4pn.tex
  - research/pde/paper/pde.tex
  - research/pde_audit/notes/stage_v2_12_stf_angular_source_map_derivation.md
  - research/pde_audit/notes/stage_v2_26_program_status_after_audit.md
last_updated: 2026-09-03
---

# Quadrupole normalization

## Current bottom line

The angular/STF representation problem is closed for the chosen grouped real
$\ell=2$ basis. The physical radial/axial source strength and outgoing-wave
normalization are not. The same unresolved normalization controls both the
2.5PN radiation-reaction coefficient and the 4PN hereditary bridge.

## What is already fixed

The real $P_{20}\oplus P_{21}\oplus P_{22}$ ports span the full real STF
quadrupole representation. The audited canonical angular source map is the
identity in its normalized basis, so its angular normalization factor is one.
This rules out a missing spherical-harmonic convention as the source of the
current normalization gap.

It does **not** fix radial overlap, axial support, mouth/source amplitude,
worldtube coupling, or outgoing transfer.

Sources:

- `research/4d_2_5pn/paper/4d_2_5pn.tex` — label `app:quad-ledger-representation`
- `research/pde_audit/notes/stage_v2_12_stf_angular_source_map_derivation.md` — headings `Source-map theorem`, `Interpretation`, and `SymPy result`

## Conditional outgoing branch

Rotational invariance and a passive outgoing isotropic closure reduce the
quadrupole sector to one scalar response. Its first odd term has the required
$+i\omega^5$ parity, and the even/odd coefficients obey the canonical branch
relations quoted in the 2.5PN paper. In the moving-throat notation,

$$
\Gamma_5=\frac{\chi_QP_0a^5}{27c_s^5},\qquad P_0=\frac{N_0}{D_0},
$$

with $\chi_Q=1$ only on the canonical compact outgoing branch. These formulas
say what a successful branch must export; they do not prove that such a branch
exists.

Sources:

- `research/4d_2_5pn/paper/4d_2_5pn.tex` — labels `app:quad-ledger-outgoing`, `app:quad-ledger-matching`, and `eq:qnorm-canonical-invariants`
- `research/pde/paper/pde.tex` — labels `eq:outgoing-Gamma5`, `eq:outgoing-factorized-normalization`, and `eq:canonical-chiQ`

## The one normalization gate

Equivalent forms of the open target are

$$
\hat m_0^2\Gamma_5=\frac{2G}{5c^5}
$$

and, in the factorized moving-throat variables,

$$
\hat m_0^2\chi_QP_0=\frac{54Gc_s^5}{5a^5c^5}.
$$

Here $\hat m_0$ is the dimensionful source-scale carrier. It must not be
confused with the separate dimensionless far-field profile
$\hat m=1+O(a^2/r^2)$.

Sources:

- `research/4d_2_5pn/paper/4d_2_5pn.tex` — labels `sec:conditional-theorem-gap-normalization` and `app:qnorm-target`
- `research/pde/paper/pde.tex` — labels `eq:outgoing-BT-target`, `eq:outgoing-natural-source-map`, and `eq:outgoing-canonical-normalization`
- `research/4d_4pn/paper/4d_4pn.tex` — labels `eq:fixed-open-gamma-target` and `eq:fixed-open-P0-target`

## Relation to 3PN and 4PN

The 3PN grouped-$P_2$ construction tests conservative carrier structure and
isotropy. It explicitly leaves the outgoing 2.5PN normalization separate.
At 4PN, the hereditary coefficient is identified as
$GM/(2c^3)$ times the same $\gamma_{\rm quad}^{\rm eff}$, so there is no new
quadrupole-normalization constant inside the canonical bridge. A model-side
transport factor such as $\Theta_{\rm tail}$ must still take its canonical
value or be derived equivalently.

Sources:

- `research/4d_3pn/paper/4d_3pn.tex` — labels `sec:interface-isotropy` and `sec:interface-open`
- `research/4d_4pn/paper/4d_4pn.tex` — labels `eq:4pn-tail-exact-bridge`, `eq:4pn-tail-toy-transport`, and `eq:4pn-tail-theta-one`

## Present branch evidence

The audit finds a coherent target surface, but current frozen reduced and
manufactured branch families miss it, while the physical nonlinear exporter is
not complete. That is compatible with the clean angular and algebraic results:
the former test representation and interface formulas; the failed candidates
test whether particular branch data land on them.

Source:

- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — headings `Why the clean algebra and failed simulations are not contradictory`, `What the audit says is missing`, and `Current status statement`

Related pages:

- `memory/topics/post-newtonian-ladder.md`
- `memory/topics/moving-throat-dynamics.md`
