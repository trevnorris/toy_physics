---
title: Paper Em Fields
type: source
status: current
sources:
- research/em_fields/paper/em_fields.tex
last_updated: '2026-09-02'
---

## Purpose and scope

### source-paper-em-fields--purpose-and-scope — Emergent electromagnetic sector

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper extends a compressible-superfluid defect toy model by assigning effective electromagnetic potentials, charge, and forces to throat-like defects. Its claimed scope is a classical proof-of-principle treatment in weak-field, long-wavelength, quasi-static, and nearly flat regimes, rather than a complete microscopic or empirically fitted theory.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{sec:intro}`
- `research/em_fields/paper/em_fields.tex` — `\label{sec:discussion}`

## Source-unit map

- Entry point and sole member: `research/em_fields/paper/em_fields.tex`
- Role: `primary`
- Read mode: `semantic`
- Unit shape: `file`

## Key statements

### source-paper-em-fields--hydrodynamic-em-dictionary — Potentials and homogeneous Maxwell identities

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The paper defines \(\phi_{\mathrm{EM}}=\lambda(h+v^2/2)\) and \(\mathbf A=\lambda\mathbf v\), giving \(\mathbf B=\lambda\nabla\times\mathbf v=\lambda\boldsymbol{\omega}\) and an electric field proportional to the negative Euler-acceleration combination. With the usual potential definitions, \(\nabla\cdot\mathbf B=0\) and Faraday’s law follow as vector identities for sufficiently smooth fields, without using the fluid equations of motion.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{eq:EM-potentials-def}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:EB-def}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:B-vorticity}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:E-euler}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:divB-zero}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:Faraday-identity}`

### source-paper-em-fields--inhomogeneous-maxwell — Ampère–Maxwell equation and continuity

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

In the weak-field, nearly flat treatment, the paper assumes the sourced acoustic wave equation \(\Box A^\mu=-\mu_0J^\mu\) and the Lorenz gauge. From these premises it derives the Ampère–Maxwell equation; taking the four-divergence of the wave equation under the gauge condition yields current conservation. This result does not independently fix the Coulomb or Gauss-law normalization.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{eq:wave-A}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:lorenz-gauge}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:ampere-maxwell}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:continuity}`

### source-paper-em-fields--breathing-mode-coulomb-field — Exterior radial solution

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a static, spherically symmetric configuration sufficiently far from the throat for the acoustic metric to be approximately flat, the paper replaces the acoustic wave operator by the flat-space Poisson operator. Outside the core it obtains the radial Laplace solution \(A+B/r\); removing the additive gauge constant gives a \(1/r\) potential, a \(1/r^2\) radial electric field, and radius-independent flux. This is an exterior point-source treatment, not a derivation of finite-core microphysics or of the coefficient \(B\).

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{subsec:gauss-breathing}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:poisson-phi}`
- `research/em_fields/paper/em_fields.tex` — `\label{app:breathing-gauss}`

### source-paper-em-fields--gauss-law-matching — Gauss-law normalization closure

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

At the effective-theory level, the paper fixes the exterior solution’s coefficient by demanding that its flux equal \(q/\epsilon_0\). The resulting Coulomb normalization and differential Gauss-law form use \(\epsilon_0=1/(\mu_0c^2)\) and are matching closures, not independent derivations of the normalization from the fluid equations or resolved core dynamics.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{eq:gauss-flux}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:gauss-differential}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:app-gauss-integral}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:eps-mu-relation}`

### source-paper-em-fields--charge-mass-hierarchy — Geometric mass, charge, and force hierarchy

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper adopts \(m_G=\kappa_m\rho_0\pi a^2L\), with \(\kappa_m\) packaging unresolved model dependence, and defines \(q=\kappa_q\rho_0\pi a^2\Gamma\). Its cavity-enthalpy extremization gives the preferred aspect ratio \(L/a=\sqrt{2}\pi/x_{01}\simeq1.85\), but a nonzero-radius extremum exists in this model only for negative \(P_{\mathrm{vac}}\). At fixed aspect ratio, the resulting identical-defect force ratio scales as \(F_{\mathrm{elec}}/F_{\mathrm{grav}}\propto\Gamma^2/a^2\). These are phenomenological estimates and normalization choices, not proofs of microscopic particle properties, measured force hierarchies, or SI-normalized couplings.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{eq:mG-def}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:q-def}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:aspect-ratio}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:bvac-condition}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:force-ratio-scaling}`
- `research/em_fields/paper/em_fields.tex` — `\label{app:cavity}`
- `research/em_fields/paper/em_fields.tex` — `\label{app:units}`

### source-paper-em-fields--magnus-lorentz-match — Straight-core magnetic-force correspondence

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a straight vortex throat with approximately uniform force along its length, the paper decomposes the Magnus force and matches its core-velocity-dependent term to the magnetic Lorentz force per unit length. The equality requires the coarse-grained field choice \(B_0=-L/(\kappa_q\pi a^2)\).

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{eq:Magnus-decomp}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:fL-mag}`
- `research/em_fields/paper/em_fields.tex` — `\label{app:lorentz-magnus}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:app-B0-solution}`

### source-paper-em-fields--electric-force-closure — Schematic electric-force identification

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

For a defect small compared with the background variation scale, the paper schematically combines pressure, enthalpy, kinetic, and unsteady-flow contributions into an electric-like force. Identifying \(q=\rho_0V_{\mathrm{eff}}/\lambda\) then gives \(F_E=qE\); geometric factors, near-core coupling, breathing-mode excitation, and detailed normalization remain part of the matching closure.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{sec:lorentz-magnus}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:F-press}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:q-from-volume}`
- `research/em_fields/paper/em_fields.tex` — `\label{eq:F-electric}`

## Computed evidence represented by the source

### source-paper-em-fields--computed-evidence-not-supplied — Computed evidence availability

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit. The paper reports that companion Mathematica and SymPy scripts were used for internal-consistency checks, while describing the work as exploratory and still in need of subject-matter review; the scripts, invocations, literal outputs, and separate interpretations needed for measured evidence are absent from this source unit.

Sources:

- `research/em_fields/paper/em_fields.tex` — heading `\section*{Acknowledgments}`

## Assumptions, exclusions, and open questions

### source-paper-em-fields--bulk-magnetic-sector-open — Long-range magnetic sector absent

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

Strict bulk irrotationality confines effective vorticity and magnetic field to defect cores. The paper therefore does not provide Biot–Savart-type long-range magnetostatics or transverse electromagnetic waves; it leaves a controlled vorticity- or shear-supporting extension, including a proposed brane–bulk topology, for future work.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{sec:discussion}`

### source-paper-em-fields--microphysics-and-phenomenology-open — Present-regime exclusions and unresolved extensions

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The present construction treats electromagnetism linearly in a weak-field, quasi-static background; idealizes defects as cylindrical cavities and point sources; applies its long-wavelength results only at distances large compared with \(a\) and \(L\); and neglects radiation reaction, self-force, gravitational radiation, and strong-field effects. It also does not quantize cavity or field modes, derive photons or realistic particle spectra, or demonstrate that one equation of state and defect geometry reproduce measured \(G\), \(\epsilon_0\), and \(\mu_0\). Nonlinear and dynamical extensions, resolved core microphysics, quantization, realistic spectra, and detailed phenomenological fitting remain open.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{sec:discussion}`
- `research/em_fields/paper/em_fields.tex` — `\label{app:units}`

## Revision and supersession relationships

### source-paper-em-fields--relationship-to-papers-i-iii — Extension of the prior model

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper presents its electromagnetic construction as an extension that assumes the gravitational results of Papers I–III rather than rederiving them. It reports no formal retraction or whole-source supersession of those earlier treatments; their gravitational sector is used as background for the electromagnetic analogue.

Sources:

- `research/em_fields/paper/em_fields.tex` — `\label{sec:intro}`
- `research/em_fields/paper/em_fields.tex` — `\label{sec:discussion}`

## Related topics and scripts

- `memory/topics/charge-and-electromagnetism.md`
- `memory/conflicts.md`
