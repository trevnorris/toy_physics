---
schema_version: 2
id: source-paper-1pn-hybrid
title: Paper 1Pn Hybrid
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/1pn_hybrid/paper/1pn_hybrid.tex
content_owner: ai_generated
last_updated: '2026-08-31'
generated_from_commit: 3f7234d3ba673817864c8938be209d451da19f36
source_kind: paper
source_unit:
  id: paper-1pn-hybrid
  shape: file
  entrypoint: research/1pn_hybrid/paper/1pn_hybrid.tex
  unit_digest_sha256: 0e5bed36ee13dea3cf10197fb6916eb99b4e067e7f826f718acf54a9e4331709
  members:
  - path: research/1pn_hybrid/paper/1pn_hybrid.tex
    role: primary
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 65227f6ca9f5b64b858a8f59b32a5f9b27e5aa53
    blob_size: 136440
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-paper-1pn-hybrid--purpose-and-scope — Hybrid 1PN and strong-field program

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper studies a toy model combining a barotropic polytropic vacuum \(P\propto\rho^n\), density-dependent defect inertia \(M\propto\rho^q\), and flow, refractive-index, wake, and finite-throat sectors. Exact coefficient comparisons with GR are restricted to first post-Newtonian order. The paper treats possible 2PN deviations qualitatively, while using analytic arguments and reduced phenomenological models to explore acoustic horizons, photon-sphere-like structures, and finite-throat effects rather than supplying a complete nonlinear theory.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:intro-goals}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{sec:discussion}`

## Source-unit map

- Entrypoint: `research/1pn_hybrid/paper/1pn_hybrid.tex`
- Role: primary monolithic paper
- Read mode: semantic
- Git object: blob `65227f6ca9f5b64b858a8f59b32a5f9b27e5aa53` (`136440` bytes, mode `100644`)

## Key statements

### source-paper-1pn-hybrid--linear-mass-density-scaling — Newtonian matching fixes \(q=1\)

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Given the paper’s normalization of \(\Phi\) and its assumed weak-field density response, expansion of the variable-mass Lagrangian produces the Newtonian potential term \(-M_0q\Phi\). Matching its coefficient to \(-M_0\Phi\) fixes \(q=1\), so defect mass tracks density linearly within the stated model class.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:L-scalar-expansion}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:q-equals-one}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:match-newtonian-app}`

### source-paper-1pn-hybrid--weak-field-optical-index — Polytrope controls the weak-field index

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Expanding \(c_s^2=Kn\rho^{n-1}\) with the assumed linear density response gives \(N(\Phi)=1-\tfrac12(n-1)\Phi/c^2+\mathcal O(\Phi^2/c^4)\). At \(n=5\) and \(\Phi=-GM/r\), the paper obtains \(N(r)\simeq1+2GM/(rc^2)\); this result is confined to the weak-field expansion.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:N-of-Phi}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:N-of-r-n5}`

### source-paper-1pn-hybrid--n5-uniqueness-claim — Claimed 1PN selection of \(n=5\)

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper reports that matching the hybrid scalar, kinematic, and flow-dependent coefficients to the full EIH structures has the unique real solution \((n,q)=(5,1)\). The prepared paper displays the scalar and kinematic dependence but leaves the decisive flow functions \(C_i^{(\mathrm{flow})}(n)\) schematic, so the uniqueness claim cannot be checked from this unit alone.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:n5-q1-solution}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:matching-eq-n}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:C1-tot}`

### source-paper-1pn-hybrid--wake-calibration — Wake ratio is inherited rather than rederived

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The vector sector adopts \(\alpha^2=3/4\) as an inherited calibration representing the longitudinal-to-transverse wake-power balance required by the cited earlier EIH treatment. This paper explicitly does not rederive the wake kernel.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:alpha-condition}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:alpha-relation}`

### source-paper-1pn-hybrid--transonic-condition — Regular sonic-point condition

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For steady, spherically symmetric, barotropic inflow, the paper combines continuity and the Euler equation into a first-order ODE. A smooth solution at its critical point must satisfy \(u^2(r_H)=c_s^2(r_H)\) together with vanishing of the effective driving term.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:inflow-ode-basic}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:transonic-conditions-app}`

### source-paper-1pn-hybrid--acoustic-horizon-ontology — Emergent horizon on a smooth throat

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper defines the sphere where \(|u(r_H)|=c_s(r_H)\) as an acoustic horizon for linearized perturbations on the three-dimensional brane. It states that this horizon is not a fundamental singularity and that the underlying four-dimensional bulk and throat geometry remain smooth across it.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:acoustic-horizon}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:acoustic-horizon-def}`

### source-paper-1pn-hybrid--horizon-scaling — Horizon and finite-throat scaling are phenomenological

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper reports representative reduced profiles with \(r_H\gtrsim a\) and approximately \(r_H\sim\kappa GM/c^2\), where \(\kappa\) depends on throat radius and imposed flux. Its standalone inflow and impedance constructions are phenomenological models rather than complete solutions of the throat microphysics. The finite-throat treatment uses \(\Lambda=L/a\) and records the inherited Paper V calibration \(\Lambda_\star\simeq1.85\).

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:rH-a-relationship}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:rH-scaling}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:acoustic-inflow}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:inflow-numerics}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:impedance_scaling}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:aspect-ratio}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:intro-params}`

### source-paper-1pn-hybrid--photon-sphere-estimates — Photon-sphere and shadow estimates

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

For the assumed optical metric, circular-ray conditions depend on extrema of \(N^4(r)r^2\). The paper’s weak-field-matched, finite-size profiles yield photon-sphere-like radii and critical impact parameters proportional to \(GM/c^2\), with model-dependent coefficients. It approximates the distant-observer shadow radius by \(b_{\mathrm{sh}}\approx b_{\mathrm{ph}}\) and warns that a quantitatively reliable photon sphere requires a nonperturbative near-throat optical profile.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:photon-sphere-N}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:rph-beta}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:rph-scaling}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{eq:bph-scaling}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:photon-sphere-shadow}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:lensing-photon-shadow}`

## Computed evidence represented by the source

### source-paper-1pn-hybrid--computed-evidence-boundary — Verification chain unavailable

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit. The paper names symbolic, inflow, impedance, scaling, and lensing scripts and says its calculations were checked for internal consistency with companion scripts, but this single-file unit supplies no exact recorded invocation, readable literal transcript containing the relevant operands and residual, or separate interpretation record. The manuscript characterizes itself as exploratory and as benefiting from subject-matter review; accordingly, its script-dependent uniqueness and numerical-profile claims are not classified as measured.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:1pn-matching}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:inflow-numerics}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:impedance_scaling}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:predictions-La}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:lensing-photon-shadow}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\section*{Acknowledgments}`

## Assumptions, exclusions, and open questions

### source-paper-1pn-hybrid--open-validation-gates — Nonlinear and higher-PN validation remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper explicitly leaves unresolved a full 2PN matching and radiation-reaction theory, fully dynamical three- or four-dimensional superfluid simulations, precise acoustic-horizon and photon-sphere scaling, throat stability under time-dependent perturbations, and strong-field electromagnetic coupling. It identifies higher-PN calculations and throat-resolving simulations as required to determine these objects.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{sec:discussion}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:predictions-tests}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{app:lensing-limitations}`

## Revision and supersession relationships

### source-paper-1pn-hybrid--series-relationship — Earlier-paper results are treated as inputs

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The manuscript presents itself as the sixth paper in a series and treats several earlier-paper values as calibrated inputs, including \(\beta=3\), \(n=5\), \(\alpha^2=3/4\), and \(\Lambda_\star\simeq1.85\). It does not explicitly declare that it corrects, retracts, or supersedes Papers I–V, so this prepared unit establishes no formal supersession relationship.

Sources:

- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:intro-background}`
- `research/1pn_hybrid/paper/1pn_hybrid.tex` — `\label{subsec:intro-params}`

## Related topics and scripts

No related memory pages or script domains were supplied by the transaction.