---
title: Paper 1Pn Spin And Nbody
type: source
status: current
sources:
- research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex
last_updated: '2026-08-31'
---

## Purpose and scope

### source-paper-1pn-spin-and-nbody--purpose-and-regime — Conditional 1PN completion of the toy model

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper extends a phenomenological model of sink-like defects in a compressible superfluid to spin-induced gravitomagnetism and conservative \(N\)-body dynamics at first post-Newtonian order. Subject to the existence of its assumed long-range translational wake, it claims that the three-paper construction reproduces GR’s scalar, optical, spin, and vector 1PN observables with \(\beta_{\mathrm{PPN}}=\gamma_{\mathrm{PPN}}=1\). This observational matching does not establish a unique bare metric: the model permits distinct acoustic and hydrodynamic projections for light and matter, calibrated to agree at 1PN. Its regime is weak-field and slow-motion, excluding a strong-field, radiative, electromagnetic, or cosmological completion.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:intro-motivation}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:intro-roadmap}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:discussion-summary}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:discussion-limitations}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{app:eih-metric}`

## Source-unit map

- Entrypoint and sole member: `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: monolithic TeX paper

## Key statements

### source-paper-1pn-spin-and-nbody--dyon-spin-calibration — Phenomenological dyon flow is calibrated to Kerr

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper models a spinning defect as a sink bound to a localized vortex ring and imposes, rather than derives from microscopic hydrodynamics, a far-zone physical azimuthal speed proportional to \(D\sin\theta/r^2\). Using its acoustic-metric dictionary, it matches the resulting off-diagonal metric to weak-field Kerr and fixes \(D=4GJ/c^2\), up to circulation and spin-orientation sign conventions. The underlying proportionality and microscopic or topological origin of this flow remain unspecified.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:spin-dyon}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:vphi-dyon}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:acoustic-metric}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:D-calibration}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{app:dyon-flow}`

### source-paper-1pn-spin-and-nbody--lense-thirring-observables — Calibrated metric is asserted to reproduce spin observables

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper states that inserting the already calibrated dyon-induced \(g_{0i}\) into the standard 1PN precession formula recovers the Lense–Thirring vector. For a nearly circular orbit of radius \(r\) around a central dyon whose spin \(J\) is aligned with the \(z\)-axis, it gives the nodal-precession rate \(2GJ/(c^2r^3)\) to leading order in \(J\). Once \(J\) is specified, the construction provides no independent continuous adjustment of the gravitomagnetic strength.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:Omega-LT-GR}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:spin-observables}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:spin-falsifiability}`

### source-paper-1pn-spin-and-nbody--spin-flow-not-translational-wake — Rotational and translational flows are distinct

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The paper separates the rotation-driven dyon flow used for spin observables from the translation-driven wakes used for the EIH velocity cross terms. A line vortex has the wrong radial and angular behavior for Lense–Thirring precession, while overlap of rigid-sphere dipolar backflows decays too quickly to generate the required \(1/r\) EIH interaction.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:spin-mismatch}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:backflow}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-vector}`

### source-paper-1pn-spin-and-nbody--static-eih-sector — Density-dependent mass supplies the static three-body structure

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper expands a density-dependent defect mass in the local potential and obtains \(G^2/c^2\) terms with the characteristic three-body distance products of the static EIH sector. It asserts that the inherited calibration \(\kappa_\rho=1\) gives the exact EIH coefficient, but the cited Mathematica analysis and its literal output are not present in this prepared unit, so the exact-coefficient claim remains provisional here.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-static}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:mass-density-dep}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{app:static-derivation}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:app-V3-EIH-form}`

### source-paper-1pn-spin-and-nbody--wake-tensor-match — Isotropic wake decomposition matches the EIH cross tensor

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Assuming a translational wake with Fourier amplitude proportional to \(1/k\), the paper evaluates its overlap using transverse, longitudinal, and optional helical projector components. Matching the resulting coefficients to the EIH cross-term targets \((-7/2,-1/2)\) gives \(\alpha^2=\tfrac34(1-a_H^2)\); the minimal choice \(a_H=0\) fixes \(\alpha^2=3/4\) and \(K=2/\pi^2\).

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:wake-ansatz}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-derivation}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:Cpara-wake}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:CL-wake}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:alpha-squared}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{eq:K-wake}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{app:alpha-tuning}`

### source-paper-1pn-spin-and-nbody--wake-stability-boundary — Real tuning removes the earlier imaginary coupling

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper attributes an earlier negative \(\alpha^2\) to an incomplete transverse wake basis and states that the real value \(\alpha^2=3/4\) makes the quadratic wake functional positive-definite in the usual Euclidean hydrodynamic sense. This stability interpretation is not accompanied in the prepared unit by a separate verification record.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-derivation}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-signature}`

## Computed evidence represented by the source

### source-paper-1pn-spin-and-nbody--computed-evidence-boundary — Verification chain is absent

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit. The paper mentions an accompanying Mathematica script for symbolic coefficient checks, but the prepared unit contains neither that script nor an exact invocation, readable literal output, and separate interpretation chain. Consequently, this capsule does not classify any result as measured.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-derivation}`

## Assumptions, exclusions, and open questions

### source-paper-1pn-spin-and-nbody--long-range-wake-gate — Existence of the required wake remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The EIH cross-term construction is conditional on moving sink defects producing a long-range translational wake with \(u(\mathbf{k})\propto1/k\), equivalently \(u(\mathbf{x})\propto1/r^2\), so that wake overlap produces a \(1/r\) interaction. Whether the proposed stiff superfluid naturally generates this response is explicitly deferred to PDE-level simulation.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-signature}`

### source-paper-1pn-spin-and-nbody--microscopic-origin-open — Microscopic defect dynamics are unresolved

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The microscopic or topological origin of the imposed dyon far field is deferred. The paper also leaves open whether throat boundary conditions or a more microscopic vacuum model naturally produce the selected thermodynamic and wake parameters rather than requiring them as phenomenological inputs.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:spin-dyon}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:discussion-limitations}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{app:dyon-flow}`

### source-paper-1pn-spin-and-nbody--excluded-physics — Beyond-1PN sectors are not supplied

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The construction is restricted to conservative weak-field, slow-motion dynamics. Strong-field behavior, gravitational radiation and radiation reaction, electromagnetic modes, and embedding in an expanding cosmological background are explicitly unbuilt or unexplored.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:discussion-limitations}`

## Revision and supersession relationships

### source-paper-1pn-spin-and-nbody--series-continuation — Continuation with a scoped wake-basis correction

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper presents itself as a continuation of earlier orbital and optical treatments and does not declare that it supersedes or retracts those works. It does explicitly correct the earlier interpretation requiring \(\alpha^2<0\): within the translational-wake tensor match, that result is attributed to an incomplete transverse basis and replaced by a full isotropic-projector treatment with real parameters.

Sources:

- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:intro-summary-papers12}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:intro-roadmap}`
- `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` — `\label{subsec:nbody-derivation}`

## Related topics and scripts

No related memory pages or script domains were supplied by the task.