---
title: Paper 1Pn Optics
type: source
status: current
sources:
- research/1pn_optics/paper/1pn_optics.tex
last_updated: '2026-08-31'
---

## Purpose and scope

### source-paper-1pn-optics--optical-clock-extension — Optical and clock extension

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper extends a superfluid-defect toy model from its previously assumed orbital sector to weak-field, 1PN gravitational optics and clock effects. It treats light and massive defects as probing distinct optical and hydrodynamically dressed projections of one fluid-based construction; strong fields, non-spherical or time-dependent sources, and complete electromagnetic or cosmological treatments remain outside its scope.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:intro-roadmap}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:discussion-future}`

## Source-unit map

- Entry point: `research/1pn_optics/paper/1pn_optics.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: monolithic TeX paper (`file`)

## Key statements

### source-paper-1pn-optics--newtonian-orbital-calibration — Newtonian potential and calibrated kinetic dressing

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For the assumed effective orbital Lagrangian, the paper takes the static scalar-sector potential to be exactly Newtonian at 1PN and derives a perihelion coefficient \(2\beta\) from the kinetic prefactor \(\sigma(r)=\beta GM/(c^2r)\). Matching the GR coefficient \(6\) fixes this model parameter to \(\beta=3\).

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:paper1-L-eff}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:paper1-Phi-eff}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:paper1-precession-total}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:paper1-beta-3}`

### source-paper-1pn-optics--vacuum-profile — Pressure, density, and refractive-index profile

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Assuming a barotropic polytrope, weak perturbations, spherical quasi-static hydrostatic balance, and the constitutive Newtonian potential \(\Phi=-GM/r\), the paper derives \(\Delta P=-GM\rho_0/r\), \(\Delta\rho/\rho_0=-GM/(c_0^2r)\), and
\[
N(r)\simeq1+\frac{n-1}{2}\frac{GM}{c_0^2r}.
\]
For \(n=5\), this becomes \(N(r)\simeq1+2GM/(c_0^2r)\).

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:hydrostatic-balance}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:deltaP}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:deltarho}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:N-general-n}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:N-n5}`

### source-paper-1pn-optics--light-bending — Weak-field light bending

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Using geometric optics, a straight unperturbed path, and the \(n=5\) index profile, the paper evaluates the transverse-index-gradient integral and obtains
\[
\Delta\theta=\frac{4GM}{bc^2}.
\]
Comparison with \(2(1+\gamma)GM/(bc^2)\) gives the observational optical value \(\gamma_{\mathrm{obs}}=1\).

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:lensing-deflection-general}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:lensing-deflection-integral}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:lensing-dtheta-result}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:lensing-ppn}`

### source-paper-1pn-optics--shapiro-delay — Weak-field Shapiro delay

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For endpoints far from the mass and a straight-path approximation, the same index profile yields the one-way delay
\[
\Delta t\simeq\frac{2GM}{c^3}
\ln\!\left(\frac{4r_{\mathrm E}r_{\mathrm R}}{b^2}\right),
\]
matching the leading GR expression and again corresponding to \(\gamma_{\mathrm{obs}}=1\).

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:shapiro-t-integral}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:shapiro-dt-integral}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:shapiro-I-exact}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:shapiro-dt-result}`

### source-paper-1pn-optics--redshift-scaling — Density-based clock redshift

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Under the additional model assumptions that a defect’s volume and mode structure remain fixed, \(m\propto\rho\), and the clock frequency obeys \(\omega\propto m\), the density deficit gives
\[
\frac{\delta\omega}{\omega_0}
=\frac{\delta m}{m_0}
=\frac{\delta\rho}{\rho_0}
=-\frac{GM}{rc^2}+O\!\left(\frac{G^2M^2}{c^4r^2}\right).
\]
This establishes the paper’s GR-like redshift only for clocks satisfying those scaling assumptions.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:redshift-deltarho}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:redshift-dm}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:redshift-omega-m}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:redshift-full-chain}`

### source-paper-1pn-optics--clock-universality — Clock universality is conditional

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For generalized scalings \(m\propto\rho^q\) and \(L\propto\rho^s\), equal density dependence of atomic and cavity-clock frequencies requires
\[
s=\frac{n-1}{2}-q.
\]
The calibrated \(n=5,\ q=1\) branch therefore requires \(L\propto\rho\); universal GR-like redshift across clock types is a microphysical consistency condition, not an automatic consequence of the optical profile.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:clock-universality-s}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:rod-scaling-n5}`

### source-paper-1pn-optics--bimetric-boundary — Optical and soliton metrics must not be conflated

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper distinguishes the pure optical metric, whose \(g_{tt}^{(\mathrm{opt})}=-1\) packages the bending effect into its spatial sector, from the hydrodynamically dressed projection used for massive solitons. Treating defects as test particles of the optical metric would give a perihelion coefficient \(10\) instead of \(6\); the paper attributes the corrected orbital result to the separately calibrated kinetic dressing.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:lensing-ppn}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:ppn-equivalence}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:trilemma}`

### source-paper-1pn-optics--n5-selection — Restricted \(n=5\) selection claim

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Within the paper’s restricted family of spherically symmetric flux-tube defects, barotropic polytropes, a shared Newtonian limit, and pure refraction, matching the GR bending coefficient selects \(n=5\). The paper prefers this branch over tuned drag or split \(n=3\) constructions; this is not a uniqueness result beyond that model class.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:degeneracies-menu}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:degeneracies-n5}`

## Computed evidence represented by the source

### source-paper-1pn-optics--computed-evidence-absent — Paper-cited numerical checks are absent from the unit

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper cites companion implementations at `scripts/ray_tracing.py` and `scripts/shapiro_delay.py`, describes numerical ratios close to unity, and provides intended optics algorithms and convergence procedures. Those scripts and their literal outputs are not members of the prepared unit, so the required recorded command, readable output, and separate interpretation chain is **not supplied in the prepared unit**; the numerical claims cannot be classified as measured here.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:lensing-numerics}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:shapiro-numerics}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{app:numerics-optics}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:app-convergence}`

## Assumptions, exclusions, and open questions

### source-paper-1pn-optics--microphysics-open — Pressure-volume microphysics remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper obtains \(\kappa_\rho=1\) from the static density-deficit contribution and derives the spherical added-mass value \(\kappa_{\mathrm{add}}=1/2\). In contrast, it treats the pressure-volume breathing coefficient \(\kappa_{\mathrm{PV}}\) phenomenologically, fixes it to \(3/2\) through orbital calibration so that \(\beta=3\), and explicitly defers its first-principles microphysical derivation.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:deltarho}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{eq:paper1-beta-decomp}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{app:dipole-derivation}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:paper1-beta}`

## Revision and supersession relationships

### source-paper-1pn-optics--supersession-not-established — Extension without formal supersession

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The source presents itself as an extension of the earlier orbital treatment and calls its optical-versus-GR discrepancy calculation a revised version of the discussion in its own main text. It does not explicitly declare that the earlier paper or any complete prior source unit is corrected, retracted, or superseded; no formal supersession relationship is therefore established by the prepared source.

Sources:

- `research/1pn_optics/paper/1pn_optics.tex` — `\label{subsec:intro-roadmap}`
- `research/1pn_optics/paper/1pn_optics.tex` — `\label{app:ppn-metric}`

## Related topics and scripts

No related memory pages or script domains were supplied by the task.