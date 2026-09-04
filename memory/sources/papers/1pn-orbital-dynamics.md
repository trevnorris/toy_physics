---
title: Paper 1Pn Orbital Dynamics
type: source
status: current
sources:
- research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex
last_updated: '2026-08-31'
---

## Purpose and scope

### source-paper-1pn-orbital-dynamics--purpose-and-regime — Superfluid-defect 0PN and 1PN model

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper proposes a scalar superfluid-defect toy model for Newtonian gravity and leading perihelion precession in the weak-field, slow-motion, test-mass central-field regime. It separates an instantaneous Poisson constraint potential from a finite-speed lag potential. Vortical and electromagnetic-like components are turned off in the gravitational field sector, while entering indirectly in the paper’s interpretation of hydrodynamic inertia.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:intro}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:model}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_beta}`

## Source-unit map

- Entry point and sole member: `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex`
- Role: `primary`
- Read mode: `semantic`

## Key statements

### source-paper-1pn-orbital-dynamics--static-limit-theorem — Static sources reduce to the Poisson sector

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a strictly time-independent density, late-time relaxation, boundary conditions sending the lag potential to zero at spatial infinity, and regularity away from the source, the paper derives that the lag field solves the homogeneous Laplace equation and is removable up to an irrelevant constant. The total potential therefore equals the Newtonian Poisson solution; for a point sink it is \(\Phi=-\mu/r\) with acceleration \(-\mu/r^2\).

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:static_limit_theorem}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:newton_potential_mu}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:newton_accel}`

### source-paper-1pn-orbital-dynamics--scalar-null-precession — Static scalar lag gives no 1PN precession

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

In the rest frame of a static central defect and in the negligible-test-mass limit, the paper specializes the retarded scalar solution. Vanishing source velocity reduces the Liénard–Wiechert denominator to unity, so the retarded potential equals \(-\mu/r\), the lag residual vanishes, and the scalar-only orbit remains a closed Kepler ellipse with \(\Delta\varphi_{\mathrm{scalar}}=0\).

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:retarded_scalar_1pn}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:scalar_lw_general}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:scalar_static_phi}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:Delta_phi_scalar_result}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{app:precession_scalar_only}`

### source-paper-1pn-orbital-dynamics--inertia-driven-precession — Kinetic prefactor produces a general 1PN advance

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Keeping the scalar potential Newtonian and adopting the kinetic-prefactor ansatz \(\sigma(r)=\beta\mu/(c_s^2r)\), the paper expands the effective Lagrangian to first order in \(\mu/(c_s^2r)\) and derives
\[
\Delta\varphi_{\mathrm{tot}}
=(2\beta)\frac{\pi\mu}{c_s^2a(1-e^2)}.
\]

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:sigma_ansatz}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{app:precession_with_sigma}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{app:precession_beta_result}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:Delta_phi_total_beta}`

### source-paper-1pn-orbital-dynamics--beta-gr-matching — Schwarzschild matching fixes beta phenomenologically

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper imposes agreement with the Schwarzschild test-body perihelion advance after identifying \(\mu=GM\) and \(c_s=c\). Comparing the modeled coefficient \(2\beta\) with the Schwarzschild coefficient \(6\) fixes \(\beta=3\) as a phenomenological matching condition, not as a hydrodynamic derivation.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:beta}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:Delta_phi_total_beta}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_kappa_PV}`

### source-paper-1pn-orbital-dynamics--hydrodynamic-coefficients — Proposed decomposition of beta

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper interprets
\[
\beta=\kappa_\rho+\kappa_{\mathrm{add}}+\kappa_{\mathrm{PV}}
\]
as contributions from density-dependent cavitation mass, entrained-fluid added mass, and compressible pressure–volume inertia. It assigns \(\kappa_\rho=1\) and \(\kappa_{\mathrm{add}}=1/2\) using simplified hydrodynamic arguments, but the companion notebooks cited for additional support are not members of the prepared unit.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{eq:beta_decomp}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_kappa_rho}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_kappa_add}`

### source-paper-1pn-orbital-dynamics--reported-numerical-checks — Paper-reported numerical behavior

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper reports grid and reduced-orbit experiments in which the static lag residual decays toward truncation error, scalar-only orbits show null precession, and reduced orbits with \(\beta=3\) approach the modeled GR-like 1PN advance under time-step and resolution refinement. These are source assertions only: the prepared unit contains neither the executable notebooks nor a recorded command/output/interpretation chain.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:numerics}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:static_numeric}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{app:numerics}`

## Computed evidence represented by the source

### source-paper-1pn-orbital-dynamics--computed-evidence-not-supplied — Computed evidence chain absent from prepared unit

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:numerics}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{app:numerics}`

## Assumptions, exclusions, and open questions

### source-paper-1pn-orbital-dynamics--applicability-boundaries — Restricted effective regime

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The construction assumes weak fields and slow motion, a flux-neutral local frame, a static heavy central source with a negligible-mass test body for its central result, and late-time and spatial-infinity conditions that remove the relaxed lag field. Uniform-drift invariance is postulated for the conservative 1PN treatment rather than demonstrated by a supplied verification record. Comparable-mass dynamics, full PPN behavior, radiation reaction, and vortical or spin sectors remain outside the demonstrated scope.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:udi}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:flux_neutrality}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:static_limit_theorem}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:retarded_scalar_1pn}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:point_sink_static}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:discussion_future}`

### source-paper-1pn-orbital-dynamics--pressure-volume-open-problem — Pressure–volume coefficient is not derived

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper does not derive \(\kappa_{\mathrm{PV}}\) from throat microphysics. It sets \(\kappa_{\mathrm{PV}}=3/2\) as the remainder required by the phenomenological value \(\beta=3\) after the proposed density and added-mass contributions. Whether a concrete superfluid equation of state and throat geometry produce this coefficient is explicitly left for future work.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_kappa_PV}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_status}`
- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:discussion_future}`

### source-paper-1pn-orbital-dynamics--exploratory-review-status — Exploratory work awaiting expert review

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper states that companion symbolic scripts were used to check the arguments and calculations for internal consistency, while explicitly characterizing the work as exploratory and indicating that it would benefit from review by subject-matter experts. The scripts and their execution records are not supplied in the prepared unit.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\section*{Acknowledgments}`

## Revision and supersession relationships

### source-paper-1pn-orbital-dynamics--revision-status — No explicit supersession relationship supplied

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper refers to “corrected 1PN accounting” when motivating the nonzero pressure–volume term, but it does not identify a prior source or statement that it formally corrects, revises, or supersedes. No scoped supersession relationship can therefore be established from the prepared unit.

Sources:

- `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` — `\label{sec:hydro_kappa_PV}`

## Related topics and scripts

No related memory pages or script domains were supplied by the task.
