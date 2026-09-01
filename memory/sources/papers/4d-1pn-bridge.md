---
schema_version: 2
id: source-paper-4d-1pn-bridge
title: Paper 4D 1Pn Bridge
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
content_owner: ai_generated
last_updated: '2026-08-25'
generated_from_commit: e15f5a358d03e0f5ca0061c9316e690758e7e625
source_kind: paper
source_unit:
  id: paper-4d-1pn-bridge
  shape: file
  entrypoint: research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
  unit_digest_sha256: bb493bcfba08a21a861b35b8574ffd2184bcf3f46d75c2f05b85271681f4d15e
  members:
  - path: research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
    role: primary
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 3b57094fb54c0d4a42e9913cb188d056a1494363
    blob_size: 183816
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-paper-4d-1pn-bridge--bridge-purpose — Bridge from the unified 4D model to effective coefficients

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper presents itself as a derivation supplement connecting a unified four-spatial-dimensional brane–bulk throat model to coefficients used in earlier weak-field and first-post-Newtonian descriptions. It separates coefficients claimed to follow from specified topology, optical matching, thermodynamics, or wave kinematics from internal-response quantities that remain protocol-dependent.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:intro}`

## Source-unit map

- Entrypoint: `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: monolithic paper file

## Key statements

### source-paper-4d-1pn-bridge--added-mass-topology — Added mass distinguishes the throat and bubble geometries

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For quasi-static, low-Mach, incompressible, irrotational exterior flow with a no-penetration boundary, the paper derives \(\kappa_{\mathrm{add}}(d)=1/(d-1)\). A throat uniform in the extra coordinate reduces slice-by-slice to the three-dimensional sphere problem and has \(\kappa_{\mathrm{add}}=1/2\); a counterfactual compact four-dimensional bubble has \(\kappa_{\mathrm{add}}=1/3\). Compressibility, radiation, leakage, and nonuniform extra-dimensional structure lie outside this coefficient’s stated regime.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:added_mass}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{app:added_mass:d_ball}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{app:added_mass:throat}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{app:added_mass:bubble}`

### source-paper-4d-1pn-bridge--optical-eos-exponent — Weak-field optical matching selects \(n=5\)

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Under the barotropic closure \(P=K_{\rm EOS}\rho^n\), the acoustic-optics identification \(N=c_0/c_s\), linear weak-field hydrostatic balance, and \(c_0=c\) in the final comparison, the paper derives \(N(r)\simeq1+[(n-1)/2]GM/(c^2r)\). Matching the stipulated GR optical coefficient \(2GM/(c^2r)\) fixes \(n=5\), while leaving the dimensional normalization \(K_{\rm EOS}\) unfixed. The paper limits this conclusion to the optical sector and warns against promoting the resulting optical metric to the complete massive-body metric.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:N_r_general_n}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:N_GR_target}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:n_equals_5}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:optics_eos}`

### source-paper-4d-1pn-bridge--eih-family — Cross-tensor matching yields a one-parameter wake family

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Matching the paper’s wake-overlap cross coefficients to the EIH targets \(C_\parallel=-7/2\) and \(C_L=-1/2\) gives \(\alpha^2=\tfrac34(1-a_H^2)\) and \(K_{\rm vec}=2/[\pi^2(1-a_H^2)]\), with \(0\leq a_H^2<1\) when \(K_{\rm vec}>0\) and \(\alpha^2\geq0\). The cross-tensor match alone therefore does not select a unique point.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:EIH_target_cross}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:family_alphaK}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{app:eih-matching-algebra}`

### source-paper-4d-1pn-bridge--thermodynamic-collapse — Thermodynamic closure collapses the wake family

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Using the paper’s definition of \(\alpha^2\) as an available-wake fraction together with the \(n=5\) thermodynamic relation yields \(\alpha^2=3/4\). Substitution into the EIH family forces \(a_H^2=0\) and \(K_{\rm vec}=2/\pi^2\). These values depend on both the stated wake-basis convention and the paper’s thermodynamic identification.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:1pn-interaction-thermo}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:alpha_value}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:alphaK_final}`

### source-paper-4d-1pn-bridge--wave-supported-energy — Wave-supported rest energy gives the \(3/8\) coefficient

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a trapped excitation with linear dispersion and the stated comoving group-velocity condition, the paper derives \(E(v)=\gamma E_0\), hence \(E=E_0+\tfrac12(E_0/c^2)v^2+\tfrac38(E_0/c^2)v^4/c^2+\cdots\). This establishes the \(3/8\) coefficient within that wave-supported kinematic model, but does not establish Lorentz invariance of the full interacting theory.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:E_gamma_E0}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:E_series}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{app:standing_wave_v4}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:sr_wave_mass}`

### source-paper-4d-1pn-bridge--projected-open-system — Projection remains open and Poisson behavior is regime-limited

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Projecting the bulk conservation laws with a fixed weight gives an exact brane continuity equation with an extra-dimensional leakage source. The projected momentum balance also contains leakage terms and a Reynolds-type stress because projection does not commute with nonlinear products, so a closed three-dimensional perfect-fluid description requires additional closure assumptions. A Poisson-like equation is not postulated: for the longitudinal velocity potential, linearized low-Mach and quasi-static assumptions give \(\rho_0\nabla_3^2\phi_3\approx S_\rho\). The paper treats this only as a diagnostic limit, spoiled when leakage or unresolved stress is important, transverse or vortical dynamics are active, or quantum pressure is non-negligible.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:brane_continuity_open}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:brane_leakage_source}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:Rij_def}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:brane_euler_open}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:poisson_candidate}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:poisson_candidate}`

### source-paper-4d-1pn-bridge--adiabatic-pv-closure — A declared adiabatic closure fixes \(\kappa_{\rm PV}\)

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For the paper’s one-degree-of-freedom, fixed-aspect-ratio, adiabatic throat closure, and using the Newtonian ledger value \(\kappa_\rho=1\) to convert \(d\ln F/d\ln\rho\) into \(\kappa_{\rm PV}\), equilibrium and density-response matching give \(x=2/11\), \(E_w:E_f:E_{\rm PV}=11:2:5\), \(\kappa_{\rm PV}=3/2\), and \(d\ln a/d\ln\rho=-57/64\). The result assumes any density-independent stabilizer energy is subdominant; a comparable contribution can shift or eliminate \(\kappa_{\rm PV}=3/2\). The absolute throat scale and nonadiabatic response remain uncalibrated.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:kappaPV:adiabatic1dof}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:kappaPV_closure_target_slope}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:kappaPV_closure_x_solution}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:kappaPV_closure_partition}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{eq:kappaPV_closure_dlnA_n5}`

## Computed evidence represented by the source

### source-paper-4d-1pn-bridge--verification-evidence-boundary — Verification chain is incomplete in the prepared unit

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper names Mathematica and Python verification programs and records commands intended to run them, but the prepared unit contains neither those executable members nor readable literal transcripts with operands and residuals. The complete command/output/interpretation chain is not supplied in the prepared unit, and no result in this capsule is classified as measured.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — heading `\subsection*{Companion verification bundle}`

## Assumptions, exclusions, and open questions

### source-paper-4d-1pn-bridge--response-closure-open — General internal response remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

Outside the explicit adiabatic example, geometry, breathing, leakage, reservoir, localization, and mouth-response coefficients remain dependent on a declared drive, constraint, stabilization model, and frequency regime. The paper leaves their calibration and the extraction of a stable low-frequency response operator to future calculations.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:undetermined:philosophy}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:undetermined:Zeff}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:undetermined:program}`

### source-paper-4d-1pn-bridge--preferred-frame-open — Full preferred-frame suppression is not established

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The wave-clock and wave-ruler construction supplies a controlled kinematic mechanism for time dilation, length contraction, and a Michelson–Morley-type cancellation when dispersion is effectively linear, confinement remains phase-locked, and leakage, drag, anisotropy, and open-system effects are negligible. The paper explicitly leaves interaction-level Lorentz covariance and residual preferred-frame signatures unresolved.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:preferred_frame_violations}`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{app:standing_wave_preferred_frame}`

## Revision and supersession relationships

### source-paper-4d-1pn-bridge--downstream-ontology — Downstream notation follows corrected matter–gauge ontology

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper declares itself downstream of corrected Paper 7/Paper 8 conventions: electric-charge sign is assigned to puncture orientation, circulation belongs to the magnetic or vortical sector, the zero-mode Maxwell reduction is only a controlled far-field brane limit, and the historical scalar mass-dressing symbol \(q\) is distinguished from electric charge and renamed \(\kappa_\rho\). The prepared unit does not contain the referenced earlier papers, so this capsule records the relationship only as the present paper states it and does not infer whole-page supersession.

Sources:

- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` — `\label{sec:intro}`

## Related topics and scripts

- Related memory topics: none supplied by the task.
- Related script domains: none supplied by the task.
- The paper names companion verification programs, but they are not members of this sealed source unit.