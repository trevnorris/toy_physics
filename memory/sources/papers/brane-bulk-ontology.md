---
schema_version: 2
id: source-paper-brane-bulk-ontology
title: Paper Brane Bulk Ontology
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
content_owner: ai_generated
last_updated: '2026-09-01'
generated_from_commit: 60e31504b8425840cf6b5e814a0c518d68ab2db6
source_kind: paper
source_unit:
  id: paper-brane-bulk-ontology
  shape: file
  entrypoint: research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
  unit_digest_sha256: 8e775db4fd931759fccc2b22655db50a2602a14dbaa1dc8de2e5dd4f9bd91a05
  members:
  - path: research/brane_bulk_ontology/paper/brane_bulk_ontology.tex
    role: primary
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 39a2d6427f820112d341f30bf49a9dc4a7b0b6c1
    blob_size: 134349
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-paper-brane-bulk-ontology--unified-throat-ontology — Unifying geometric proposal

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper proposes that a defect is a throat connecting an observable three-dimensional brane at \(w=0\) to a four-dimensional superfluid bulk. Its nearly spherical mouth and approximately cylindrical interior are presented as different projections of this geometry: the gravitational sector probes effective brane fields near the mouth, while the electromagnetic sector probes bulk cavity modes. The paper also quotes the Paper IV result that enthalpy minimization at fixed charge selects \(L/a=\sqrt{2}\pi/x_{01}\approx1.85\), reinterpreting it as a geometric property of the throat rather than rederiving the variational result.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:intro_ontology}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:aspect_ratio}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:aspect_ratio_result}`

### source-paper-brane-bulk-ontology--toy-model-scope — Toy-model scope

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The construction is an ontology for a weak-field, slow-motion toy universe, not a microscopic theory or a literal description of our universe. The paper aims to organize earlier gravitational and electromagnetic results and constrain future calculations rather than replace general relativity or electromagnetism.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:intro_scope}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:relation_other_models}`

## Source-unit map

- Entry point and sole member: `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: monolithic TeX paper

## Key statements

### source-paper-brane-bulk-ontology--dimensional-reduction — Gaussian dimensional-reduction result

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For the paper’s separable Gaussian throat profile with \(K(w)=1\), the supplied integrations give \(M=\pi^2La^3\rho_0\) and \(Q_{20}/M=(3/10)\varepsilon a^2\). Combining the derived multipole scaling with the stated far-field expansion yields an anisotropic potential correction of order \(\varepsilon(a/r)^2P_2(\cos\theta)\) for \(r\gg a,L\). The coefficient depends on the profile and normalization conventions, and the Gaussian is only a caricature of the physical throat.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_mass_result}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_q_over_m}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:phi_farfield_scaling}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{app:gaussian_multipoles}`

### source-paper-brane-bulk-ontology--one-pn-spherical-sink — Spherical-sink retrieval at 1PN

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper interprets monopole dominance and \((a/r)^2\)-suppressed multipoles as explaining why the far-field throat appears as a spherical sink at 1PN even though its bulk interior is cylindrical. It further asserts that the leading Newtonian and 1PN results of Papers I–III are insensitive to detailed throat geometry at fixed mass; this connection is interpretive rather than a new derivation of those earlier 1PN results.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:connection_papers1to3}`

### source-paper-brane-bulk-ontology--throat-mode-spectrum — Cylindrical throat modes

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Within the paper’s straight-cylinder, axisymmetric, thin-throat approximation and Dirichlet conditions at \(r=a\) and \(w=0,L\), separation of the linear acoustic equation yields Bessel radial modes and axial standing waves. The spectrum is \(\omega_{mn}^2=c_s^2[(x_{0m}/a)^2+(n\pi/L)^2]\), with fundamental profile \(J_0(x_{01}r/a)\sin(\pi w/L)\).

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:throat_modes}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_wave_rw}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_full_mode}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_mode_spectrum}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_fundamental_mode}`

### source-paper-brane-bulk-ontology--mass-charge-projections — Mass and charge as distinct projections

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

In the proposed dictionary, mass is associated with drainage and localized volume or density deficit. Charge is proportional to circulation around a closed brane loop encircling the throat mouth; in the absence of vorticity creation or reconnection, the circulation is treated as an integer-valued conserved input label. By Stokes’ theorem it equals vorticity flux through a surface bounded by that loop, and it approximates throat-threading flux when the throat carries the only significant vorticity near the loop. These are interpretive identifications within the toy model, not microscopic derivations.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:charge_flux}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:mass_flux}`

### source-paper-brane-bulk-ontology--brane-magnetostatics — Postulated brane-wake magnetostatics

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper postulates a Coulomb-gauge vector Poisson equation for a localized moving throat and identifies \(\mathbf A\) with a brane-confined transverse wake. It then quotes the corresponding Green-function solution, \(\mathbf A\propto\mathcal Q\mathbf u/r\), whose curl has the Biot–Savart form \(\mathbf B\propto\mathcal Q(\mathbf u\times\mathbf r)/r^3\). This dictionary avoids requiring far-zone bulk vorticity, but induction and radiation remain deferred to a dynamical completion.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:brane_magnetostatics}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:vector_poisson_A}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:A_from_vT}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:biot_savart_point}`

### source-paper-brane-bulk-ontology--wake-mixing-weight — Positive wake-mixing interpretation

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper adopts the corrected value \(\alpha^2=3/4\) as a real longitudinal/transverse weighting in the 1PN vector kernel. Its two-mode illustration shows that a positive weight is compatible with a positive-definite quadratic form, but does not repeat the earlier EIH matching that fixes the value.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{sec:wake_mixing_constraint}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:alpha_squared_value}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:app_two_mode_energy}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{app:two_mode_details}`

## Computed evidence represented by the source

### source-paper-brane-bulk-ontology--computed-evidence-availability — Rounded-funnel numerical result

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper reports numerical quadrature for a hard-bounded rounded funnel with \(a=1\), \(L/a=2\), \(\rho_0=1\), and \(\varepsilon=0.1\), obtaining \(M\approx9.98\), \(Q\approx0.143\), and \(\alpha_2=Q/(Ma^2)\approx0.014\), quoted in the main text as \(\alpha_2\simeq1.4\times10^{-2}\). The prepared unit contains the asserted setup and results, but the command/output/interpretation chain is not supplied in the prepared unit: there is no exact recorded invocation, literal stdout, or separate interpretation record. The result therefore remains provisional rather than measured.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:rounded_funnel}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{app:physical_throat}`

## Assumptions, exclusions, and open questions

### source-paper-brane-bulk-ontology--phenomenological-closures — Adopted EOS and calibrated coupling

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The bulk framework adopts the polytropic equation of state \(P=K\rho^n\) with \(n=5\), attributed to the earlier Paper II optical matching rather than rederived here. Its effective three-dimensional Poisson equation uses \(G_{\mathrm{eff}}\), which the paper treats as calibrated by matching Newtonian gravity and 1PN corrections in Papers I–III rather than deriving the coupling from first principles.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:bulk_fields}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:poly_eos}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:brane_hypersurface}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{eq:3d_poisson}`

### source-paper-brane-bulk-ontology--two-pn-status — Full 2PN dynamics remain open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper does not derive a complete 2PN effective Lagrangian. It proposes a future double expansion in \(\epsilon\sim GM/(rc^2)\) and \(\delta\sim a/r\). If \(a\) is of order the gravitational radius, then \((a/r)^2\sim\epsilon^2\) in usual PN counting; the paper adds that somewhat larger \(a\) can still place these finite-size effects beyond leading 1PN terms at sufficiently large separations. The microscopic origin and stabilization of the brane, the throat-bottom boundary condition, multi-throat dynamics, strong-field evolution, mergers, and a complete electromagnetic sector also remain unspecified or outside scope.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:intro_scope}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:2pn_meaning}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:2pn_interpretation}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:bottom_throat}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:brane_stability}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{subsec:multi_defects}`

## Revision and supersession relationships

### source-paper-brane-bulk-ontology--wake-basis-correction — Earlier wake-basis interpretation corrected

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper states that earlier exploratory wake decompositions using a truncated basis could suggest a negative longitudinal weight, whereas the completed isotropic projector basis gives the corrected positive value \(\alpha^2=3/4\). The source describes this correction of the wake-basis interpretation but does not declare a formal page- or paper-level supersession relationship.

Sources:

- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{sec:wake_mixing_constraint}`
- `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` — `\label{app:two_mode_details}`

## Related topics and scripts

No related topic pages or script domains were supplied by the task.