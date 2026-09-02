---
schema_version: 2
id: source-paper-4d-em-fields
title: Paper 4D Em Fields
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/4d_em_fields/paper/4d_em_fields.tex
content_owner: ai_generated
last_updated: '2026-09-01'
generated_from_commit: 60e31504b8425840cf6b5e814a0c518d68ab2db6
source_kind: paper
source_unit:
  id: paper-4d-em-fields
  shape: file
  entrypoint: research/4d_em_fields/paper/4d_em_fields.tex
  unit_digest_sha256: 74a9584271b5fda78b24083619ce202e8d8a86154103edbda3e6cd42174caddf
  members:
  - path: research/4d_em_fields/paper/4d_em_fields.tex
    role: primary
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 582224346384289b734db22779a08ad3b8fc8c62
    blob_size: 86841
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-paper-4d-em-fields--localized-maxwell-scope — Localized Maxwell sector

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper isolates a \(4+1\)-dimensional gauge sector whose Maxwell kinetic term is weighted by a transverse localization profile \(Z(w)\). It seeks a controlled long-distance, low-frequency \(3+1\) brane limit and quantifies departures caused by transverse modes. The Gaussian profile is the worked analytic case; microscopic spin and mixed-twist defect-core dynamics are outside scope.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:intro}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:conclusion}`

## Source-unit map

- Entry point: `research/4d_em_fields/paper/4d_em_fields.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: monolithic paper file

## Key statements

### source-paper-4d-em-fields--bulk-equations-and-identities — Bulk equations and consistency identities

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Varying the stated localized action gives
\(\partial_M(ZF^{MN})+\xi^{-1}\partial^N(\partial\!\cdot\!A)=\mu_0J^N\).
The potential definition supplies the off-shell Bianchi identities, while antisymmetry makes the double divergence of \(ZF^{MN}\) vanish. Consequently, the divergence of the field equation relates \(\Box(\partial\!\cdot\!A)\) to \(\partial_NJ^N\), enforcing current conservation in Lorenz gauge.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:action_em}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:5d_maxwell_eom_full}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:bianchi}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:div_identity_full}`

### source-paper-4d-em-fields--controlled-brane-reduction — Controlled brane Maxwell reduction

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Integrating the brane components over \(w\) gives an equation containing the transverse boundary term and, after division by \(Z_{\mathrm{int}}\), the paper’s exact weighted brane equation. Standard \(3+1\) Maxwell form follows only after imposing vanishing transverse boundary flux, axial gauge, Lorenz gauge, a brane-dominant \(w\)-independent zero mode, and no transverse current. Within that scope,
\(\mu_0^{\mathrm{eff}}=\mu_0/Z_{\mathrm{int}}\); for \(Z=e^{-w^2/\lambda^2}\), this becomes \(\mu_0/(\lambda\sqrt{\pi})\).

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:integrated_exact}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:brane_maxwell_exact_weighted}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:zero_mode_ansatz}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:brane_maxwell_gf}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:mu0_eff}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:eom_w}`

### source-paper-4d-em-fields--thickness-normalized-charge — Thickness-normalized brane charge

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a fixed defect branch \(q_\star=\eta_Qe_\star\), canonical normalization of the reduced zero mode gives
\(q_{\mathrm{eff}}=q_\star/\sqrt{Z_{\mathrm{int}}}\). The branch sign remains fixed while the magnitude of the canonically normalized reduced coupling depends on localization thickness.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:canonical_charge}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:qeff_canonical}`

### source-paper-4d-em-fields--gaussian-kk-spectrum — Gaussian KK spectrum and brane selection rule

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The Gaussian transverse Sturm–Liouville equation reduces to the Hermite equation, yielding \(f_n(w)=H_n(w/\lambda)\) and \(m_n^2=2n/\lambda^2\). Weighted normalization gives \(\|f_n\|^2=\lambda\sqrt{\pi}\,2^n n!\). A centered brane source has \(c_{2m+1}=0\), so odd modes decouple from that source without implying that microscopic odd-\(w\) structure is absent.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:sl_problem}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:hermite_spectrum}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:hermite_norm}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:odd_vanish}`

### source-paper-4d-em-fields--coulomb-yukawa-correction — Coulomb plus Yukawa correction

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a static point source on a fixed branch, the KK sum produces a Coulomb zero mode plus massive Yukawa terms. At \(r\gg\lambda\), the leading relative correction is
\(\tfrac12 e^{-2r/\lambda}\), followed at order \(e^{-2\sqrt{2}r/\lambda}\). This coefficient and range follow from the Gaussian profile rather than being independent fit parameters.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:A0_KK_sum}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:ratio_cn}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:leading_yukawa}`

### source-paper-4d-em-fields--retarded-brane-response — Covariant causal brane response

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For conserved time-dependent brane sources, the \(n=0\) amplitude obeys the massless Maxwell wave equation and the \(n\geq1\) amplitudes obey massive wave equations. The summed retarded propagator depends on the brane Lorentz scalar \(k^2\), together with the retarded prescription. Massive modes add response inside the forward light cone while leaving the wavefront causal.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:mode_wave_eq}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:Deff_retarded}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:KG_retarded}`

### source-paper-4d-em-fields--matter-current-consistency — Matter current closes Maxwell consistency

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The representative minimally coupled brane matter model has a \(U(1)\) Noether current with \(j^0=|\psi|^2\) and the stated gauge-covariant spatial components. Its off-shell identity implies \(\partial_\mu j^\mu=0\) on shell. Embedding \(J_\psi^\mu=q_\star j^\mu\) as \(J^\mu(x,w)=J_\psi^\mu(x)\delta(w)\) with \(J^w=0\) supplies the conserved bulk source required by the Maxwell divergence identity.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:j0}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:ji_covariant}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:charge_current_def}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:noether_identity}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:J_bulk_from_brane}`

## Computed evidence represented by the source

### source-paper-4d-em-fields--verification-evidence-boundary — Verification artifacts unavailable

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper describes six Wolfram Language checks and a standard `wolframscript -file <file.wl>` invocation pattern, but the scripts, recorded invocations, literal transcripts, and separate interpretation records required to assess their residuals or pass/fail output are not supplied in the prepared unit.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{app:wl_suite}`

## Assumptions, exclusions, and open questions

### source-paper-4d-em-fields--maxwell-regime-boundary — Maxwell regime boundary

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The reduced Maxwell law requires zero-mode dominance, \(J^w=0\), sufficiently decaying transverse boundary flux, and a conserved brane current. Significant \(w\)-dependence, transverse flux, or mixed \((A_w,J^w,F_{\mu w})\) structure lies outside that pure-Maxwell reduction. These conditions suppress the mixed channels as a controlled far-field approximation; they do not erase them from the microscopic defect ontology. The paper presents \(r\gg\lambda\) and \(\omega\ll2/\lambda\) as sufficient regime estimates rather than exact universal boundaries.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:intro}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:integrated_exact}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:eom_w}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:discussion}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:conclusion}`

### source-paper-4d-em-fields--profile-microphysics-open — Localization microphysics remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The paper treats \(Z(w)\) and \(\lambda\) phenomenologically. Deriving or constraining them from throat or defect microphysics, and developing spin or mixed-twist defect-core dynamics, is explicitly left for future work.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:discussion}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:conclusion}`

## Revision and supersession relationships

### source-paper-4d-em-fields--charge-ontology-revision — Revised charge ontology without formal supersession target

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper explicitly adopts an updated post–Paper VII ontology: charge is labeled by the fixed puncture orientation \(\eta_Q\) and microscopic coupling \(q_\star\), while localization controls \(q_{\mathrm{eff}}\). It rejects circulation, throat radius, and breathing variables as charge definitions. The prepared unit does not identify a concrete earlier statement as formally superseded, so this capsule records the revision without assigning supersession.

Sources:

- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{sec:intro}`
- `research/4d_em_fields/paper/4d_em_fields.tex` — `\label{eq:qeff_canonical}`

## Related topics and scripts

No related memory pages or script domains were supplied by the task.