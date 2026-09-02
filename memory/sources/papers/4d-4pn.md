---
schema_version: 2
id: source-paper-4d-4pn
title: Paper 4D 4Pn
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/4d_4pn/paper/4d_4pn.tex
content_owner: ai_generated
last_updated: '2026-09-01'
generated_from_commit: 9f40c18eccc04f88d8c989d0f4ea091e3e619058
source_kind: paper
source_unit:
  id: paper-4d-4pn
  shape: file
  entrypoint: research/4d_4pn/paper/4d_4pn.tex
  unit_digest_sha256: 0572ab20d4ee72af2ced119ab76e56f417e0fe0664ff374fefa5b4144bcd4b05
  members:
  - path: research/4d_4pn/paper/4d_4pn.tex
    role: primary
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 52c0daa38baf1e69eba9b17ac1f2770944397074
    blob_size: 297747
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-paper-4d-4pn--conditional-conservative-scope — Conditional conservative 4PN assembly

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper addresses the nonspinning conservative two-body sector through order \(c^{-8}\) in a unified 4D toy-model program. It presents the result as a local instantaneous sector plus a hereditary tail sector within a declared reduction and closure hierarchy, not as an assumption-free consequence of a solved moving-throat bulk PDE.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:intro}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:intro-nonclaims}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:intro-taxonomy}`

## Source-unit map

- Entry point: `research/4d_4pn/paper/4d_4pn.tex`
- Role: `primary`
- Read mode: `semantic`
- Shape: monolithic paper file

## Key statements

### source-paper-4d-4pn--one-body-gate — Exact Schwarzschild one-body gate

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

Expansion of the exact isotropic Schwarzschild worldline Lagrangian through \(c^{-8}\) fixes the test-mass 4PN coefficient. Matching that gate requires the repair ledger
\[
\mu_{\rho,4}=\frac18,\qquad d_4=\frac{205}{16},\qquad
s_{34}=-\frac{15}{32},\qquad s_{26}=-\frac1{16},
\]
so the carried 3PN one-body package cannot be continued using only one new denominator coefficient.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-schwarzschild-exact-worldline}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-onebody-lagrangian-gate}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-onebody-repair-ledger}`

### source-paper-4d-4pn--quartic-legendre-compiler — Quartic chart compiler

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

For a perturbative Lagrangian \(L=L_0+\varepsilon L_1+\cdots+\varepsilon^4L_4\), with quadratic \(L_0\) and a constant Newtonian mass matrix, the paper derives the Legendre-transform coefficients through \(H_4\). This supplies the ordinary/Hamiltonian translation used in the local 4PN construction when the lower-order ledger is fixed.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:onebody-quartic-compiler}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-legendre-h4}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{app:quartic-compiler}`

### source-paper-4d-4pn--hamiltonian-first-firewall — Hamiltonian-first degree-ceiling firewall

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

The imported local target has symmetric-mass-ratio degree ceilings
\[
(4,3,3,2,2)
\]
across the \(G/r\) through \(G^5/r^5\) blocks in the Hamiltonian chart, while its ordinary-chart translation has ceilings
\[
(4,4,4,4,2).
\]
The extra ordinary-chart \(\nu^4\) terms in the \(G^2/r^2\), \(G^3/r^3\), and \(G^4/r^4\) blocks are translation effects. They must not be interpreted as terms that the original constant-coefficient comparable-mass scaffold should span directly; the prescribed order is Hamiltonian lift first and ordinary translation afterward.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-hamiltonian-degree-ceilings}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-ordinary-degree-ceilings}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-hamiltonian-first-scheme}`

### source-paper-4d-4pn--hamiltonian-span — Hamiltonian-first local span

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper reports that the exchange-symmetric generic-frame Hamiltonian scaffold reaches all 47 comparable-mass interaction coefficients and has a 92-dimensional null family. It characterizes that family as COM-blind algebraic freedom and reports a canonical sparse representative. The sealed paper states the rank counts and describes an ideal-membership audit, but the matrices and audit transcripts needed to independently establish those computations are not present in this source unit.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-interaction-rank-target}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-hamiltonian-total-rank-nullity}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:local-canonical-slice}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{app:local-referee-reconstruction-status}`

### source-paper-4d-4pn--ordinary-seed-alignment — Ordinary-chart seed alignment

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

With the lower-order ledger fixed, the paper derives the canonical residual translation
\[
\Delta L_{4,\mathrm{loc}}^{\mathrm{can}}
=-\Delta H_{4,\mathrm{loc}}^{\mathrm{can}}.
\]
It assigns the remaining ordinary-chart mismatch to an aligned-seed correction, specifies enlarged structured seed spaces, and states their maximal ranks. The in-paper theorem supplies the structured-lift claim, while the matrices and audit transcripts that would independently establish the reported ranks are not included in the prepared unit; the combined closure claim therefore remains provisional.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-sign-flip}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{thm:app-ordinary-translation-sign-flip}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-seed-alignment-correction}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-structured-seed-spaces}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-structured-seed-ranks}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{thm:app-structured-aligned-seed-lift}`

### source-paper-4d-4pn--local-reconstruction — Reported local reconstruction

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Within its declared closure hierarchy and fixed chart, the paper analytically assembles a natural generic-frame seed, a structured aligned-seed correction, and a translated canonical residual. It states as an analytic theorem that exact COM reduction of this candidate reproduces the complete fixed-chart local 4PN target block by block. The separately described referee suite is presented as verification of, not a replacement for, that analytic argument; because its executable artifacts and literal outputs are not members of the prepared unit, the reconstruction claim remains provisional here.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:local-final}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-generic-ordinary-candidate}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-referee-identity}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-local-final-theorem}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{thm:app-local-referee-reconstruction}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:verification-does-not-verify}`

### source-paper-4d-4pn--hereditary-bridge — Hereditary coefficient bridge

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper identifies the compatible conservative hereditary coefficient on the canonically normalized STF quadrupole branch as
\[
C_{\mathrm{tail}}=\frac{GM}{2c^3}\gamma_{\mathrm{quad}}^{\mathrm{eff}},
\]
so its closure scheme introduces no second quadrupole-normalization datum at 4PN. Separately, it derives the GR arithmetic identity
\[
\frac{GM}{2c^3}\frac{2G}{5c^5}
=\frac{G^2M}{5c^8}.
\]
The general branch identification is a closure-level assertion, while the GR substitution is an explicit arithmetic derivation; the combined statement is conservatively classified as provisional. Whether the model realizes the required quadrupole value is recorded separately as open.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-tail-exact-bridge}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-tail-burke-thorne-gr}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-tail-gr-bridge}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:tail-no-new-datum}`

## Computed evidence represented by the source

### source-paper-4d-4pn--referee-archive-unavailable — Referee computations are described but not supplied

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

not supplied in the prepared unit

The paper names a staged SymPy archive and a grouped Mathematica mirror intended to check the one-body gate, compiler, Hamiltonian lift, aligned-seed lift, local reconstruction, and hereditary bridge. None of those scripts, invocations, or literal outputs is a member of this prepared source unit, so their reported checks cannot be classified here as measured evidence.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{tab:verification-4pn-artifacts}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:verification-verifies}`

## Assumptions, exclusions, and open questions

### source-paper-4d-4pn--quadrupole-normalization-open — Moving-throat normalization and transport remain open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The passive/outgoing moving-throat branch has not been shown to yield
\[
\gamma_{\mathrm{quad}}^{\mathrm{eff}}=\frac{2G}{5c^5},
\]
or the equivalent stated \(P_0\) normalization. The grouped conservative operator moments and outgoing-transfer moments for the \(P_{20}\oplus P_{21}\oplus P_{22}\) branch also remain to be obtained from a completed moving-throat PDE. In the paper’s toy hereditary parameterization, a scalar transport factor \(\Theta_{\mathrm{tail}}\) remains as a possible model-side mismatch; the canonical branch requires \(\Theta_{\mathrm{tail}}=1\).

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:fixed-open-still-open}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:fixed-open-gamma-target}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:fixed-open-P0-target}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:fixed-open-grouped-operator}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:fixed-open-grouped-transfer}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-tail-toy-transport}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{eq:4pn-tail-theta-one}`

### source-paper-4d-4pn--scope-exclusions — Excluded regimes and interpretations

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The result excludes spin, generic many-body response, strong-field nonperturbative completion, and higher-PN sectors. It does not establish a first-principles moving-throat origin for all grouped conservative data, a final uniqueness or interpretation theorem for aligned-seed ordinary representatives, or a new dissipative or radiative theorem beyond the narrowed quadrupole route.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:intro-nonclaims}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:verification-does-not-verify}`

## Revision and supersession relationships

### source-paper-4d-4pn--lower-order-continuity — Continuity with the lower-order program

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper carries forward its declared reduction hierarchy and lower-order Newtonian through 3PN ledger. It treats the surviving 2.5PN STF quadrupole-normalization gate as the same unresolved interface controlling the 4PN tail. It presents this as continuity and extension of the program rather than an explicit supersession of the earlier papers as complete source units.

Sources:

- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:intro-motivation}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:interface-25pn-pde}`
- `research/4d_4pn/paper/4d_4pn.tex` — `\label{sec:fixed-open-frozen}`

## Related topics and scripts

No related memory pages or script domains were supplied by the task.