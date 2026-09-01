---
schema_version: 2
id: source-pde-ledger-v3-s10
title: Pde Ledger V3 S10
type: source_capsule
lifecycle: current
memory_review: ai_draft
sources:
- research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl
- research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out
- research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex
- research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.premises
- research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py
- research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
- research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
- research/pde_ledger_v3/scripts/S10_exports.py
- research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out
- research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
- research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
- research/pde_ledger_v3/steps/S10_PREREGISTERED_PREDICTION.md
- research/pde_ledger_v3/steps/S10_two_transverse_photons.md
content_owner: ai_generated
last_updated: '2026-08-25'
generated_from_commit: e15f5a358d03e0f5ca0061c9316e690758e7e625
source_kind: step_family
source_unit:
  id: pde-ledger-v3-s10
  shape: bundle
  entrypoint: research/pde_ledger_v3/steps/S10_PREREGISTERED_PREDICTION.md
  unit_digest_sha256: 3ef1a7b11e3847e3b1b4724855e6910186feab2b30e236ad30a8f5427230c00d
  members:
  - path: research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl
    role: engine
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 7ce0af3743467194d7b4efcc1a644953738f5d2c
    blob_size: 54590
  - path: research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out
    role: measurement
    read_mode: excerpt
    mode: '100644'
    object_type: blob
    blob_oid: 8bf3cc587cf908703c51e1d3b20b5502f5a46174
    blob_size: 660892
  - path: research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex
    role: paper_card
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 51735f91c16d045dae599491fd1ee7c553c8f6ae
    blob_size: 25027
  - path: research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.premises
    role: premises
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 3b648d92f9cacd38fcb8c4e2bec87e92923cf241
    blob_size: 1450
  - path: research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py
    role: engine
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: b8c2d1979f5c7574f9758384bf0503f10b91c01f
    blob_size: 95219
  - path: research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py
    role: comparator
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 1b4e965ba48e328f3222be74f21caea1c1f30bac
    blob_size: 37011
  - path: research/pde_ledger_v3/scripts/S10_cross_engine_comparator_repoint_ablation.py
    role: control
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 21478246bae694c6be0f1e2d3bb68131793a4667
    blob_size: 2292
  - path: research/pde_ledger_v3/scripts/S10_exports.py
    role: dependency
    read_mode: identity_only
    mode: '100644'
    object_type: blob
    blob_oid: fefd528ac1294a239251fa7bc8e8cabbc784cefd
    blob_size: 1284971
  - path: research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out
    role: measurement
    read_mode: excerpt
    mode: '100644'
    object_type: blob
    blob_oid: 132a8588506eecfce1a248921833efa13a678a12
    blob_size: 725576
  - path: research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out
    role: measurement
    read_mode: excerpt
    mode: '100644'
    object_type: blob
    blob_oid: de8939c2a7b3411f151a28c2e7a12bbf964621d5
    blob_size: 717815
  - path: research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator_repoint_ablation.out
    role: measurement
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: a400fde5fef7426418275f77f4a77a98ba32c486
    blob_size: 1607
  - path: research/pde_ledger_v3/steps/S10_PREREGISTERED_PREDICTION.md
    role: step_record
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: a52610aa8c01fb7ef7ad648f8d4c087f54652368
    blob_size: 2108
  - path: research/pde_ledger_v3/steps/S10_two_transverse_photons.md
    role: step_record
    read_mode: semantic
    mode: '100644'
    object_type: blob
    blob_oid: 7fcc87fa6b391ef2bf1921073ff8e821dec526fe
    blob_size: 49255
extractor_version: 1
---

> Generated capsule. Refresh from the source unit; do not hand-edit.

## Purpose and scope

### source-pde-ledger-v3-s10--conditional-main-spectrum — Conditional MAIN spectrum

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

For the supplied curl-only quadratic in-plane action, positive coefficients, nonzero real wavevector, and generic strata, the sources report roots \(0\) and \(\mu_R|k|^2/\rho_{\rm br}\), with \(D-1\) exactly transverse directions at the nonzero root in MAIN cases \(D=2,3,4,5\). S10 does not select physical \(D=3\). No exact recorded invocation is supplied, so the transcript-reported MAIN result remains provisional here.

Sources:

- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — `\label{step:S10}`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Record boundary`

### source-pde-ledger-v3-s10--longitudinal-zero-mode — Longitudinal zero mode remains

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The sources report one longitudinal direction at the zero root: curl-only stiffness removes its restoring stiffness rather than the degree of freedom. Its disposition is deferred to S11, and S10 does not establish a complete Maxwell light sector. This transcript-backed result remains provisional because no recorded invocation is supplied.

Sources:

- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — `\label{step:S10}`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Measured MAIN spectrum and count`

## Source-unit map

- Step records: preregistration and current interpretation (`semantic`).
- Paper card: derivation and publication-facing scope (`semantic`).
- Engines: SymPy and Wolfram constructions (`semantic`).
- Premises: frozen SymPy premise surface (`semantic`).
- Comparator and repoint control: joined-name comparison and sensitivity test (`semantic`).
- Measurements: bounded engine and comparator excerpts plus the semantic repoint-ablation transcript.
- Export dependency: identity only; its contents were not read.

## Key statements

### source-pde-ledger-v3-s10--preregistered-predictions — Every preregistered prediction held

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The interpretive record reports that every preregistered prediction held in both transcripts: MAIN’s \(D-1\) pattern, FULLGRAD’s single root with total nullity \(D\), coefficient-rescaled roots with unchanged nullities, the stated curl normalization, and nondegenerate \(D=2\) behavior. Later controls narrow the headline interpretation and the generality of exact transversality; they do not retract the factual FULLGRAD longitudinal-polarization prediction.

Sources:

- `research/pde_ledger_v3/steps/S10_PREREGISTERED_PREDICTION.md` — heading `## Main computation`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Prospective pre-registration`

### source-pde-ledger-v3-s10--control-boundary — Controls bound exact transversality

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The sources report that FULLGRAD shares MAIN’s nonzero root and \(D-1\) transverse nullity at \(D=3,4\) while also propagating the longitudinal direction; DIVONLY reverses the sectors. ANISO preserves total propagating nullity \(D-1\) but generically reduces exact transversality to \(D-2\), restored on its allowed axial stratum. These controls bound the role of stiffness and inertia, but the result remains provisional without a recorded invocation.

Sources:

- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — `\label{step:S10}`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Controls measured from their actions`

## Computed evidence represented by the source

### source-pde-ledger-v3-s10--evidence-chain — Prepared computation chain

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The readable engines construct Euler–Lagrange systems, equation-of-motion and quadratic-form matrices, spectra, rank/nullity objects, dimensional checks, and tagged output. Excerpts expose literal engine tags and comparator totals, and the step record supplies interpretation. Exact engine and comparator invocations are not supplied, so the conditional MAIN spectrum, longitudinal zero mode, control boundary, and MAIN dimensional result are not promoted to `measured` or `derived`.

Sources:

- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — function `main`
- `research/pde_ledger_v3/mathematica/S10_brane_mode_spectrum_mathematica_audit.wl` — identifier `WL_S10_RUN_PAIRS`
- `research/pde_ledger_v3/scripts/S10_cross_engine_comparator.py` — function `main`
- `research/pde_ledger_v3/scripts/out/S10_brane_mode_spectrum_sympy_audit.out` — tag `PY_S10_MAIN_D2_Q3_ROOT_COUNT`
- `research/pde_ledger_v3/mathematica/out/S10_brane_mode_spectrum_mathematica_audit.out` — tag `WL_S10_MAIN_D2_Q3_ROOT_COUNT`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## What was supplied and what was computed`

### source-pde-ledger-v3-s10--comparator-coverage — Comparator coverage is partial and negative overall

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The comparator transcript reports 388 agreements, 164 disagreements, and `FINAL_GUARD: FAIL`. There is no shared action-payload name; all 13 Euler–Lagrange rows fail representational comparison; all 13 route-account rows fail token-content comparison; and the underlying route constructions are not compared. The two matrix routes share one action and remain a coding-consistency control, not independent physical derivations.

Sources:

- `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out` — tag `CATEGORY: SUMMARY`
- `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out` — tag `FINAL_GUARD`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Live comparator: measured scope, not a blanket verdict`

### source-pde-ledger-v3-s10--main-dimensions — Conditional dimensional homogeneity

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

Under supplied \([u]=(1,0,0)\), the sources report MAIN values \([\rho_{\rm br}]=(-D,0,1)\), \([\mu_R]=(2-D,-2,1)\), and difference \((2,-2,0)\). This is action homogeneity under a shared, currently unfalsifiable field-dimension premise, not independent dimensional determination or whole-layer agreement. It remains provisional because no recorded invocation is supplied.

Sources:

- `research/pde_ledger_v3/paper/steps/S10_two_transverse_photons.tex` — `\label{step:S10}`
- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.premises` — distinctive text `The in-plane displacement has supplied dimension`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Dimensional result`

### source-pde-ledger-v3-s10--export-guard-boundary — S9-to-S10 guard has blind spots

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The S9-to-S10 guard compares three dependent SymPy dimension records. The interpretive record reports that it caught a form mutation but missed a dimension-preserving coefficient mutation. It is blind to common-mode premise changes, can lose assertions under optimized Python, and is neither exhaustive nor a second cross-engine check.

Sources:

- `research/pde_ledger_v3/scripts/S10_brane_mode_spectrum_sympy_audit.py` — function `main`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Export chain and what its guard can detect`

## Assumptions, exclusions, and open questions

### source-pde-ledger-v3-s10--physical-assumptions — Physical and substrate delivery gates remain open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The curl-only functional is supplied rather than derived from the substrate. The substrate still owes the joint quadratic operator on in-plane \(u\) and out-of-plane \(h\), including whether \(h\) couples to or is degenerate with the transverse pair; failure can change the headline mode count. Also open are equilibrium of the unstrained state, a stable positive stiffness sign, whether curl-only continuum mechanics requires internal angular momentum or couple stress, and whether the rotational-stiffness frame is admissible rather than an external preferred-orientation frame.

Sources:

- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Claims this step still does not establish`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Open requirements on interpreting the count as light`

## Revision and supersession relationships

### source-pde-ledger-v3-s10--name-coverage-history — Declared-name coverage superseded

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

Former declared-name coverage is superseded for comparator coverage by mechanical matching without an exception map. The current D12 worklist retains unrefuted \(D\leftrightarrow\)`braneDimension` and \(s\leftrightarrow\)`coefficientScale` pairs; injectivity remains open.

Sources:

- `research/pde_ledger_v3/scripts/out/S10_cross_engine_comparator.out` — tag `CATEGORY: D12_NONMECHANICAL_SYMBOL_WORKLIST`
- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Claims this step still does not establish`

### source-pde-ledger-v3-s10--q7-history — Historical Q7 counts are unreproducible

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

Historical before/after Q7 failure counts are unreproducible from the live tree because the old harness, pinned outputs, and configuration are absent. No historical count is carried forward.

Sources:

- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Claims this step still does not establish`

### source-pde-ledger-v3-s10--registry-history — Reduction registry replaced by distributed records

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The former reduction registry is absent. Its live replacement is distributed across the S9 structural export, S10 import, indexed S10 records, and cross-step dimension operands, while physical selection of \(D=3\) remains owed. The historical provenance defect concerned the absent registry and is not evidence that the physical dimension has been selected.

Sources:

- `research/pde_ledger_v3/steps/S10_two_transverse_photons.md` — heading `## Registry disposition`

## Related topics and scripts