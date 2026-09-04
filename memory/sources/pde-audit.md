---
title: PDE audit
type: source
status: current
sources:
  - research/pde_audit/CLAIM_CHECK_MATRIX.md
  - research/pde_audit/run_all.sh
  - research/pde_audit/scripts/run_all_audits.sh
  - research/pde_audit/mathematica/run_all_audits.sh
  - research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md
  - research/pde_audit/notes/stage_v2_26_program_status_after_audit.md
last_updated: 2026-09-03
---


## Purpose and scope

### source-pde-audit--verification-scope — Verification project scope

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

This bundle is a referee-facing verification stack for the PDE V2 program. It combines symbolic Python audits, selected Mathematica mirrors, branch-extraction adapters and validators, negative controls, reduced solver prototypes, simulation protocols, and status notes. The Python scripts are the primary executable audits; the Mathematica mirrors add coverage for selected algebra and fixture contracts but are not independent derivations. Its physical object is a finite brane-bulk throat or puncture with an open conduit and interior structure, not a surface dimple, depression, or hard-capped pocket. The bundle checks conditional algebra and interfaces rather than establishing a completed nonlinear moving-throat theory.

Sources:

- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `# PDE V2 Claim-Check Matrix`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- `research/pde_audit/notes/stage_v2_28_physical_picture_and_ontology.md` — heading `## 1. Core ontology`

## Source-unit map

### source-pde-audit--source-unit-workflow — Runner and interpretation workflow

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The top-level referee harness coordinates fixture verification, Python audits, the simulation suite, selected Mathematica mirrors, and combined summary reporting. The Python runner captures stage output and writes summaries and JSON artifacts; the Mathematica runner loops over the selected audit scripts and invokes them with `math -script`. The notes and claim-check matrix own paper-facing interpretation. This is implemented workflow, not measured scientific evidence.

Sources:

- `research/pde_audit/run_all.sh` — nearby identifier `run_top_check`
- `research/pde_audit/scripts/run_all_audits.sh` — nearby identifier `OUTPUT_DIR`
- `research/pde_audit/mathematica/run_all_audits.sh` — nearby identifier `math -script`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Reproducibility Artifacts`

## Key statements

### source-pde-audit--parent-wall-status — Parent wall dynamics remain split

Status: `lifecycle=current` · `evidence=derived` · `memory_review=ai_draft`

With geometry present only through the confinement potential, variation supplies a wall force but no autonomous wall PDE. Adding the quadratic distributed wall action supplies the stated linear wall equation, boundary terms, modal split, and positivity gates. The effective linear closure is supported within its assumptions, while strict parent-level moving-throat dynamics still require promotion of a wall or throat action.

Sources:

- `research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md` — heading `## Executive result`
- `research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md` — heading `## 2. Variation of the confinement-only parent term`
- `research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md` — heading `## 3. Variation after adding the quadratic wall action`

### source-pde-audit--freeze-before-target — Target comparison requires a frozen branch

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The branch protocol requires geometry, gauge and boundary class, support basis, stability gates, source physics, extraction formulas, and source/protocol hashes to be fixed before target residuals are evaluated. Post-residual adjustment of support, normalization, boundary data, or moment-shape coefficients is excluded from a target-blind branch claim. The associated derivation explains why refitting adjustable coefficients after comparison destroys the test, but the complete freeze checklist is protocol-defined.

Sources:

- `research/pde_audit/notes/stage_v2_16_branch_freeze_no_refit_derivation.md` — heading `## 4. Why the no-refit rule is mandatory`
- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — heading `## Frozen Inputs`
- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — heading `## Non-Rescue Rules`

### source-pde-audit--handoff-validation — Solver packets are validated before evaluation

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The V2-22B validator checks solver-handoff schema and open-throat conditions, including grid and profile consistency, pre-target freeze, target blindness, gauge and boundary protocol, solver metadata, and target-leakage exclusions. V2-22A adapts accepted profile manifests to V2-21 coefficient packets, while V2-21 extracts frozen branch quantities and residuals. These are implemented guards and transformations; no readable recorded invocation and literal output chain in the packet promotes their run outcomes to measured evidence.

Sources:

- `research/pde_audit/scripts/stage_v2_22b_solver_handoff_validator.py` — function `validate_solver_output`
- `research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py` — function `adapt_manifest`
- `research/pde_audit/scripts/stage_v2_22a_profile_to_coefficient_adapter.py` — function `adapt_profile_branch_to_v21`
- `research/pde_audit/scripts/stage_v2_21_branch_extraction_fixture.py` — function `extract_branch`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`

### source-pde-audit--reduced-branch-result — Reduced branch prototype misses the target

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The interpretive records report that the minimal V2-23 open-throat linear-FEM branch is open and stable under its declared reduced gates but fails the full target packet. Separately, the V2-23 matrix row reports stable freeze identity, decreasing pairwise mesh-refinement deltas, and converging reduced solver values; the mesh audit implements those checks. This is a reduced prototype and an honest target miss, not a completed nonlinear PDE export or a falsification of every possible physical branch.

Sources:

- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- `research/pde_audit/notes/stage_v2_23_minimal_branch_solver_derivation.md` — heading `# Stage V2-23 — Minimal open-throat branch solver and first real residual extraction`
- `research/pde_audit/scripts/stage_v2_23_minimal_open_throat_branch_solver.py` — function `evaluate_gates`
- `research/pde_audit/scripts/stage_v2_23_mesh_convergence_audit.py` — function `main`

### source-pde-audit--target-blind-simulation-miss — Current frozen candidate families miss the target

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The paper-facing status records report that none of 192 reduced frozen candidates and none of three manufactured nonlinear candidates passed the target packet. This target-blind miss rejects those tested frozen families; it is distinct from the still-unbuilt physical exporter and does not exclude untested physical branches.

Sources:

- `research/pde_audit/notes/stage_v2_25_actual_branch_protocol_derivation.md` — heading `## 1. Status after the target-blind simulation miss`
- `research/pde_audit/notes/stage_v2_25_actual_branch_protocol_derivation.md` — heading `## Current claim status`
- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — heading `## 3. Why the clean algebra and failed simulations are not contradictory`
- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — heading `## Current status statement`

### source-pde-audit--physical-export-blocked — Full physical export remains open

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The program-status and protocol sources leave a strict nonlinear moving-throat exporter unbuilt. Before such export is permitted, the project must freeze or derive the parent throat dynamics, coupled physical residual equations, source and port laws, stability certificates, and same-branch coefficient extraction. Manufactured nonlinear readiness and reduced sweeps do not close this gate.

Sources:

- `research/pde_audit/notes/stage_v2_26_program_status_after_audit.md` — heading `# Stage V2-26 - Program Status After the PDE Audit`
- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — heading `## Required Packet`
- `research/pde_audit/simulation/README.md` — heading `# PDE Audit Simulation Bundle`

### source-pde-audit--material-closure-open — Superfluid material closure is incomplete

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

Candidate equations of state and sound-speed relations occur in reduced or frozen form, but the audited branch exporter does not yet determine density, phase/current, sound speed, any density-dependent effective light speed, open-system feedback, and the resulting audit coefficients on one self-consistent branch. Consequently these quantities remain branch inputs rather than jointly derived outputs.

Sources:

- `research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md` — heading `## 1. What is already present`
- `research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md` — heading `## 2. What is missing`
- `research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md` — heading `## 3. Why this can block the final PDE`

### source-pde-audit--electromagnetic-boundary — Electromagnetic ontology and gauge-localization boundary

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The current ontology contains a localized Maxwell field and mixed brane-bulk gauge-invariant components. Electric branch sign is assigned to the oriented puncture label, with localization dressing applied to observed charge, while circulation and fluxoid quantization belong to the distinct magnetic or vortical sector. For noncompact zero-mode reduction, the unweighted gauge-fixing term is acceptable only if Lorenz gauge is imposed before reduction; otherwise a finite localized gauge-fixing weight is required. These supported structures do not supply the still-missing closed recirculation law.

Sources:

- `research/pde_audit/notes/stage_v2_02_maxwell_gauge_localization_derivation.md` — heading `## 5. Stage verdict`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — heading `## 1. What EM structure is present`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — heading `## 2. Electric charge variables`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — heading `## 4. What is reduced versus still open`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — heading `## 6. Minimal claim boundary`

## Computed evidence represented by the source

### source-pde-audit--evidence-chain — Executable evidence boundary

Status: `lifecycle=current` · `evidence=provisional` · `memory_review=ai_draft`

The supplied sources define the top-level invocation `bash research/pde_audit/run_all.sh`, stage engines, fixture hashing, output capture, JSON-artifact generation, negative-control rejection paths, and combined referee-summary construction. The Mathematica runner invokes selected mirror scripts with `math -script`; those mirrors compare load-bearing algebra and fixture contracts but are not independent derivations. No tracked literal stdout or result transcript is included in the selected source files, so stage pass labels and reported residual outcomes remain provisional rather than measured.

Sources:

- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Reproducibility Artifacts`
- `research/pde_audit/run_all.sh` — nearby identifier `run_top_check`
- `research/pde_audit/scripts/run_all_audits.sh` — nearby identifier `OUTPUT_DIR`
- `research/pde_audit/mathematica/run_all_audits.sh` — nearby identifier `math -script`
- `research/pde_audit/scripts/stage_v2_24_negative_controls.py` — function `main`

## Assumptions, exclusions, and open questions

### source-pde-audit--scope-limitations — Conditional checks do not establish branch realization

Status: `lifecycle=current` · `evidence=open` · `memory_review=ai_draft`

The audit explicitly leaves universal compact-source strength, full nonlinear branch realization, negative-Krein or ghost sectors, material-sector completion, the future actual-branch exporter, and a closed electromagnetic recirculation/plumbing law unresolved. A stable target-failing reduced packet is a falsification of that frozen candidate, not a pipeline failure and not permission for post-hoc rescue.

Sources:

- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Summary`
- `research/pde_audit/CLAIM_CHECK_MATRIX.md` — heading `## Current Referee Gaps`
- `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md` — heading `## 4. What is reduced versus still open`
- `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md` — heading `## Status`

## Revision and supersession relationships

## Related topics and scripts

- `memory/sources/papers/moving-throat-pde.md`
- `memory/sources/software/stage1-solver.md`
- `memory/sources/software/stage1-verdict.md`
