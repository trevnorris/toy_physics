---
id: FILE_CIRCULATION_PACKAGE
title: notes/circulation derivation package
type: source_file
layer: file_anchor
status: note_audit_package
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Six-step circulation pair-force audit package covering fluxoid winding, no-universal-force status, coaxial current-loop closure, facing-mouth swirl labels, vector dipole orienta...
source_kind: future_paper_note
future_paper_needed: true
source_files:
- notes/circulation/README.md
- notes/circulation/step_01_fluxoid_firewall.md
- notes/circulation/step_01_fluxoid_firewall_sympy.py
- notes/circulation/step_01_fluxoid_firewall_output.txt
- notes/circulation/step_02_no_universal_force.md
- notes/circulation/step_02_no_universal_force_sympy.py
- notes/circulation/step_02_no_universal_force_output.txt
- notes/circulation/step_03_current_loop_closure_mutual_inductance.md
- notes/circulation/step_03_current_loop_closure_mutual_inductance_sympy.py
- notes/circulation/step_03_current_loop_closure_mutual_inductance_output.txt
- notes/circulation/step_04_facing_mouth_swirl_sign.md
- notes/circulation/step_04_facing_mouth_swirl_sign_sympy.py
- notes/circulation/step_04_facing_mouth_swirl_sign_output.txt
- notes/circulation/step_05_3d_orientation_and_finite_size.md
- notes/circulation/step_05_3d_orientation_and_finite_size_sympy.py
- notes/circulation/step_05_3d_orientation_and_finite_size_output.txt
- notes/circulation/step_06_mixed_sector_plumbing_closure.md
- notes/circulation/step_06_mixed_sector_plumbing_closure_sympy.py
- notes/circulation/step_06_mixed_sector_plumbing_closure_output.txt
claim_ids:
- CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
- CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION
- CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
- CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
outgoing_edges:
- target: DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP
  relation: DOCUMENTS
  status: exact_within_fixed_current_far_field_closure
  note: Circulation package Step 03 documents the fixed-current Neumann loop derivation.
- target: DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK
  relation: DOCUMENTS
  status: exact_dipole_far_field_crosscheck
  note: Circulation package Step 05 documents the vector dipole orientation cross-check.
- target: CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_fixed_current_orientation_closure
  note: Circulation package Step 04 anchors the facing-mouth local-swirl sign rule.
- target: CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_identity_audit
  note: Circulation package Step 01 anchors this quantized phase-winding audit.
- target: CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: conditional_open_plumbing
  note: Circulation package Step 06 anchors the Lambda_A-conditional mixed plumbing status.
- target: CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_negative_within_3d_closure_analysis
  note: Circulation package Step 02 anchors the no-universal-force scope result.
tags:
- atlas/node
- atlas/sources
- layer/file_anchor
- status/note_audit_package
- topic/maxwell
- topic/moving_throat
- type/source_file
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# notes/circulation derivation package

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `FILE_CIRCULATION_PACKAGE`
> **Status:** `note_audit_package`
> **Layer:** `file_anchor`
> **Type:** `source_file`

## Summary

Six-step circulation pair-force audit package covering fluxoid winding, no-universal-force status, coaxial current-loop closure, facing-mouth swirl labels, vector dipole orientation, and mixed-sector plumbing status.

> [!note] Future paper needed
> This node is intentionally anchored to a notes/derivation source until a maintained paper exists.

## Physical Meaning

Six-step circulation pair-force audit package covering fluxoid winding, no-universal-force status, coaxial current-loop closure, facing-mouth swirl labels, vector dipole orientation, and mixed-sector plumbing status.

## Mathematical Role

- Layer: `file_anchor`
- Type: `source_file`
- Status: `note_audit_package`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]]
- [[CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION]]
- [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]]
- [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP]] | Circulation package Step 03 documents the fixed-current Neumann loop derivation. |
| `DOCUMENTS` | [[DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK]] | Circulation package Step 05 documents the vector dipole orientation cross-check. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]] | Circulation package Step 04 anchors the facing-mouth local-swirl sign rule. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION]] | Circulation package Step 01 anchors this quantized phase-winding audit. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]] | Circulation package Step 06 anchors the Lambda_A-conditional mixed plumbing status. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]] | Circulation package Step 02 anchors the no-universal-force scope result. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/circulation/README.md`
- `notes/circulation/step_01_fluxoid_firewall.md`
- `notes/circulation/step_01_fluxoid_firewall_sympy.py`
- `notes/circulation/step_01_fluxoid_firewall_output.txt`
- `notes/circulation/step_02_no_universal_force.md`
- `notes/circulation/step_02_no_universal_force_sympy.py`
- `notes/circulation/step_02_no_universal_force_output.txt`
- `notes/circulation/step_03_current_loop_closure_mutual_inductance.md`
- `notes/circulation/step_03_current_loop_closure_mutual_inductance_sympy.py`
- `notes/circulation/step_03_current_loop_closure_mutual_inductance_output.txt`
- `notes/circulation/step_04_facing_mouth_swirl_sign.md`
- `notes/circulation/step_04_facing_mouth_swirl_sign_sympy.py`
- `notes/circulation/step_04_facing_mouth_swirl_sign_output.txt`
- `notes/circulation/step_05_3d_orientation_and_finite_size.md`
- `notes/circulation/step_05_3d_orientation_and_finite_size_sympy.py`
- `notes/circulation/step_05_3d_orientation_and_finite_size_output.txt`
- `notes/circulation/step_06_mixed_sector_plumbing_closure.md`
- `notes/circulation/step_06_mixed_sector_plumbing_closure_sympy.py`
- `notes/circulation/step_06_mixed_sector_plumbing_closure_output.txt`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
