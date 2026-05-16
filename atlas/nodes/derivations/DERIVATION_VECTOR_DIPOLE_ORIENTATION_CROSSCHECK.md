---
id: DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK
title: Vector dipole orientation law and coaxial recovery
type: symbolic_audit
layer: derivation
status: exact_dipole_far_field_crosscheck
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Step 05 builds m_A=pi*R_A^2*I_A*s_A as actual 3-vectors, computes m1.m2 and m_i.dhat from components, recovers the Step-03 coaxial leading force, and treats finite-size terms as...
future_paper_needed: false
source_files:
- notes/circulation/step_05_3d_orientation_and_finite_size.md
- notes/circulation/step_05_3d_orientation_and_finite_size_sympy.py
- notes/circulation/step_05_3d_orientation_and_finite_size_output.txt
source_links:
- '[[FILE_CIRCULATION_PACKAGE]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
claim_ids:
- CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
source_ids:
- FILE_CIRCULATION_PACKAGE
outgoing_edges:
- target: DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP
  relation: RECOVERS_COAXIAL_LIMIT
  status: exact_dipole_far_field_crosscheck
  note: Vector dipole audit recovers the Step-03 leading coaxial force.
- target: CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
  relation: SUPPLIES_ORIENTATION_LAW
  status: exact_dipole_far_field_crosscheck
  note: 3D vector dipole orientation law supplies the geometric sign rule.
incoming_edges:
- source: FILE_CIRCULATION_PACKAGE
  relation: DOCUMENTS
  status: exact_dipole_far_field_crosscheck
  note: Circulation package Step 05 documents the vector dipole orientation cross-check.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/exact_dipole_far_field_crosscheck
- type/symbolic_audit
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Vector dipole orientation law and coaxial recovery

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK`
> **Status:** `exact_dipole_far_field_crosscheck`
> **Layer:** `derivation`
> **Type:** `symbolic_audit`

## Summary

Step 05 builds m_A=pi*R_A^2*I_A*s_A as actual 3-vectors, computes m1.m2 and m_i.dhat from components, recovers the Step-03 coaxial leading force, and treats finite-size terms as asymptotic.

## Physical Meaning

Step 05 builds m_A=pi*R_A^2*I_A*s_A as actual 3-vectors, computes m1.m2 and m_i.dhat from components, recovers the Step-03 coaxial leading force, and treats finite-size terms as asymptotic.

## Mathematical Role

- Layer: `derivation`
- Type: `symbolic_audit`
- Status: `exact_dipole_far_field_crosscheck`

## Equation

$$
U=-mu0/(4*pi*d^3)*(3*(m1.dhat)*(m2.dhat)-m1.m2)
$$

$$
F_d=-dU/dd
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_CIRCULATION_PACKAGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `RECOVERS_COAXIAL_LIMIT` | [[DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP]] | Vector dipole audit recovers the Step-03 leading coaxial force. |
| `SUPPLIES_ORIENTATION_LAW` | [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]] | 3D vector dipole orientation law supplies the geometric sign rule. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_CIRCULATION_PACKAGE]] | Circulation package Step 05 documents the vector dipole orientation cross-check. |

## Source Anchors

### Source anchor notes
- [[FILE_CIRCULATION_PACKAGE]]

### Source files
- `notes/circulation/step_05_3d_orientation_and_finite_size.md`
- `notes/circulation/step_05_3d_orientation_and_finite_size_sympy.py`
- `notes/circulation/step_05_3d_orientation_and_finite_size_output.txt`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
