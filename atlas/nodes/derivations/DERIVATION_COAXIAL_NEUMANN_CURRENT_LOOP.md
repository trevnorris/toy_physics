---
id: DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP
title: Coaxial current-loop mutual inductance finite-mouth expansion
type: symbolic_audit
layer: derivation
status: exact_within_fixed_current_far_field_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Step 03 derives the reduced Neumann integral from 3-vector loop parameterizations and symbolically expands the coaxial fixed-current mutual inductance and force through the disp...
future_paper_needed: false
source_files:
- notes/circulation/step_03_current_loop_closure_mutual_inductance.md
- notes/circulation/step_03_current_loop_closure_mutual_inductance_sympy.py
- notes/circulation/step_03_current_loop_closure_mutual_inductance_output.txt
source_links:
- '[[FILE_CIRCULATION_PACKAGE]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
claim_ids:
- CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
- CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
source_ids:
- FILE_CIRCULATION_PACKAGE
outgoing_edges:
- target: CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
  relation: SUPPLIES_CLOSURE_EXAMPLE
  status: exact_within_fixed_current_far_field_closure
  note: Fixed-current loop law is a named closure, not a consequence of fluxoid alone.
- target: CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
  relation: SUPPLIES_FIXED_CURRENT_FORCE
  status: exact_within_fixed_current_far_field_closure
  note: Facing-mouth sign table imports the fixed-current coaxial force sign.
incoming_edges:
- source: FILE_CIRCULATION_PACKAGE
  relation: DOCUMENTS
  status: exact_within_fixed_current_far_field_closure
  note: Circulation package Step 03 documents the fixed-current Neumann loop derivation.
- source: DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK
  relation: RECOVERS_COAXIAL_LIMIT
  status: exact_dipole_far_field_crosscheck
  note: Vector dipole audit recovers the Step-03 leading coaxial force.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/exact_within_fixed_current_far_field_closure
- type/symbolic_audit
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Coaxial current-loop mutual inductance finite-mouth expansion

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP`
> **Status:** `exact_within_fixed_current_far_field_closure`
> **Layer:** `derivation`
> **Type:** `symbolic_audit`

## Summary

Step 03 derives the reduced Neumann integral from 3-vector loop parameterizations and symbolically expands the coaxial fixed-current mutual inductance and force through the displayed finite-mouth terms.

## Physical Meaning

Step 03 derives the reduced Neumann integral from 3-vector loop parameterizations and symbolically expands the coaxial fixed-current mutual inductance and force through the displayed finite-mouth terms.

## Mathematical Role

- Layer: `derivation`
- Type: `symbolic_audit`
- Status: `exact_within_fixed_current_far_field_closure`

## Equation

$$
M(d)=mu0*pi*R1^2*R2^2/(2*d^3)*(1-3*(R1^2+R2^2)/(2*d^2)+15*(R1^4+3*R1^2*R2^2+R2^4)/(8*d^4)+...)
$$

$$
F_d=-3*mu0*pi*I1*I2*R1^2*R2^2/(2*d^4)*(1-5*(R1^2+R2^2)/(2*d^2)+35*(R1^4+3*R1^2*R2^2+R2^4)/(8*d^4)+...)
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
- [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_CIRCULATION_PACKAGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `SUPPLIES_CLOSURE_EXAMPLE` | [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]] | Fixed-current loop law is a named closure, not a consequence of fluxoid alone. |
| `SUPPLIES_FIXED_CURRENT_FORCE` | [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]] | Facing-mouth sign table imports the fixed-current coaxial force sign. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_CIRCULATION_PACKAGE]] | Circulation package Step 03 documents the fixed-current Neumann loop derivation. |
| `RECOVERS_COAXIAL_LIMIT` | [[DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK]] | Vector dipole audit recovers the Step-03 leading coaxial force. |

## Source Anchors

### Source anchor notes
- [[FILE_CIRCULATION_PACKAGE]]

### Source files
- `notes/circulation/step_03_current_loop_closure_mutual_inductance.md`
- `notes/circulation/step_03_current_loop_closure_mutual_inductance_sympy.py`
- `notes/circulation/step_03_current_loop_closure_mutual_inductance_output.txt`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
