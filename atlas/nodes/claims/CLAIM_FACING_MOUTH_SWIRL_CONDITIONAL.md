---
id: CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
title: Facing-mouth opposite local swirl attraction is closure-conditional
type: claim
layer: claim_theorem
status: exact_within_fixed_current_orientation_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: With a fixed-current/current-like coaxial closure, a real tangential swirl circulation, and oriented area vectors, facing local normals give sigma1*sigma2=-1, so opposite local ...
future_paper_needed: false
source_links:
- '[[FILE_CIRCULATION_PACKAGE]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
claim_ids:
- CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
- CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
source_ids:
- FILE_CIRCULATION_PACKAGE
outgoing_edges:
- target: PHYS_MAGNETIC_VORTICAL_CIRCULATION
  relation: MAPS_LOCAL_LABEL_TO_GLOBAL_AXIS
  status: exact_within_fixed_current_orientation_closure
  note: Oriented loop area maps local swirl labels to the global current/dipole axis.
- target: CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
  relation: RESPECTS_SCOPE_LIMIT
  status: exact_within_fixed_current_orientation_closure
  note: Facing-mouth attraction statement remains closure-conditional and does not override the no-universal-force result.
incoming_edges:
- source: CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
  relation: IMPORTS_SIGN_TABLE_CONDITIONALLY
  status: conditional_open_plumbing
  note: Mixed-sector statement imports the facing-mouth sign rule only after a current-like plumbing map is assumed.
- source: FILE_CIRCULATION_PACKAGE
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_fixed_current_orientation_closure
  note: Circulation package Step 04 anchors the facing-mouth local-swirl sign rule.
- source: DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP
  relation: SUPPLIES_FIXED_CURRENT_FORCE
  status: exact_within_fixed_current_far_field_closure
  note: Facing-mouth sign table imports the fixed-current coaxial force sign.
- source: DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK
  relation: SUPPLIES_ORIENTATION_LAW
  status: exact_dipole_far_field_crosscheck
  note: 3D vector dipole orientation law supplies the geometric sign rule.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_within_fixed_current_orientation_closure
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Facing-mouth opposite local swirl attraction is closure-conditional

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL`
> **Status:** `exact_within_fixed_current_orientation_closure`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

With a fixed-current/current-like coaxial closure, a real tangential swirl circulation, and oriented area vectors, facing local normals give sigma1*sigma2=-1, so opposite local swirl labels attract; this is not a fluxoid-only theorem.

## Claim

With a fixed-current/current-like coaxial closure, a real tangential swirl circulation, and oriented area vectors, facing local normals give sigma1*sigma2=-1, so opposite local swirl labels attract; this is not a fluxoid-only theorem.

## Physical Meaning

With a fixed-current/current-like coaxial closure, a real tangential swirl circulation, and oriented area vectors, facing local normals give sigma1*sigma2=-1, so opposite local swirl labels attract; this is not a fluxoid-only theorem.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_within_fixed_current_orientation_closure`
- Outputs: `CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID`

## Atlas Links

### Related physical nodes
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]]
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
| `MAPS_LOCAL_LABEL_TO_GLOBAL_AXIS` | [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]] | Oriented loop area maps local swirl labels to the global current/dipole axis. |
| `RESPECTS_SCOPE_LIMIT` | [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]] | Facing-mouth attraction statement remains closure-conditional and does not override the no-universal-force ... |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `IMPORTS_SIGN_TABLE_CONDITIONALLY` | [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]] | Mixed-sector statement imports the facing-mouth sign rule only after a current-like plumbing map is assumed. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_CIRCULATION_PACKAGE]] | Circulation package Step 04 anchors the facing-mouth local-swirl sign rule. |
| `SUPPLIES_FIXED_CURRENT_FORCE` | [[DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP]] | Facing-mouth sign table imports the fixed-current coaxial force sign. |
| `SUPPLIES_ORIENTATION_LAW` | [[DERIVATION_VECTOR_DIPOLE_ORIENTATION_CROSSCHECK]] | 3D vector dipole orientation law supplies the geometric sign rule. |

## Source Anchors

### Source anchor notes
- [[FILE_CIRCULATION_PACKAGE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
