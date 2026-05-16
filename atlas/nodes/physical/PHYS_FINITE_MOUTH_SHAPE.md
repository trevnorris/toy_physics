---
id: PHYS_FINITE_MOUTH_SHAPE
title: Finite-mouth shape physics
type: geometry_response
layer: physical_ontology
status: reduced_conditional
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Finite throat size makes Hessian/tidal mouth loading and headless P22 shape physics relevant.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/atom_work.md
- pde_audit_full.md
- atom_work.md
legacy_sources:
- pde_audit_full.md:V2-28
- atom_work.md
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_PACKET_A_PACKET_B_SPLIT
- CLAIM_PARENT_WALL_STATUS_SPLIT
- CLAIM_STAGE024_O3_ISOTROPY
- CLAIM_STAGE025_031_SELECTED_BRANCH
- CLAIM_STAGE1_GEOMETRY_LIFT
open_gate_ids:
- TARGET_PACKET_B
outgoing_edges:
- target: TARGET_PACKET_B
  relation: AFFECTS
  status: reduced/open
  note: P22/orbit-lock/mouth bracing enters weak-axisymmetric packet.
- target: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: FEEDS
  status: reduced
  note: Finite-mouth Hessian/tidal loading drives P0/P2 response.
- target: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_bundle_open_branch
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_branch_packets
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: GROUNDS_PHYSICAL_MEANING
  status: strict_parent_fail_effective_wall_pass
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE024_O3_ISOTROPY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_angular_reduced
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE025_031_SELECTED_BRANCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_selected_reduced_branch
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: GROUNDS_PHYSICAL_MEANING
  status: effective_geometry_lift
  note: Physical ontology object grounded by this claim.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/reduced_conditional
- topic/atom
- topic/moving_throat
- topic/quadrupole
- type/geometry_response
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Finite-mouth shape physics

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_FINITE_MOUTH_SHAPE`
> **Status:** `reduced_conditional`
> **Layer:** `physical_ontology`
> **Type:** `geometry_response`

## Summary

Finite throat size makes Hessian/tidal mouth loading and headless P22 shape physics relevant.

## Physical Meaning

Finite throat size makes Hessian/tidal mouth loading and headless P22 shape physics relevant.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `geometry_response`
- Status: `reduced_conditional`

## Equation

$$
P0 trace channel
$$

$$
P2 traceless channel
$$

$$
P22 mouth ellipse
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]
- [[CLAIM_STAGE024_O3_ISOTROPY]]
- [[CLAIM_STAGE025_031_SELECTED_BRANCH]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- [[TARGET_PACKET_B]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `AFFECTS` | [[TARGET_PACKET_B]] | P22/orbit-lock/mouth bracing enters weak-axisymmetric packet. |
| `FEEDS` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Finite-mouth Hessian/tidal loading drives P0/P2 response. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE024_O3_ISOTROPY]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE025_031_SELECTED_BRANCH]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Physical ontology object grounded by this claim. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/atom_work.md`
- `pde_audit_full.md`
- `atom_work.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
