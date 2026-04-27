---
id: PHYS_RESPONSE_READOUTS
title: Compressed response readouts
type: readout_interpretation
layer: physical_ontology
status: not_physical_object
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: D0, C, P0, N2, N4 are low-order projected response readouts of the throat branch, not the throat itself.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
source_links:
- '[[FILE_PDE_AUDIT]]'
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_G2_COMMON_QUOTIENT
- CLAIM_PACKET_A_PACKET_B_SPLIT
- CLAIM_PARENT_WALL_STATUS_SPLIT
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE5_P0_TARGET
- CLAIM_STAGE6_FULL_BUNDLE_RATIO
status_firewall_ids:
- FIREWALL_READOUTS_NOT_THROAT
source_ids:
- FILE_PDE_AUDIT
outgoing_edges:
- target: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_bundle_open_branch
  note: Physical ontology object grounded by this claim.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_reduced_residual
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_branch_packets
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: GROUNDS_PHYSICAL_MEANING
  status: strict_parent_fail_effective_wall_pass
  note: Physical ontology object grounded by this claim.
- target: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: GROUNDS_PHYSICAL_MEANING
  status: paper_facing_ontology_discipline
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE5_P0_TARGET
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_grouped_bridge
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE6_FULL_BUNDLE_RATIO
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_bundle
  note: Physical ontology object grounded by this claim.
- target: READOUT_D0_C_P0_N2_N4
  relation: INTERPRETS
  status: status
  note: Low-order readouts are compressed projected quantities, not the object.
incoming_edges:
- source: FILE_PDE_AUDIT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: FIREWALL_READOUTS_NOT_THROAT
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/not_physical_object
- topic/moving_throat
- topic/quadrupole
- type/readout_interpretation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Compressed response readouts

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_RESPONSE_READOUTS`  
> **Status:** `not_physical_object`  
> **Layer:** `physical_ontology`  
> **Type:** `readout_interpretation`

## Summary

D0, C, P0, N2, N4 are low-order projected response readouts of the throat branch, not the throat itself.

## Physical Meaning

D0, C, P0, N2, N4 are low-order projected response readouts of the throat branch, not the throat itself.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `readout_interpretation`
- Status: `not_physical_object`

## Equation

$$
D0
$$

$$
C
$$

$$
P0
$$

$$
N2
$$

$$
N4
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_G2_COMMON_QUOTIENT]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE5_P0_TARGET]]
- [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_READOUTS_NOT_THROAT]]

### Source anchors
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_G2_COMMON_QUOTIENT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE5_P0_TARGET]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]] | Physical ontology object grounded by this claim. |
| `INTERPRETS` | [[READOUT_D0_C_P0_N2_N4]] | Low-order readouts are compressed projected quantities, not the object. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_PDE_AUDIT]] | File anchor documents this node. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_READOUTS_NOT_THROAT]] | Firewall preserves this correct status boundary. |

## Source Anchors

### Source anchor notes
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
