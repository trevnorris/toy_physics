---
id: PHYS_MOUTH_CROSS_SECTION
title: Brane-side mouth
type: geometry_feature
layer: physical_ontology
status: physical_ontology
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The brane entrance cross-section of the throat; not the entire defect.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_notes_full.md
- pde_audit_full.md
- moving_throat_output_full.md
legacy_sources:
- pde_audit_full.md:V2-28
- moving_throat_output_full.md
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
- PHYS_MOUTH_INTERIOR_FIREWALL
claim_ids:
- CLAIM_STAGE2_AL_RECOVERY
status_firewall_ids:
- PHYS_MOUTH_INTERIOR_FIREWALL
outgoing_edges:
- target: CLAIM_STAGE2_AL_RECOVERY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_wall_action
  note: Physical ontology object grounded by this claim.
- target: PHYS_MOUTH_INTERIOR_FIREWALL
  relation: REQUIRES_FIREWALL
  status: paper_facing
  note: Do not identify the mouth with the whole defect.
incoming_edges:
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: HAS_PART
  status: ontology
  note: Mouth is brane entrance cross-section.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: HAS_PART
  status: physical
  note: Mouth is the brane-side cross-section.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/physical_ontology
- topic/moving_throat
- type/geometry_feature
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Brane-side mouth

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_MOUTH_CROSS_SECTION`  
> **Status:** `physical_ontology`  
> **Layer:** `physical_ontology`  
> **Type:** `geometry_feature`

## Summary

The brane entrance cross-section of the throat; not the entire defect.

## Physical Meaning

The brane entrance cross-section of the throat; not the entire defect.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `geometry_feature`
- Status: `physical_ontology`

## Equation

$$
R(Omega,0,t)
$$

$$
a(t)=average_S2 R(Omega,0,t)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]
- [[PHYS_MOUTH_INTERIOR_FIREWALL]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_STAGE2_AL_RECOVERY]]

### Open gates
- none

### Status firewalls
- [[PHYS_MOUTH_INTERIOR_FIREWALL]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE2_AL_RECOVERY]] | Physical ontology object grounded by this claim. |
| `REQUIRES_FIREWALL` | [[PHYS_MOUTH_INTERIOR_FIREWALL]] | Do not identify the mouth with the whole defect. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `HAS_PART` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Mouth is brane entrance cross-section. |
| `HAS_PART` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Mouth is the brane-side cross-section. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_notes_full.md`
- `pde_audit_full.md`
- `moving_throat_output_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
