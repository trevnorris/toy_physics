---
id: MT_V2_28_ONTOLOGY_CHECKLIST
title: V2-28 physical picture checklist
type: ontology_status
layer: physical_ontology
status: paper_facing
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Defines finite brane-bulk throat/puncture picture and misreadings to avoid.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-28
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
outgoing_edges:
- target: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: ANCHORS
  status: paper-facing
  note: Physical picture anchors the atlas ontology layer.
- target: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: PATCHES_LANGUAGE
  status: paper_facing
  note: Paper-facing physical picture must be used consistently.
incoming_edges:
- source: READOUT_D0_C_P0_N2_N4
  relation: CONSTRAINED_BY
  status: status
  note: Checklist warns readouts are not physical ontology variables.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/paper_facing
- topic/moving_throat
- topic/quadrupole
- type/ontology_status
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-28 physical picture checklist

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_28_ONTOLOGY_CHECKLIST`
> **Status:** `paper_facing`
> **Layer:** `physical_ontology`
> **Type:** `ontology_status`

## Summary

Defines finite brane-bulk throat/puncture picture and misreadings to avoid.

## Physical Meaning

Defines finite brane-bulk throat/puncture picture and misreadings to avoid.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `ontology_status`
- Status: `paper_facing`

## Equation

$$
D0,C,P0,N2,N4 are response readouts, not ontology variables
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Physical picture anchors the atlas ontology layer. |
| `PATCHES_LANGUAGE` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Paper-facing physical picture must be used consistently. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONSTRAINED_BY` | [[READOUT_D0_C_P0_N2_N4]] | Checklist warns readouts are not physical ontology variables. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
