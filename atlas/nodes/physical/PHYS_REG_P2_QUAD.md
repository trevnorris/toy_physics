---
id: PHYS_REG_P2_QUAD
title: Grouped P2/quadrupole physical spine
type: register
layer: physical_register
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Grouped P2 is a physical mouth/support shape-response bundle tied to the orbital/worldtube STF quadrupole, not a symbolic coefficient list.
future_paper_needed: false
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
math_ids:
- MATH_STF_SOURCE_MAP
claim_ids:
- CLAIM_STAGE7_O3_ISOTROPY
outgoing_edges:
- target: CLAIM_STAGE7_O3_ISOTROPY
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: MATH_STF_SOURCE_MAP
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
incoming_edges:
- source: VIEW_EXTENSION_WORKBENCH
  relation: USES_REGISTER
  status: active
  note: Physical register for extension workflow.
tags:
- atlas/node
- atlas/physical
- layer/physical_register
- status/active
- topic/quadrupole
- type/register
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Grouped P2/quadrupole physical spine

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_REG_P2_QUAD`  
> **Status:** `active`  
> **Layer:** `physical_register`  
> **Type:** `register`

## Summary

Grouped P2 is a physical mouth/support shape-response bundle tied to the orbital/worldtube STF quadrupole, not a symbolic coefficient list.

## Physical Meaning

Grouped P2 is a physical mouth/support shape-response bundle tied to the orbital/worldtube STF quadrupole, not a symbolic coefficient list.

## Mathematical Role

- Layer: `physical_register`
- Type: `register`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]

### Related math nodes
- [[MATH_STF_SOURCE_MAP]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE7_O3_ISOTROPY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LINKS_TO` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[MATH_STF_SOURCE_MAP]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Physical register entry links to graph object. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `USES_REGISTER` | [[VIEW_EXTENSION_WORKBENCH]] | Physical register for extension workflow. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
