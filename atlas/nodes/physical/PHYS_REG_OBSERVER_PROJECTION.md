---
id: PHYS_REG_OBSERVER_PROJECTION
title: Observer/projection physical spine
type: register
layer: physical_register
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Brane observation is projection with W(w), while brane effective laws are reductions under extra assumptions.
future_paper_needed: false
physical_ids:
- PHYS_BRANE_OBSERVER
math_ids:
- MATH_W_PROJECTION
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
outgoing_edges:
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: MATH_W_PROJECTION
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_BRANE_OBSERVER
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
- topic/projection
- type/register
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Observer/projection physical spine

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_REG_OBSERVER_PROJECTION`  
> **Status:** `active`  
> **Layer:** `physical_register`  
> **Type:** `register`

## Summary

Brane observation is projection with W(w), while brane effective laws are reductions under extra assumptions.

## Physical Meaning

Brane observation is projection with W(w), while brane effective laws are reductions under extra assumptions.

## Mathematical Role

- Layer: `physical_register`
- Type: `register`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_OBSERVER]]

### Related math nodes
- [[MATH_W_PROJECTION]]

### Related equations
- none

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LINKS_TO` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[MATH_W_PROJECTION]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_BRANE_OBSERVER]] | Physical register entry links to graph object. |

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
