---
id: ATLAS_RULE_PHYS_TO_MATH
title: 'Atlas rule: every equation has physical meaning'
type: graph_rule
layer: atlas_meta
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Every equation node should link to a physical object, and every physical object should link to its mathematical representation.
future_paper_needed: false
source_files:
- derived_from_current_task
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
outgoing_edges:
- target: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: APPLIES_TO
  status: active
  note: Start with physical defect, then map to Sigma/R, wall action, response readouts.
- target: VIEW_CLAIM_LAYER
  relation: GOVERNS
  status: active
  note: Every claim node must link physical meaning, math representation, proof status, and file anchor.
tags:
- atlas/meta
- atlas/node
- layer/atlas_meta
- status/active
- type/graph_rule
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Atlas rule: every equation has physical meaning

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `ATLAS_RULE_PHYS_TO_MATH`  
> **Status:** `active`  
> **Layer:** `atlas_meta`  
> **Type:** `graph_rule`

## Summary

Every equation node should link to a physical object, and every physical object should link to its mathematical representation.

## Physical Meaning

Every equation node should link to a physical object, and every physical object should link to its mathematical representation.

## Mathematical Role

- Layer: `atlas_meta`
- Type: `graph_rule`
- Status: `active`

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
| `APPLIES_TO` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Start with physical defect, then map to Sigma/R, wall action, response readouts. |
| `GOVERNS` | [[VIEW_CLAIM_LAYER]] | Every claim node must link physical meaning, math representation, proof status, and file anchor. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `derived_from_current_task`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
