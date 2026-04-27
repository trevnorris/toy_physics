---
id: PHYS_REG_BRANCH_EXPORT
title: Actual branch physical spine
type: register
layer: physical_register
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Final theorem work requires a target-blind branch exporter that returns physical profile, material, support, mixed, and readout packets before target comparison.
future_paper_needed: false
physical_ids:
- PHYS_BRANCH_EXPORTER
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
outgoing_edges:
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_BRANCH_EXPORTER
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
- topic/maxwell
- type/register
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Actual branch physical spine

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_REG_BRANCH_EXPORT`  
> **Status:** `active`  
> **Layer:** `physical_register`  
> **Type:** `register`

## Summary

Final theorem work requires a target-blind branch exporter that returns physical profile, material, support, mixed, and readout packets before target comparison.

## Physical Meaning

Final theorem work requires a target-blind branch exporter that returns physical profile, material, support, mixed, and readout packets before target comparison.

## Mathematical Role

- Layer: `physical_register`
- Type: `register`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANCH_EXPORTER]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LINKS_TO` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_BRANCH_EXPORTER]] | Physical register entry links to graph object. |

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
