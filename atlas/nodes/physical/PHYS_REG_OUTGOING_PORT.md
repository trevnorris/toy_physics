---
id: PHYS_REG_OUTGOING_PORT
title: Outgoing port physical spine
type: register
layer: physical_register
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The radiative/outgoing route is a passive l=2 port coupled through the mixed sector and measured by the P0/N0/D0 response readout.
future_paper_needed: false
physical_ids:
- PHYS_OUTGOING_QUADRUPOLE_PORT
math_ids:
- MATH_COMPACT_L2_OUTGOING_FINGERPRINT
claim_ids:
- CLAIM_STAGE5_P0_TARGET
outgoing_edges:
- target: CLAIM_STAGE5_P0_TARGET
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: MATH_COMPACT_L2_OUTGOING_FINGERPRINT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_OUTGOING_QUADRUPOLE_PORT
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
- topic/quadrupole
- type/register
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Outgoing port physical spine

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_REG_OUTGOING_PORT`  
> **Status:** `active`  
> **Layer:** `physical_register`  
> **Type:** `register`

## Summary

The radiative/outgoing route is a passive l=2 port coupled through the mixed sector and measured by the P0/N0/D0 response readout.

## Physical Meaning

The radiative/outgoing route is a passive l=2 port coupled through the mixed sector and measured by the P0/N0/D0 response readout.

## Mathematical Role

- Layer: `physical_register`
- Type: `register`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_OUTGOING_QUADRUPOLE_PORT]]

### Related math nodes
- [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE5_P0_TARGET]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LINKS_TO` | [[CLAIM_STAGE5_P0_TARGET]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_OUTGOING_QUADRUPOLE_PORT]] | Physical register entry links to graph object. |

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
