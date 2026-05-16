---
id: PHYS_REG_CHARGE_EM
title: Charge/EM physical spine
type: register
layer: physical_register
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Electric sign, brane charge magnitude, magnetic/vortical circulation, and mixed-core transport are separate objects.
future_paper_needed: false
physical_ids:
- PHYS_CHARGE_BRANCH
- PHYS_MIXED_EM_CORE
math_ids:
- MATH_QSTAR_QEFF
claim_ids:
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
outgoing_edges:
- target: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: MATH_QSTAR_QEFF
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_CHARGE_BRANCH
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_MIXED_EM_CORE
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
- topic/charge
- topic/maxwell
- topic/moving_throat
- type/register
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Charge/EM physical spine

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_REG_CHARGE_EM`
> **Status:** `active`
> **Layer:** `physical_register`
> **Type:** `register`

## Summary

Electric sign, brane charge magnitude, magnetic/vortical circulation, and mixed-core transport are separate objects.

## Physical Meaning

Electric sign, brane charge magnitude, magnetic/vortical circulation, and mixed-core transport are separate objects.

## Mathematical Role

- Layer: `physical_register`
- Type: `register`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_CHARGE_BRANCH]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- [[MATH_QSTAR_QEFF]]

### Related equations
- none

### Related claims
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LINKS_TO` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[MATH_QSTAR_QEFF]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_CHARGE_BRANCH]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_MIXED_EM_CORE]] | Physical register entry links to graph object. |

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
