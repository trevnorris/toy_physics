---
id: WORKFLOW_STEP_1_PHYSICAL_OBJECT
title: 1. Identify physical object
type: workflow_step
layer: extension_workflow_step
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: 'Start with the physical object: throat, projection, charge, P2 response, mixed core, outgoing port, material branch, or lepton/anomaly object.'
future_paper_needed: false
physical_ids:
- PHYS_REG_PARTICLE_THROAT
outgoing_edges:
- target: PHYS_REG_PARTICLE_THROAT
  relation: MAY_START_WITH
  status: active
  note: Example physical register.
- target: WORKFLOW_STEP_2_MATH_REPRESENTATION
  relation: NEXT_STEP
  status: active
  note: Extension workflow order.
incoming_edges:
- source: VIEW_EXTENSION_WORKBENCH
  relation: CONTAINS_STEP
  status: active
  note: Future extension workflow step.
tags:
- atlas/extensions
- atlas/node
- layer/extension_workflow_step
- status/active
- topic/charge
- topic/g2
- topic/lepton
- topic/maxwell
- topic/projection
- type/workflow_step
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 1. Identify physical object

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `WORKFLOW_STEP_1_PHYSICAL_OBJECT`
> **Status:** `active`
> **Layer:** `extension_workflow_step`
> **Type:** `workflow_step`

## Summary

Start with the physical object: throat, projection, charge, P2 response, mixed core, outgoing port, material branch, or lepton/anomaly object.

## Physical Meaning

Start with the physical object: throat, projection, charge, P2 response, mixed core, outgoing port, material branch, or lepton/anomaly object.

## Mathematical Role

- Layer: `extension_workflow_step`
- Type: `workflow_step`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_REG_PARTICLE_THROAT]]

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
| `MAY_START_WITH` | [[PHYS_REG_PARTICLE_THROAT]] | Example physical register. |
| `NEXT_STEP` | [[WORKFLOW_STEP_2_MATH_REPRESENTATION]] | Extension workflow order. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_STEP` | [[VIEW_EXTENSION_WORKBENCH]] | Future extension workflow step. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
