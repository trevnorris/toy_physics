---
id: WORKFLOW_STEP_2_MATH_REPRESENTATION
title: 2. Identify math representation
type: workflow_step
layer: extension_workflow_step
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Attach the target to its math objects and equation anchors before adding new derivation nodes.
future_paper_needed: false
outgoing_edges:
- target: WORKFLOW_STEP_3_STATUS_CHECK
  relation: NEXT_STEP
  status: active
  note: Extension workflow order.
incoming_edges:
- source: VIEW_EXTENSION_WORKBENCH
  relation: CONTAINS_STEP
  status: active
  note: Future extension workflow step.
- source: WORKFLOW_STEP_1_PHYSICAL_OBJECT
  relation: NEXT_STEP
  status: active
  note: Extension workflow order.
tags:
- atlas/extensions
- atlas/node
- layer/extension_workflow_step
- status/active
- type/workflow_step
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 2. Identify math representation

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `WORKFLOW_STEP_2_MATH_REPRESENTATION`
> **Status:** `active`
> **Layer:** `extension_workflow_step`
> **Type:** `workflow_step`

## Summary

Attach the target to its math objects and equation anchors before adding new derivation nodes.

## Physical Meaning

Attach the target to its math objects and equation anchors before adding new derivation nodes.

## Mathematical Role

- Layer: `extension_workflow_step`
- Type: `workflow_step`
- Status: `active`

## Atlas Links

### Related physical nodes
- none

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
| `NEXT_STEP` | [[WORKFLOW_STEP_3_STATUS_CHECK]] | Extension workflow order. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_STEP` | [[VIEW_EXTENSION_WORKBENCH]] | Future extension workflow step. |
| `NEXT_STEP` | [[WORKFLOW_STEP_1_PHYSICAL_OBJECT]] | Extension workflow order. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
