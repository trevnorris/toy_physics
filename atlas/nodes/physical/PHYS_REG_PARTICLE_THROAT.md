---
id: PHYS_REG_PARTICLE_THROAT
title: Particle/throat physical spine
type: register
layer: physical_register
status: active
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Particle-like object = finite brane-bulk throat/mouth/interior/open-conduit system; math stack should not collapse it to a point too early.
future_paper_needed: false
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
math_ids:
- MATH_SIGMA_R_FIELD
claim_ids:
- CLAIM_STAGE1_GEOMETRY_LIFT
outgoing_edges:
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: MATH_SIGMA_R_FIELD
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- target: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
incoming_edges:
- source: WORKFLOW_STEP_1_PHYSICAL_OBJECT
  relation: MAY_START_WITH
  status: active
  note: Example physical register.
- source: VIEW_EXTENSION_WORKBENCH
  relation: USES_REGISTER
  status: active
  note: Physical register for extension workflow.
tags:
- atlas/node
- atlas/physical
- layer/physical_register
- status/active
- type/register
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Particle/throat physical spine

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_REG_PARTICLE_THROAT`
> **Status:** `active`
> **Layer:** `physical_register`
> **Type:** `register`

## Summary

Particle-like object = finite brane-bulk throat/mouth/interior/open-conduit system; math stack should not collapse it to a point too early.

## Physical Meaning

Particle-like object = finite brane-bulk throat/mouth/interior/open-conduit system; math stack should not collapse it to a point too early.

## Mathematical Role

- Layer: `physical_register`
- Type: `register`
- Status: `active`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]

### Related math nodes
- [[MATH_SIGMA_R_FIELD]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LINKS_TO` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[MATH_SIGMA_R_FIELD]] | Physical register entry links to graph object. |
| `LINKS_TO` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Physical register entry links to graph object. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `MAY_START_WITH` | [[WORKFLOW_STEP_1_PHYSICAL_OBJECT]] | Example physical register. |
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
