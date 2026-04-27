---
id: MT_V2_22C_SMOKE_PIPELINE
title: V2-22C branch-realization smoke pipeline
type: pipeline_gate
layer: status_audit
status: implemented_control_and_smoke
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Runs end-to-end branch realization smoke pipeline, tolerance budget, and calibrated direct-coefficient control.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
outgoing_edges:
- target: MT_V2_23_OPEN_BRANCH_SOLVER
  relation: PRECEDES
  status: implemented
  note: Smoke pipeline precedes first real residual extraction.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/implemented_control_and_smoke
- topic/moving_throat
- type/pipeline_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-22C branch-realization smoke pipeline

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_22C_SMOKE_PIPELINE`  
> **Status:** `implemented_control_and_smoke`  
> **Layer:** `status_audit`  
> **Type:** `pipeline_gate`

## Summary

Runs end-to-end branch realization smoke pipeline, tolerance budget, and calibrated direct-coefficient control.

## Physical Meaning

Runs end-to-end branch realization smoke pipeline, tolerance budget, and calibrated direct-coefficient control.

## Mathematical Role

- Layer: `status_audit`
- Type: `pipeline_gate`
- Status: `implemented_control_and_smoke`

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
| `PRECEDES` | [[MT_V2_23_OPEN_BRANCH_SOLVER]] | Smoke pipeline precedes first real residual extraction. |

## Incoming Edges

- none

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
