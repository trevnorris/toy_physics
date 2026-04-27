---
id: MT_V2_25_ACTUAL_BRANCH_PROTOCOL
title: V2-25 actual branch protocol
type: protocol_gate
layer: status_audit
status: open_actual_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Defines post-miss coefficient map, outgoing prefactor, moment-shape conditions, Packet A/B, notes intake, and no-refit rules.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: REFINES
  status: open
  note: Defines Packet A/B and no-refit actual branch protocol.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/open_actual_branch
- topic/moving_throat
- type/protocol_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-25 actual branch protocol

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_25_ACTUAL_BRANCH_PROTOCOL`  
> **Status:** `open_actual_branch`  
> **Layer:** `status_audit`  
> **Type:** `protocol_gate`

## Summary

Defines post-miss coefficient map, outgoing prefactor, moment-shape conditions, Packet A/B, notes intake, and no-refit rules.

## Physical Meaning

Defines post-miss coefficient map, outgoing prefactor, moment-shape conditions, Packet A/B, notes intake, and no-refit rules.

## Mathematical Role

- Layer: `status_audit`
- Type: `protocol_gate`
- Status: `open_actual_branch`

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
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REFINES` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Defines Packet A/B and no-refit actual branch protocol. |

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
