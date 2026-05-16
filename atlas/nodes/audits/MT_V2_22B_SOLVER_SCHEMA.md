---
id: MT_V2_22B_SOLVER_SCHEMA
title: V2-22B numerical branch profile schema
type: solver_schema_gate
layer: status_audit
status: implemented
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Defines solver export schema, validation gates, and conversion into profile manifest; rejects hard-cap negative control.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
outgoing_edges:
- target: MT_V2_22A_PROFILE_ADAPTER
  relation: FEEDS
  status: implemented
  note: Solver export schema validates and feeds profile adapter.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/implemented
- topic/moving_throat
- type/solver_schema_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-22B numerical branch profile schema

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_22B_SOLVER_SCHEMA`
> **Status:** `implemented`
> **Layer:** `status_audit`
> **Type:** `solver_schema_gate`

## Summary

Defines solver export schema, validation gates, and conversion into profile manifest; rejects hard-cap negative control.

## Physical Meaning

Defines solver export schema, validation gates, and conversion into profile manifest; rejects hard-cap negative control.

## Mathematical Role

- Layer: `status_audit`
- Type: `solver_schema_gate`
- Status: `implemented`

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
| `FEEDS` | [[MT_V2_22A_PROFILE_ADAPTER]] | Solver export schema validates and feeds profile adapter. |

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
