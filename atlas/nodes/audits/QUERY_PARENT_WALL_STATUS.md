---
id: QUERY_PARENT_WALL_STATUS
title: What changes if S_eta or S_Sigma is not promoted?
type: query_test
layer: query_validation
status: v06_seed_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Confinement-only parent gives wall force/source but not autonomous wall PDE; parent promotion remains an open gate.
future_paper_needed: false
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
open_gate_ids:
- OPEN_PARENT_PROMOTION_S_SIGMA
outgoing_edges:
- target: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v06_seed_query
- topic/moving_throat
- type/query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# What changes if S_eta or S_Sigma is not promoted?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `QUERY_PARENT_WALL_STATUS`
> **Status:** `v06_seed_query`
> **Layer:** `query_validation`
> **Type:** `query_test`

## Summary

Confinement-only parent gives wall force/source but not autonomous wall PDE; parent promotion remains an open gate.

## Physical Meaning

Confinement-only parent gives wall force/source but not autonomous wall PDE; parent promotion remains an open gate.

## Mathematical Role

- Layer: `query_validation`
- Type: `query_test`
- Status: `v06_seed_query`
- Target: `OPEN_PARENT_PROMOTION_S_SIGMA`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- [[OPEN_PARENT_PROMOTION_S_SIGMA]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `EXPECTS_TARGET` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Query validation expected target node. |
| `STARTS_AT` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Query validation start node. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
