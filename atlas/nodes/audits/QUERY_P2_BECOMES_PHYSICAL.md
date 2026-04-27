---
id: QUERY_P2_BECOMES_PHYSICAL
title: Where does grouped P2 become physical geometry?
type: query_test
layer: query_validation
status: v06_seed_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Stage 1 makes grouped P2 a real wall/support mode; later stages add BdG, mixed/outgoing bridge, full bundle, and O(3) isotropy.
future_paper_needed: false
claim_ids:
- CLAIM_STAGE1_GEOMETRY_LIFT
- CLAIM_STAGE7_O3_ISOTROPY
outgoing_edges:
- target: CLAIM_STAGE7_O3_ISOTROPY
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v06_seed_query
- topic/maxwell
- topic/moving_throat
- type/query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Where does grouped P2 become physical geometry?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `QUERY_P2_BECOMES_PHYSICAL`  
> **Status:** `v06_seed_query`  
> **Layer:** `query_validation`  
> **Type:** `query_test`

## Summary

Stage 1 makes grouped P2 a real wall/support mode; later stages add BdG, mixed/outgoing bridge, full bundle, and O(3) isotropy.

## Physical Meaning

Stage 1 makes grouped P2 a real wall/support mode; later stages add BdG, mixed/outgoing bridge, full bundle, and O(3) isotropy.

## Mathematical Role

- Layer: `query_validation`
- Type: `query_test`
- Status: `v06_seed_query`
- Target: `CLAIM_STAGE7_O3_ISOTROPY`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]
- [[CLAIM_STAGE7_O3_ISOTROPY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `EXPECTS_TARGET` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Query validation expected target node. |
| `STARTS_AT` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Query validation start node. |

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
