---
id: QUERY_4PN_TAIL_BLOCKER
title: What blocks full 4PN tail completion?
type: query_test
layer: query_validation
status: v06_seed_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: 4PN local is closed within hierarchy; full tail inherits the 2.5PN passive/outgoing STF quadrupole normalization gate.
future_paper_needed: false
claim_ids:
- CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
outgoing_edges:
- target: OPEN_QUAD_NORMALIZATION
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- target: CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v06_seed_query
- topic/pn_chain
- topic/quadrupole
- type/query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# What blocks full 4PN tail completion?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `QUERY_4PN_TAIL_BLOCKER`  
> **Status:** `v06_seed_query`  
> **Layer:** `query_validation`  
> **Type:** `query_test`

## Summary

4PN local is closed within hierarchy; full tail inherits the 2.5PN passive/outgoing STF quadrupole normalization gate.

## Physical Meaning

4PN local is closed within hierarchy; full tail inherits the 2.5PN passive/outgoing STF quadrupole normalization gate.

## Mathematical Role

- Layer: `query_validation`
- Type: `query_test`
- Status: `v06_seed_query`
- Target: `OPEN_QUAD_NORMALIZATION`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `EXPECTS_TARGET` | [[OPEN_QUAD_NORMALIZATION]] | Query validation expected target node. |
| `STARTS_AT` | [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]] | Query validation start node. |

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
