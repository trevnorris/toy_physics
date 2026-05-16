---
id: QUERY_MIXED_SECTOR
title: Why can the outgoing bridge use the mixed sector?
type: query_test
layer: query_validation
status: v06_seed_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Mixed fields are microscopic gauge-invariant channels; Stage 4 uses localized Maxwell/mixed coordinates to carry the outgoing quadrupole bridge.
future_paper_needed: false
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
outgoing_edges:
- target: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v06_seed_query
- topic/maxwell
- topic/quadrupole
- type/query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Why can the outgoing bridge use the mixed sector?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `QUERY_MIXED_SECTOR`
> **Status:** `v06_seed_query`
> **Layer:** `query_validation`
> **Type:** `query_test`

## Summary

Mixed fields are microscopic gauge-invariant channels; Stage 4 uses localized Maxwell/mixed coordinates to carry the outgoing quadrupole bridge.

## Physical Meaning

Mixed fields are microscopic gauge-invariant channels; Stage 4 uses localized Maxwell/mixed coordinates to carry the outgoing quadrupole bridge.

## Mathematical Role

- Layer: `query_validation`
- Type: `query_test`
- Status: `v06_seed_query`
- Target: `CLAIM_PROJECTED_EM_OUTGOING_BRIDGE`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `EXPECTS_TARGET` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Query validation expected target node. |
| `STARTS_AT` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Query validation start node. |

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
