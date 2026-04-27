---
id: NEG_QUERY_HYDROGEN_FULL_MOVING_THROAT_THEOREM
title: Does the hydrogenic zero-mode reduction prove full atomic structure from the solved moving-throat PDE?
type: negative_query_test
layer: query_validation
status: v07_negative_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Does the hydrogenic zero-mode reduction prove full atomic structure from the solved moving-throat PDE?
future_paper_needed: false
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
status_firewall_ids:
- FIREWALL_ATOM_REDUCED_SECTOR
outgoing_edges:
- target: INVALID_HYDROGEN_FULL_PDE_THEOREM
  relation: FORBIDS_TARGET
  status: v07
  note: This query is designed to block this invalid inference.
- target: FIREWALL_ATOM_REDUCED_SECTOR
  relation: PROTECTED_BY
  status: v07
  note: Negative query is protected by a status-firewall rule.
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: STARTS_AT
  status: v07
  note: Negative query starts from CLAIM_ATOMIC_HYDROGEN_ZERO_MODE.
incoming_edges:
- source: ATLAS_PHASE_7_NEGATIVE_VALIDATION
  relation: CONTAINS_NEGATIVE_QUERY
  status: v07
  note: Negative query belongs to the v0.7 failure-mode suite.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v07_negative_query
- topic/atom
- topic/maxwell
- topic/moving_throat
- topic/projection
- type/negative_query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Does the hydrogenic zero-mode reduction prove full atomic structure from the solved moving-throat PDE?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `NEG_QUERY_HYDROGEN_FULL_MOVING_THROAT_THEOREM`  
> **Status:** `v07_negative_query`  
> **Layer:** `query_validation`  
> **Type:** `negative_query_test`

## Summary

Does the hydrogenic zero-mode reduction prove full atomic structure from the solved moving-throat PDE?

## Physical Meaning

No additional physical interpretation is recorded beyond the graph metadata.

## Mathematical Role

- Layer: `query_validation`
- Type: `negative_query_test`
- Status: `v07_negative_query`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_ATOM_REDUCED_SECTOR]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FORBIDS_TARGET` | [[INVALID_HYDROGEN_FULL_PDE_THEOREM]] | This query is designed to block this invalid inference. |
| `PROTECTED_BY` | [[FIREWALL_ATOM_REDUCED_SECTOR]] | Negative query is protected by a status-firewall rule. |
| `STARTS_AT` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Negative query starts from CLAIM_ATOMIC_HYDROGEN_ZERO_MODE. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_NEGATIVE_QUERY` | [[ATLAS_PHASE_7_NEGATIVE_VALIDATION]] | Negative query belongs to the v0.7 failure-mode suite. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
