---
id: NEG_QUERY_25PN_UNCONDITIONAL
title: Is the 2.5PN package an unconditional solved moving-throat PDE theorem?
type: negative_query_test
layer: query_validation
status: v07_negative_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Is the 2.5PN package an unconditional solved moving-throat PDE theorem?
future_paper_needed: false
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
status_firewall_ids:
- FIREWALL_25PN_CONDITIONAL
outgoing_edges:
- target: INVALID_25PN_UNCONDITIONAL_THEOREM
  relation: FORBIDS_TARGET
  status: v07
  note: This query is designed to block this invalid inference.
- target: FIREWALL_25PN_CONDITIONAL
  relation: PROTECTED_BY
  status: v07
  note: Negative query is protected by a status-firewall rule.
- target: CLAIM_25PN_QUAD_NARROWING
  relation: STARTS_AT
  status: v07
  note: Negative query starts from CLAIM_25PN_QUAD_NARROWING.
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
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/negative_query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Is the 2.5PN package an unconditional solved moving-throat PDE theorem?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `NEG_QUERY_25PN_UNCONDITIONAL`  
> **Status:** `v07_negative_query`  
> **Layer:** `query_validation`  
> **Type:** `negative_query_test`

## Summary

Is the 2.5PN package an unconditional solved moving-throat PDE theorem?

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
- [[CLAIM_25PN_QUAD_NARROWING]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_25PN_CONDITIONAL]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FORBIDS_TARGET` | [[INVALID_25PN_UNCONDITIONAL_THEOREM]] | This query is designed to block this invalid inference. |
| `PROTECTED_BY` | [[FIREWALL_25PN_CONDITIONAL]] | Negative query is protected by a status-firewall rule. |
| `STARTS_AT` | [[CLAIM_25PN_QUAD_NARROWING]] | Negative query starts from CLAIM_25PN_QUAD_NARROWING. |

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
