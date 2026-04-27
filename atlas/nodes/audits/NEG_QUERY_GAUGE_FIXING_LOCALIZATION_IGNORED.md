---
id: NEG_QUERY_GAUGE_FIXING_LOCALIZATION_IGNORED
title: Can the unlocalized Maxwell gauge-fixing term be ignored without a safe interpretation or patch?
type: negative_query_test
layer: query_validation
status: v07_negative_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Can the unlocalized Maxwell gauge-fixing term be ignored without a safe interpretation or patch?
future_paper_needed: false
math_ids:
- MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
status_firewall_ids:
- FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED
outgoing_edges:
- target: INVALID_GAUGE_FIXING_LOCALIZATION_IGNORED
  relation: FORBIDS_TARGET
  status: v07
  note: This query is designed to block this invalid inference.
- target: FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED
  relation: PROTECTED_BY
  status: v07
  note: Negative query is protected by a status-firewall rule.
- target: MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MATH_MAXWELL_WEIGHTED_GAUGE_FIXING.
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
- topic/maxwell
- topic/moving_throat
- type/negative_query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Can the unlocalized Maxwell gauge-fixing term be ignored without a safe interpretation or patch?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `NEG_QUERY_GAUGE_FIXING_LOCALIZATION_IGNORED`  
> **Status:** `v07_negative_query`  
> **Layer:** `query_validation`  
> **Type:** `negative_query_test`

## Summary

Can the unlocalized Maxwell gauge-fixing term be ignored without a safe interpretation or patch?

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
- [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- [[FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FORBIDS_TARGET` | [[INVALID_GAUGE_FIXING_LOCALIZATION_IGNORED]] | This query is designed to block this invalid inference. |
| `PROTECTED_BY` | [[FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED]] | Negative query is protected by a status-firewall rule. |
| `STARTS_AT` | [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]] | Negative query starts from MATH_MAXWELL_WEIGHTED_GAUGE_FIXING. |

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
