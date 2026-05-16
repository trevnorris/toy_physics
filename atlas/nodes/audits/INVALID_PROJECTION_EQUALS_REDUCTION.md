---
id: INVALID_PROJECTION_EQUALS_REDUCTION
title: 'Invalid inference: projection equals reduction'
type: invalid_inference
layer: status_audit
status: forbidden_by_v07_firewall
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: This is a deliberately forbidden inference used by the negative query-validation suite.
future_paper_needed: false
status_firewall_ids:
- FIREWALL_PROJECTION_NOT_REDUCTION
incoming_edges:
- source: NEG_QUERY_PROJECTION_EQUALS_REDUCTION
  relation: FORBIDS_TARGET
  status: v07
  note: This query is designed to block this invalid inference.
- source: FIREWALL_PROJECTION_NOT_REDUCTION
  relation: GUARDS_AGAINST
  status: v07
  note: Projection is an exact measurement definition; reduction is a controlled dynamical simplification.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/forbidden_by_v07_firewall
- topic/moving_throat
- topic/projection
- type/invalid_inference
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Invalid inference: projection equals reduction

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `INVALID_PROJECTION_EQUALS_REDUCTION`
> **Status:** `forbidden_by_v07_firewall`
> **Layer:** `status_audit`
> **Type:** `invalid_inference`

## Summary

This is a deliberately forbidden inference used by the negative query-validation suite.

## Physical Meaning

This is a deliberately forbidden inference used by the negative query-validation suite.

## Mathematical Role

- Layer: `status_audit`
- Type: `invalid_inference`
- Status: `forbidden_by_v07_firewall`

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
- [[FIREWALL_PROJECTION_NOT_REDUCTION]]

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FORBIDS_TARGET` | [[NEG_QUERY_PROJECTION_EQUALS_REDUCTION]] | This query is designed to block this invalid inference. |
| `GUARDS_AGAINST` | [[FIREWALL_PROJECTION_NOT_REDUCTION]] | Projection is an exact measurement definition; reduction is a controlled dynamical simplification. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
