---
id: QUERY_ZERO_MODE_DOWNSTREAM
title: What depends on zero-mode Maxwell?
type: query_test
layer: query_validation
status: v06_seed_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Hydrogen/atomic reduced sector and brane Maxwell/Coulomb branch depend on controlled zero-mode Maxwell, while mixed core remains microscopic.
future_paper_needed: false
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
outgoing_edges:
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v06_seed_query
- topic/atom
- topic/maxwell
- topic/projection
- type/query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# What depends on zero-mode Maxwell?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `QUERY_ZERO_MODE_DOWNSTREAM`
> **Status:** `v06_seed_query`
> **Layer:** `query_validation`
> **Type:** `query_test`

## Summary

Hydrogen/atomic reduced sector and brane Maxwell/Coulomb branch depend on controlled zero-mode Maxwell, while mixed core remains microscopic.

## Physical Meaning

Hydrogen/atomic reduced sector and brane Maxwell/Coulomb branch depend on controlled zero-mode Maxwell, while mixed core remains microscopic.

## Mathematical Role

- Layer: `query_validation`
- Type: `query_test`
- Status: `v06_seed_query`
- Target: `CLAIM_ATOMIC_HYDROGEN_ZERO_MODE`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `EXPECTS_TARGET` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Query validation expected target node. |
| `STARTS_AT` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Query validation start node. |

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
