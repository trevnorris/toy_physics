---
id: QUERY_G2_XI1
title: How does the anomaly branch connect to the moving-throat quotient?
type: query_test
layer: query_validation
status: v06_seed_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The anomaly/g-2 branch uses q_tr,q_nt,q_eta and the Xi1=P1/P0 prefactor slope on the weak-axisymmetric branch.
future_paper_needed: false
math_ids:
- MATH_XI1_PREF_SLOPE
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
outgoing_edges:
- target: MATH_XI1_PREF_SLOPE
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/audits
- atlas/node
- layer/query_validation
- status/v06_seed_query
- topic/g2
- topic/quadrupole
- type/query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# How does the anomaly branch connect to the moving-throat quotient?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `QUERY_G2_XI1`
> **Status:** `v06_seed_query`
> **Layer:** `query_validation`
> **Type:** `query_test`

## Summary

The anomaly/g-2 branch uses q_tr,q_nt,q_eta and the Xi1=P1/P0 prefactor slope on the weak-axisymmetric branch.

## Physical Meaning

The anomaly/g-2 branch uses q_tr,q_nt,q_eta and the Xi1=P1/P0 prefactor slope on the weak-axisymmetric branch.

## Mathematical Role

- Layer: `query_validation`
- Type: `query_test`
- Status: `v06_seed_query`
- Target: `MATH_XI1_PREF_SLOPE`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_XI1_PREF_SLOPE]]

### Related equations
- none

### Related claims
- [[CLAIM_G2_COMMON_QUOTIENT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `EXPECTS_TARGET` | [[MATH_XI1_PREF_SLOPE]] | Query validation expected target node. |
| `STARTS_AT` | [[CLAIM_G2_COMMON_QUOTIENT]] | Query validation start node. |

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
