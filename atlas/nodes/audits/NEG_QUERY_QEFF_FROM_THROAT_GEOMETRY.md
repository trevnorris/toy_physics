---
id: NEG_QUERY_QEFF_FROM_THROAT_GEOMETRY
title: Can q_eff be inferred from throat radius, length, or breathing variables?
type: negative_query_test
layer: query_validation
status: v07_negative_query
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Can q_eff be inferred from throat radius, length, or breathing variables?
future_paper_needed: false
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
status_firewall_ids:
- FIREWALL_QEFF_THICKNESS_NOT_BREATHING
outgoing_edges:
- target: INVALID_QEFF_FROM_GEOMETRY_BREATHING
  relation: FORBIDS_TARGET
  status: v07
  note: This query is designed to block this invalid inference.
- target: FIREWALL_QEFF_THICKNESS_NOT_BREATHING
  relation: PROTECTED_BY
  status: v07
  note: Negative query is protected by a status-firewall rule.
- target: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: STARTS_AT
  status: v07
  note: Negative query starts from PHYS_BRANE_BULK_THROAT_DEFECT.
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
- topic/charge
- topic/moving_throat
- type/negative_query_test
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Can q_eff be inferred from throat radius, length, or breathing variables?

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `NEG_QUERY_QEFF_FROM_THROAT_GEOMETRY`
> **Status:** `v07_negative_query`
> **Layer:** `query_validation`
> **Type:** `negative_query_test`

## Summary

Can q_eff be inferred from throat radius, length, or breathing variables?

## Physical Meaning

No additional physical interpretation is recorded beyond the graph metadata.

## Mathematical Role

- Layer: `query_validation`
- Type: `negative_query_test`
- Status: `v07_negative_query`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FORBIDS_TARGET` | [[INVALID_QEFF_FROM_GEOMETRY_BREATHING]] | This query is designed to block this invalid inference. |
| `PROTECTED_BY` | [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]] | Negative query is protected by a status-firewall rule. |
| `STARTS_AT` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Negative query starts from PHYS_BRANE_BULK_THROAT_DEFECT. |

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
