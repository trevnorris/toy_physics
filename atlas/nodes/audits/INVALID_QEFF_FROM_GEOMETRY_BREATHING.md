---
id: INVALID_QEFF_FROM_GEOMETRY_BREATHING
title: 'Invalid inference: q_eff is a throat-breathing/geometric variable'
type: invalid_inference
layer: status_audit
status: forbidden_by_v07_firewall
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: This is a deliberately forbidden inference used by the negative query-validation suite.
future_paper_needed: false
status_firewall_ids:
- FIREWALL_QEFF_THICKNESS_NOT_BREATHING
incoming_edges:
- source: NEG_QUERY_QEFF_FROM_THROAT_GEOMETRY
  relation: FORBIDS_TARGET
  status: v07
  note: This query is designed to block this invalid inference.
- source: FIREWALL_QEFF_THICKNESS_NOT_BREATHING
  relation: GUARDS_AGAINST
  status: v07
  note: Keep microscopic branch charge q_star and observed brane coupling q_eff separate from throat radius/length dynamics.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/forbidden_by_v07_firewall
- topic/charge
- topic/moving_throat
- type/invalid_inference
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Invalid inference: q_eff is a throat-breathing/geometric variable

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `INVALID_QEFF_FROM_GEOMETRY_BREATHING`  
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
- [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]]

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FORBIDS_TARGET` | [[NEG_QUERY_QEFF_FROM_THROAT_GEOMETRY]] | This query is designed to block this invalid inference. |
| `GUARDS_AGAINST` | [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]] | Keep microscopic branch charge q_star and observed brane coupling q_eff separate from throat radius/length ... |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
