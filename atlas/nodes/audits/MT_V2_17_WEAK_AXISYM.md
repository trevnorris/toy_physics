---
id: MT_V2_17_WEAK_AXISYM
title: V2-17 weak-axisymmetric P2 splitting audit
type: audit_gate
layer: status_audit
status: exact_first_order
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Freezes weak-axisymmetric grouped signature, hidden-even relation, outgoing-prefactor transport, and Xi1 packet.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
math_ids:
- MATH_WEAK_AXISYM_SPLITTING
outgoing_edges:
- target: MATH_WEAK_AXISYM_SPLITTING
  relation: FREEZES
  status: exact_first_order
  note: Weak-axisymmetric first-order grouped signature.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/exact_first_order
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-17 weak-axisymmetric P2 splitting audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_17_WEAK_AXISYM`  
> **Status:** `exact_first_order`  
> **Layer:** `status_audit`  
> **Type:** `audit_gate`

## Summary

Freezes weak-axisymmetric grouped signature, hidden-even relation, outgoing-prefactor transport, and Xi1 packet.

## Physical Meaning

Freezes weak-axisymmetric grouped signature, hidden-even relation, outgoing-prefactor transport, and Xi1 packet.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `exact_first_order`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WEAK_AXISYM_SPLITTING]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FREEZES` | [[MATH_WEAK_AXISYM_SPLITTING]] | Weak-axisymmetric first-order grouped signature. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
