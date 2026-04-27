---
id: MT_V2_15_25PN_4PN_INTERFACE
title: V2-15 2.5PN/4PN interface audit
type: interface_gate
layer: status_audit
status: conditional_on_quadrupole_normalization
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Shows the same quadrupole normalization controls both the 2.5PN Burke-Thorne branch and the 4PN tail coefficient.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
outgoing_edges:
- target: PN_2_5_QUAD_NARROWING
  relation: INTERFACES
  status: conditional
  note: Relates moving-throat quadrupole branch to 2.5PN channel.
- target: PN_4_LOCAL_TAIL
  relation: INTERFACES
  status: conditional
  note: Relates same quadrupole branch to 4PN tail coefficient.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/conditional_on_quadrupole_normalization
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/interface_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-15 2.5PN/4PN interface audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_15_25PN_4PN_INTERFACE`  
> **Status:** `conditional_on_quadrupole_normalization`  
> **Layer:** `status_audit`  
> **Type:** `interface_gate`

## Summary

Shows the same quadrupole normalization controls both the 2.5PN Burke-Thorne branch and the 4PN tail coefficient.

## Physical Meaning

Shows the same quadrupole normalization controls both the 2.5PN Burke-Thorne branch and the 4PN tail coefficient.

## Mathematical Role

- Layer: `status_audit`
- Type: `interface_gate`
- Status: `conditional_on_quadrupole_normalization`

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
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `INTERFACES` | [[PN_2_5_QUAD_NARROWING]] | Relates moving-throat quadrupole branch to 2.5PN channel. |
| `INTERFACES` | [[PN_4_LOCAL_TAIL]] | Relates same quadrupole branch to 4PN tail coefficient. |

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
