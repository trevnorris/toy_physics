---
id: MT_V2_16_NO_REFIT
title: V2-16 branch-freeze/no-refit protocol
type: protocol_gate
layer: status_audit
status: mandatory
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Freeze action status, gauge convention, branch class, support basis, ports, stability gates, extraction formulas before target comparison.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-16
physical_ids:
- PHYS_TARGET_BLIND_BRANCH_SELECTION
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: GOVERNS
  status: mandatory
  note: Branch exporter must freeze protocol before target comparison.
incoming_edges:
- source: PHYS_TARGET_BLIND_BRANCH_SELECTION
  relation: IMPLEMENTS
  status: mandatory
  note: No-refit protocol enforces target-blind branch freezing.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/mandatory
- topic/moving_throat
- type/protocol_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-16 branch-freeze/no-refit protocol

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_16_NO_REFIT`  
> **Status:** `mandatory`  
> **Layer:** `status_audit`  
> **Type:** `protocol_gate`

## Summary

Freeze action status, gauge convention, branch class, support basis, ports, stability gates, extraction formulas before target comparison.

## Physical Meaning

Freeze action status, gauge convention, branch class, support basis, ports, stability gates, extraction formulas before target comparison.

## Mathematical Role

- Layer: `status_audit`
- Type: `protocol_gate`
- Status: `mandatory`

## Equation

$$
no post-hoc retuning
$$

$$
do not project realized branch to target
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_TARGET_BLIND_BRANCH_SELECTION]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- none

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GOVERNS` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Branch exporter must freeze protocol before target comparison. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `IMPLEMENTS` | [[PHYS_TARGET_BLIND_BRANCH_SELECTION]] | No-refit protocol enforces target-blind branch freezing. |

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
