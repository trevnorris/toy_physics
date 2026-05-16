---
id: MT_V2_23_OPEN_BRANCH_SOLVER
title: V2-23 minimal open-throat branch solver
type: simulation_gate
layer: status_audit
status: target_miss
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: First target-blind branch families did not land on target packet.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-23
- pde_audit_full.md:V2-26
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: DOES_NOT_CLOSE
  status: open
  note: First target-blind simulated branch missed target; exporter still incomplete.
- target: MT_V2_26_STATUS
  relation: FEEDS_STATUS
  status: audit
  note: Target-blind misses inform honest status boundary.
incoming_edges:
- source: MT_V2_22C_SMOKE_PIPELINE
  relation: PRECEDES
  status: implemented
  note: Smoke pipeline precedes first real residual extraction.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/target_miss
- topic/moving_throat
- type/simulation_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-23 minimal open-throat branch solver

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_23_OPEN_BRANCH_SOLVER`
> **Status:** `target_miss`
> **Layer:** `status_audit`
> **Type:** `simulation_gate`

## Summary

First target-blind branch families did not land on target packet.

## Physical Meaning

First target-blind branch families did not land on target packet.

## Mathematical Role

- Layer: `status_audit`
- Type: `simulation_gate`
- Status: `target_miss`

## Equation

$$
0 target-passing candidates in current reduced/manufactured families
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

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
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DOES_NOT_CLOSE` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | First target-blind simulated branch missed target; exporter still incomplete. |
| `FEEDS_STATUS` | [[MT_V2_26_STATUS]] | Target-blind misses inform honest status boundary. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `PRECEDES` | [[MT_V2_22C_SMOKE_PIPELINE]] | Smoke pipeline precedes first real residual extraction. |

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
