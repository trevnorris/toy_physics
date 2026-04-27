---
id: OPEN_NONLINEAR_S_SIGMA
title: Nonlinear throat action S_Sigma
type: parent_action_gate
layer: open_gate
status: open_patch_required
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Promote or derive nonlinear throat action whose quadratic limit is S_eta^(2).
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
open_gate_ids:
- OPEN_PARENT_PROMOTION_S_SIGMA
status_firewall_ids:
- FIREWALL_WALL_COEFFS_BRANCH_DATA
incoming_edges:
- source: BACKLINK_PDE_AUDIT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_NONLINEAR_S_SIGMA.
- source: FIREWALL_WALL_COEFFS_BRANCH_DATA
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: REFINES_TO
  status: open
  note: Best patch is a nonlinear S_Sigma whose quadratic limit is S_eta^(2).
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open_patch_required
- topic/moving_throat
- type/parent_action_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Nonlinear throat action S_Sigma

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_NONLINEAR_S_SIGMA`  
> **Status:** `open_patch_required`  
> **Layer:** `open_gate`  
> **Type:** `parent_action_gate`

## Summary

Promote or derive nonlinear throat action whose quadratic limit is S_eta^(2).

## What Remains Open

Promote or derive nonlinear throat action whose quadratic limit is S_eta^(2).

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Promote or derive nonlinear throat action whose quadratic limit is S_eta^(2).

## Mathematical Role

- Layer: `open_gate`
- Type: `parent_action_gate`
- Status: `open_patch_required`

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
- [[OPEN_PARENT_PROMOTION_S_SIGMA]]

### Status firewalls
- [[FIREWALL_WALL_COEFFS_BRANCH_DATA]]

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FLAGS_OPEN_GATE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block flags open gate OPEN_NONLINEAR_S_SIGMA. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_WALL_COEFFS_BRANCH_DATA]] | Firewall preserves this correct status boundary. |
| `REFINES_TO` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Best patch is a nonlinear S_Sigma whose quadratic limit is S_eta^(2). |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_pde_program_compact.md`
- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
