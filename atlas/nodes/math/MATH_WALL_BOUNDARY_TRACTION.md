---
id: MATH_WALL_BOUNDARY_TRACTION
title: Wall boundary traction
type: boundary_term
layer: math_object
status: exact_if_S_eta_included
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Variational boundary momentum for wall modes; driven mouth/worldtube ports must use the correct traction sign.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
math_ids:
- MATH_WALL_MODAL_PDE
incoming_edges:
- source: MATH_WALL_MODAL_PDE
  relation: HAS_BOUNDARY_TERM
  status: exact_if_closure
  note: Boundary terms define driven/free/Dirichlet port conditions.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_if_s_eta_included
- topic/moving_throat
- type/boundary_term
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Wall boundary traction

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_WALL_BOUNDARY_TRACTION`
> **Status:** `exact_if_S_eta_included`
> **Layer:** `math_object`
> **Type:** `boundary_term`

## Summary

Variational boundary momentum for wall modes; driven mouth/worldtube ports must use the correct traction sign.

## Physical Meaning

Variational boundary momentum for wall modes; driven mouth/worldtube ports must use the correct traction sign.

## Mathematical Role

- Layer: `math_object`
- Type: `boundary_term`
- Status: `exact_if_S_eta_included`

## Equation

$$
p_w=-T_w q_w
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WALL_MODAL_PDE]]

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

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `HAS_BOUNDARY_TERM` | [[MATH_WALL_MODAL_PDE]] | Boundary terms define driven/free/Dirichlet port conditions. |

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
