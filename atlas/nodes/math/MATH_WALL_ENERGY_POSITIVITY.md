---
id: MATH_WALL_ENERGY_POSITIVITY
title: Wall energy positivity gate
type: stability_gate
layer: math_object
status: exact_if_S_eta_included
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Local quadratic positivity conditions for the effective wall closure.
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
  relation: IMPLIES_GATE
  status: exact_if_closure
  note: Wall PDE has local positivity gates.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_if_s_eta_included
- topic/moving_throat
- type/stability_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Wall energy positivity gate

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_WALL_ENERGY_POSITIVITY`  
> **Status:** `exact_if_S_eta_included`  
> **Layer:** `math_object`  
> **Type:** `stability_gate`

## Summary

Local quadratic positivity conditions for the effective wall closure.

## Physical Meaning

Local quadratic positivity conditions for the effective wall closure.

## Mathematical Role

- Layer: `math_object`
- Type: `stability_gate`
- Status: `exact_if_S_eta_included`

## Equation

$$
mu_eta>0
$$

$$
T_w>0
$$

$$
K_eta+l(l+1)T_Omega>=0
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
| `IMPLIES_GATE` | [[MATH_WALL_MODAL_PDE]] | Wall PDE has local positivity gates. |

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
