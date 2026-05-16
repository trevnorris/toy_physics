---
id: EQ_WALL_MODAL_PDE
title: Wall modal PDE from S_eta
type: equation
layer: equation_anchor
status: effective_linear_wall_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Autonomous wall/modal PDE exists only when S_eta/S_Sigma is included.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- research/pde_ledger/paper/stages/stage_001.tex
- pde_audit_full.md
- moving_throat_pde_stage001_geometry_lift.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_stage001_geometry_lift.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
math_ids:
- MATH_WALL_MODAL_PDE
equation_ids:
- EQ_WALL_P2_STIFFNESS
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
- CLAIM_STAGE1_GEOMETRY_LIFT
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_WALL_MODAL_PDE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_WALL_P2_STIFFNESS
  relation: SPECIALIZES_TO
  status: effective_closure
  note: l=2 specialization gives grouped P2 stiffness.
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: SUPPORTS_CLAIM
  status: strict_parent_fail_effective_wall_pass
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: SUPPORTS_CLAIM
  status: effective_geometry_lift
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WALL_MODAL_PDE.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WALL_MODAL_PDE.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/effective_linear_wall_closure
- topic/moving_throat
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Wall modal PDE from S_eta

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_WALL_MODAL_PDE`
> **Status:** `effective_linear_wall_closure`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Autonomous wall/modal PDE exists only when S_eta/S_Sigma is included.

## Physical Meaning

Autonomous wall/modal PDE exists only when S_eta/S_Sigma is included.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `effective_linear_wall_closure`
- Parent node: [[MATH_WALL_MODAL_PDE]]

## Equation

$$
μ_η q_lm,tt - ∂_w(T_w q_lm,w) + [K_η+l(l+1)T_Ω]q_lm = S_lm
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WALL_MODAL_PDE]]

### Related equations
- [[EQ_WALL_P2_STIFFNESS]]

### Related claims
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_WALL_MODAL_PDE]] | Equation anchor belongs to or formalizes this graph node. |
| `SPECIALIZES_TO` | [[EQ_WALL_P2_STIFFNESS]] | l=2 specialization gives grouped P2 stiffness. |
| `SUPPORTS_CLAIM` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_WALL_MODAL_PDE. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_WALL_MODAL_PDE. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`
- `research/pde_ledger/paper/stages/stage_001.tex`
- `pde_audit_full.md`
- `moving_throat_pde_stage001_geometry_lift.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
