---
id: EQ_WALL_P2_STIFFNESS
title: Grouped P2 stiffness gate
type: equation
layer: equation_anchor
status: effective_linear_wall_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Shows grouped real P2 is next harmonic family of the wall action.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- research/pde_ledger/paper/stages/stage_002.tex
- pde_audit_full.md
- moving_throat_pde_stage002_breathing_reduction.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_stage002_breathing_reduction.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
equation_ids:
- EQ_BDG_SCHUR_KERNEL
- EQ_WALL_MODAL_PDE
claim_ids:
- CLAIM_STAGE2_AL_RECOVERY
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_BDG_SCHUR_KERNEL
  relation: FEEDS
  status: reduced
  note: Bare P2 wall operator is then dressed by BdG support.
- target: CLAIM_STAGE2_AL_RECOVERY
  relation: SUPPORTS_CLAIM
  status: exact_within_wall_action
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WALL_P2_STIFFNESS.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WALL_P2_STIFFNESS.
- source: EQ_WALL_MODAL_PDE
  relation: SPECIALIZES_TO
  status: effective_closure
  note: l=2 specialization gives grouped P2 stiffness.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/effective_linear_wall_closure
- topic/moving_throat
- topic/projection
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Grouped P2 stiffness gate

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_WALL_P2_STIFFNESS`  
> **Status:** `effective_linear_wall_closure`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Shows grouped real P2 is next harmonic family of the wall action.

## Physical Meaning

Shows grouped real P2 is next harmonic family of the wall action.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `effective_linear_wall_closure`
- Parent node: `PHYS_GROUPED_P2_SHAPE_RESPONSE`

## Equation

$$
K_2 = ∫dw [T_w β_2′² + (K_η+6T_Ω)β_2²]
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]

### Related math nodes
- none

### Related equations
- [[EQ_BDG_SCHUR_KERNEL]]
- [[EQ_WALL_MODAL_PDE]]

### Related claims
- [[CLAIM_STAGE2_AL_RECOVERY]]

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
| `ANCHORS` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_BDG_SCHUR_KERNEL]] | Bare P2 wall operator is then dressed by BdG support. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE2_AL_RECOVERY]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_WALL_P2_STIFFNESS. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_WALL_P2_STIFFNESS. |
| `SPECIALIZES_TO` | [[EQ_WALL_MODAL_PDE]] | l=2 specialization gives grouped P2 stiffness. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`
- `research/pde_ledger/paper/stages/stage_002.tex`
- `pde_audit_full.md`
- `moving_throat_pde_stage002_breathing_reduction.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
