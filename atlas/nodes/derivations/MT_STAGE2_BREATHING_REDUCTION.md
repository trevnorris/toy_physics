---
id: MT_STAGE2_BREATHING_REDUCTION
title: Stage 2 breathing reduction
type: moving_throat_stage
layer: derivation
status: exact_within_wall_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Shows old a,L equations are lowest axisymmetric truncation of wall field.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- moving_throat_output_full.md
legacy_sources:
- moving_throat_output_full.md
math_ids:
- MATH_WALL_ACTION_S_ETA
- MATH_WALL_MODAL_PDE
claim_ids:
- CLAIM_STAGE2_AL_RECOVERY
outgoing_edges:
- target: MT_STAGE3_BDG_COUPLING
  relation: FEEDS
  status: reduced
  note: Wall modes then coupled to stable matter support modes.
incoming_edges:
- source: CLAIM_STAGE2_AL_RECOVERY
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_wall_action
  note: Claim feeds this downstream object, output, or open gate.
- source: MATH_WALL_ACTION_S_ETA
  relation: REDUCES_TO
  status: exact within closure
  note: Axisymmetric two-mode truncation recovers old a,L closure.
- source: MATH_WALL_MODAL_PDE
  relation: REDUCES_TO
  status: exact_if_closure
  note: Axisymmetric two-profile reduction recovers old a,L closure.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/exact_within_wall_closure
- topic/moving_throat
- topic/projection
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 2 breathing reduction

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE2_BREATHING_REDUCTION`  
> **Status:** `exact_within_wall_closure`  
> **Layer:** `derivation`  
> **Type:** `moving_throat_stage`

## Summary

Shows old a,L equations are lowest axisymmetric truncation of wall field.

## Physical Meaning

Shows old a,L equations are lowest axisymmetric truncation of wall field.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `exact_within_wall_closure`

## Equation

$$
M_AB ddot Q^B+K_AB Q^B=0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WALL_ACTION_S_ETA]]
- [[MATH_WALL_MODAL_PDE]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE2_AL_RECOVERY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MT_STAGE3_BDG_COUPLING]] | Wall modes then coupled to stable matter support modes. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE2_AL_RECOVERY]] | Claim feeds this downstream object, output, or open gate. |
| `REDUCES_TO` | [[MATH_WALL_ACTION_S_ETA]] | Axisymmetric two-mode truncation recovers old a,L closure. |
| `REDUCES_TO` | [[MATH_WALL_MODAL_PDE]] | Axisymmetric two-profile reduction recovers old a,L closure. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/moving_throat_notes_full.md`
- `moving_throat_output_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
