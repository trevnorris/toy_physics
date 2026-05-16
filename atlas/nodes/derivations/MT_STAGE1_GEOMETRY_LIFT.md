---
id: MT_STAGE1_GEOMETRY_LIFT
title: Stage 1 geometry lift
type: moving_throat_stage
layer: derivation
status: effective_closure_scaffold
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Promotes a,L closure to distributed shape field and l=0/l=2 modal skeleton.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- research/pde_ledger/paper/stages/stage_001.tex
- moving_throat_output_full.md
- moving_throat_pde_stage001_geometry_lift.md
legacy_sources:
- moving_throat_output_full.md
- moving_throat_pde_stage001_geometry_lift.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
tex_anchor:
  file: research/pde_ledger/paper/stages/stage_001.tex
  line: 1
  heading_level: section
  heading: 'Stage 001: Geometry Lift and Linearized PDE Skeleton'
  nearest_label:
    name: stage:001
    line: 2
  nearby_labels:
  - name: stage:001
    line: 2
  match_basis: semantic_heading_match
  match_score: 0.653
  confidence: medium
math_ids:
- MATH_SIGMA_R_FIELD
- MATH_WALL_ACTION_S_ETA
claim_ids:
- CLAIM_STAGE1_GEOMETRY_LIFT
source_ids:
- FILE_MOVING_THROAT_COMPACT
outgoing_edges:
- target: MATH_WALL_ACTION_S_ETA
  relation: INTRODUCES
  status: effective
  note: Stage 1 introduces minimal quadratic wall action as new ansatz.
incoming_edges:
- source: FILE_MOVING_THROAT_COMPACT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: FEEDS_OR_STATUS_OF
  status: effective_geometry_lift
  note: Claim feeds this downstream object, output, or open gate.
- source: MATH_SIGMA_R_FIELD
  relation: ORGANIZES
  status: effective
  note: Stage 1 uses Sigma/R to create distributed geometry field.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/effective_closure_scaffold
- topic/moving_throat
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 1 geometry lift

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE1_GEOMETRY_LIFT`
> **Status:** `effective_closure_scaffold`
> **Layer:** `derivation`
> **Type:** `moving_throat_stage`

## Summary

Promotes a,L closure to distributed shape field and l=0/l=2 modal skeleton.

## Physical Meaning

Promotes a,L closure to distributed shape field and l=0/l=2 modal skeleton.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `effective_closure_scaffold`

## Equation

$$
Sigma=r-R(Omega,w,t)
$$

$$
eta expansion in real harmonics
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_SIGMA_R_FIELD]]
- [[MATH_WALL_ACTION_S_ETA]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `INTRODUCES` | [[MATH_WALL_ACTION_S_ETA]] | Stage 1 introduces minimal quadratic wall action as new ansatz. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_MOVING_THROAT_COMPACT]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Claim feeds this downstream object, output, or open gate. |
| `ORGANIZES` | [[MATH_SIGMA_R_FIELD]] | Stage 1 uses Sigma/R to create distributed geometry field. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]

### Source files
- `notes/moving_throat_notes_full.md`
- `research/pde_ledger/paper/stages/stage_001.tex`
- `moving_throat_output_full.md`
- `moving_throat_pde_stage001_geometry_lift.md`

### TeX anchor
- File: `research/pde_ledger/paper/stages/stage_001.tex`
- Line: `1`
- Heading: Stage 001: Geometry Lift and Linearized PDE Skeleton
- Nearest label: `stage:001` at line `2`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
