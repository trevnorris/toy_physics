---
id: EQ_BULK_CONTINUITY
title: Exact bulk continuity
type: equation
layer: equation_anchor
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Exact number continuity in the 4D spatial bulk.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/moving_throat_pde_program_compact.md
- 4d_summary.md
- moving_throat_pde_program_compact.md
legacy_sources:
- 4d_summary.md
- moving_throat_pde_program_compact.md
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_MOVING_THROAT_COMPACT]]'
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 592
  heading_level: subsection
  heading: Conserved current and continuity identity
  nearest_label:
    name: eq:bulk-continuity
    line: 592
  nearby_labels:
  - name: eq:bulk-continuity
    line: 592
  match_basis: semantic_label_match
  match_score: 0.759
  confidence: high
math_ids:
- MATH_PROJECTED_CONTINUITY
equation_ids:
- EQ_PROJECTED_CONTINUITY
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
source_ids:
- FILE_4D_PARENT
- FILE_MOVING_THROAT_COMPACT
outgoing_edges:
- target: MATH_PROJECTED_CONTINUITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_PROJECTED_CONTINUITY
  relation: PROJECTS_TO
  status: exact
  note: Projection with W(w) gives the open brane identity.
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_BULK_CONTINUITY.
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_BULK_CONTINUITY.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact
- topic/moving_throat
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Exact bulk continuity

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_BULK_CONTINUITY`  
> **Status:** `exact`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Exact number continuity in the 4D spatial bulk.

## Physical Meaning

Exact number continuity in the 4D spatial bulk.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact`
- Parent node: [[MATH_PROJECTED_CONTINUITY]]

## Equation

$$
∂_t ρ + ∂_i j^i = 0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_PROJECTED_CONTINUITY]]

### Related equations
- [[EQ_PROJECTED_CONTINUITY]]

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_MOVING_THROAT_COMPACT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_PROJECTED_CONTINUITY]] | Equation anchor belongs to or formalizes this graph node. |
| `PROJECTS_TO` | [[EQ_PROJECTED_CONTINUITY]] | Projection with W(w) gives the open brane identity. |
| `SUPPORTS_CLAIM` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_BULK_CONTINUITY. |
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_BULK_CONTINUITY. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_MOVING_THROAT_COMPACT]]

### Source files
- `research/4d/paper/4d.tex`
- `notes/moving_throat_pde_program_compact.md`
- `4d_summary.md`
- `moving_throat_pde_program_compact.md`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `592`
- Heading: Conserved current and continuity identity
- Nearest label: `eq:bulk-continuity` at line `592`
- Match basis: `semantic_label_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
