---
id: EQ_GNLS_BULK
title: Bulk gauged GNLS equation
type: equation
layer: equation_anchor
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Exact matter equation from the parent action.
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
  line: 1860
  heading_level: paragraph
  heading: Matter current from the action.
  nearest_label:
    name: eq:Jpsi-variational
    line: 1867
  nearby_labels:
  - name: eq:Jpsi-variational
    line: 1867
  - name: eq:Jpsi-spatial
    line: 1875
  - name: eq:Jpsi-time
    line: 1880
  match_basis: semantic_heading_match
  match_score: 0.632
  confidence: medium
math_ids:
- MATH_GNLS_PSI
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
source_ids:
- FILE_4D_PARENT
- FILE_MOVING_THROAT_COMPACT
outgoing_edges:
- target: MATH_GNLS_PSI
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_GNLS_BULK.
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_GNLS_BULK.
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_GNLS_BULK.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact
- topic/moving_throat
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Bulk gauged GNLS equation

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_GNLS_BULK`
> **Status:** `exact`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Exact matter equation from the parent action.

## Physical Meaning

Exact matter equation from the parent action.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact`
- Parent node: [[MATH_GNLS_PSI]]

## Equation

$$
iℏ D_t ψ = [-(ℏ²/2m)D_iD_i + V_conf(X;Σ)+h(ρ)]ψ
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_GNLS_PSI]]

### Related equations
- none

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
| `ANCHORS` | [[MATH_GNLS_PSI]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references EQ_GNLS_BULK. |
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_GNLS_BULK. |
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_GNLS_BULK. |

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
- Line: `1860`
- Heading: Matter current from the action.
- Nearest label: `eq:Jpsi-variational` at line `1867`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
