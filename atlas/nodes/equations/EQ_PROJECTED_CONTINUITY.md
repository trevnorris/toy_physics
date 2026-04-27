---
id: EQ_PROJECTED_CONTINUITY
title: Projected brane continuity with leakage
type: equation
layer: equation_anchor
status: exact_after_projection_definition
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Projection with W(w) turns conservative bulk continuity into open-system brane continuity.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_plasma/paper/4d_plasma.tex
- notes/pde_audit_full.md
- 4d_summary.md
- 4d_plasma_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- 4d_plasma_summary.md
- pde_audit_full.md
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[FILE_PLASMA]]'
tex_anchor:
  file: research/4d_plasma/paper/4d_plasma.tex
  line: 933
  heading_level: paragraph
  heading: Projected (brane) continuity with leakage.
  nearest_label:
    name: eq:projected_continuity_species_repeat
    line: 939
  nearby_labels:
  - name: eq:projected_continuity_species_repeat
    line: 939
  match_basis: semantic_heading_match
  match_score: 1.0
  confidence: high
math_ids:
- MATH_PROJECTED_CONTINUITY
equation_ids:
- EQ_BULK_CONTINUITY
- EQ_LONGITUDINAL_IDENTITY
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
source_ids:
- FILE_4D_PARENT
- FILE_PDE_AUDIT
- FILE_PLASMA
outgoing_edges:
- target: MATH_PROJECTED_CONTINUITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_LONGITUDINAL_IDENTITY
  relation: FEEDS
  status: exact
  note: Longitudinal Helmholtz split of brane velocity gives identity.
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: SUPPORTS_CLAIM
  status: exact_projection_plus_controlled_hook
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_PROJECTED_CONTINUITY.
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_PROJECTED_CONTINUITY.
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_PROJECTED_CONTINUITY.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_PROJECTED_CONTINUITY.
- source: FILE_PLASMA
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_PROJECTED_CONTINUITY.
- source: EQ_BULK_CONTINUITY
  relation: PROJECTS_TO
  status: exact
  note: Projection with W(w) gives the open brane identity.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_after_projection_definition
- topic/moving_throat
- topic/projection
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Projected brane continuity with leakage

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_PROJECTED_CONTINUITY`  
> **Status:** `exact_after_projection_definition`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Projection with W(w) turns conservative bulk continuity into open-system brane continuity.

## Physical Meaning

Projection with W(w) turns conservative bulk continuity into open-system brane continuity.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_after_projection_definition`
- Parent node: [[MATH_PROJECTED_CONTINUITY]]

## Equation

$$
∂_t ρ_brane + ∇_3·j_brane = S_leak
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_PROJECTED_CONTINUITY]]

### Related equations
- [[EQ_BULK_CONTINUITY]]
- [[EQ_LONGITUDINAL_IDENTITY]]

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_PDE_AUDIT]]
- [[FILE_PLASMA]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_PROJECTED_CONTINUITY]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_LONGITUDINAL_IDENTITY]] | Longitudinal Helmholtz split of brane velocity gives identity. |
| `SUPPORTS_CLAIM` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references EQ_PROJECTED_CONTINUITY. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references EQ_PROJECTED_CONTINUITY. |
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_PROJECTED_CONTINUITY. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_PROJECTED_CONTINUITY. |
| `CONTAINS_EQUATION` | [[FILE_PLASMA]] | Source artifact contains or supports EQ_PROJECTED_CONTINUITY. |
| `PROJECTS_TO` | [[EQ_BULK_CONTINUITY]] | Projection with W(w) gives the open brane identity. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_PDE_AUDIT]]
- [[FILE_PLASMA]]

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_plasma/paper/4d_plasma.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `4d_plasma_summary.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/4d_plasma/paper/4d_plasma.tex`
- Line: `933`
- Heading: Projected (brane) continuity with leakage.
- Nearest label: `eq:projected_continuity_species_repeat` at line `939`
- Match basis: `semantic_heading_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
