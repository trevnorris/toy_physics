---
id: MATH_PROJECTED_CONTINUITY
title: Projected continuity with leakage
type: identity
layer: math_object
status: exact_after_projection_definition
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Bulk conservation becomes open-system brane balance under projection.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_plasma/paper/4d_plasma.tex
- research/4d_1pn_bridge/paper/4d_1pn_bridge.tex
- 4d_summary.md
- 4d_plasma_summary.md
- 4d_1pn_bridge_summary.md
legacy_sources:
- 4d_summary.md
- 4d_plasma_summary.md
- 4d_1pn_bridge_summary.md
source_links:
- '[[FILE_4D_PARENT]]'
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
- MATH_LONGITUDINAL_IDENTITY
- MATH_W_PROJECTION
equation_ids:
- EQ_BULK_CONTINUITY
- EQ_PROJECTED_CONTINUITY
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
source_ids:
- FILE_4D_PARENT
outgoing_edges:
- target: MATH_LONGITUDINAL_IDENTITY
  relation: DERIVES
  status: exact
  note: Helmholtz split yields exact longitudinal identity.
incoming_edges:
- source: EQ_BULK_CONTINUITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: EQ_PROJECTED_CONTINUITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: MATH_W_PROJECTION
  relation: DERIVES
  status: exact
  note: Projection of exact bulk continuity gives leakage source.
- source: FILE_4D_PARENT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_after_projection_definition
- topic/pn_chain
- topic/projection
- type/identity
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Projected continuity with leakage

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_PROJECTED_CONTINUITY`
> **Status:** `exact_after_projection_definition`
> **Layer:** `math_object`
> **Type:** `identity`

## Summary

Bulk conservation becomes open-system brane balance under projection.

## Physical Meaning

Bulk conservation becomes open-system brane balance under projection.

## Mathematical Role

- Layer: `math_object`
- Type: `identity`
- Status: `exact_after_projection_definition`

## Equation

$$
∂t ρ_brane+∇3·j_brane=S_leak
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_LONGITUDINAL_IDENTITY]]
- [[MATH_W_PROJECTION]]

### Related equations
- [[EQ_BULK_CONTINUITY]]
- [[EQ_PROJECTED_CONTINUITY]]

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DERIVES` | [[MATH_LONGITUDINAL_IDENTITY]] | Helmholtz split yields exact longitudinal identity. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_BULK_CONTINUITY]] | Equation anchor belongs to or formalizes this graph node. |
| `ANCHORS` | [[EQ_PROJECTED_CONTINUITY]] | Equation anchor belongs to or formalizes this graph node. |
| `DERIVES` | [[MATH_W_PROJECTION]] | Projection of exact bulk continuity gives leakage source. |
| `DOCUMENTS` | [[FILE_4D_PARENT]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_plasma/paper/4d_plasma.tex`
- `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex`
- `4d_summary.md`
- `4d_plasma_summary.md`
- `4d_1pn_bridge_summary.md`

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
