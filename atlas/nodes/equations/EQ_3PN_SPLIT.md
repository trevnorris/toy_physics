---
id: EQ_3PN_SPLIT
title: 3PN conservative split
type: equation
layer: equation_anchor
status: full_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: 3PN fixed ADM-chart split into pure-kinetic image, grouped P2 middle block, and geometry completion.
future_paper_needed: false
source_files:
- research/4d_3pn/paper/4d_3pn.tex
- 4d_3pn_summary.md
- 4d_3pn.tex
legacy_sources:
- 4d_3pn_summary.md
- 4d_3pn.tex
source_links:
- '[[FILE_3PN_FULL]]'
tex_anchor:
  file: research/4d_3pn/paper/4d_3pn.tex
  line: 2962
  heading_level: section
  heading: Final conservative 3PN theorem
  nearest_label:
    name: sec:final-theorem
    line: 2963
  nearby_labels:
  - name: sec:final-theorem
    line: 2963
  match_basis: semantic_heading_match
  match_score: 0.659
  confidence: medium
claim_ids:
- CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
source_ids:
- FILE_3PN_FULL
outgoing_edges:
- target: PN_3_GROUPED_P2_FULL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
  relation: SUPPORTS_CLAIM
  status: exact_assembly_within_closure
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_3PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_3PN_SPLIT.
- source: FILE_3PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_3PN_SPLIT.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/full_assembly_within_closure
- topic/pn_chain
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 3PN conservative split

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_3PN_SPLIT`  
> **Status:** `full_assembly_within_closure`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

3PN fixed ADM-chart split into pure-kinetic image, grouped P2 middle block, and geometry completion.

## Physical Meaning

3PN fixed ADM-chart split into pure-kinetic image, grouped P2 middle block, and geometry completion.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `full_assembly_within_closure`
- Parent node: `PN_3_GROUPED_P2_FULL`

## Equation

$$
ΔL_3^GR = Δl_1 v⁸ + L_P2^mid + Δl_15^(g) U⁴
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_3PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[PN_3_GROUPED_P2_FULL]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_3PN_FULL]] | Paper backlink block references EQ_3PN_SPLIT. |
| `CONTAINS_EQUATION` | [[FILE_3PN_FULL]] | Source artifact contains or supports EQ_3PN_SPLIT. |

## Source Anchors

### Source anchor notes
- [[FILE_3PN_FULL]]

### Source files
- `research/4d_3pn/paper/4d_3pn.tex`
- `4d_3pn_summary.md`
- `4d_3pn.tex`

### TeX anchor
- File: `research/4d_3pn/paper/4d_3pn.tex`
- Line: `2962`
- Heading: Final conservative 3PN theorem
- Nearest label: `sec:final-theorem` at line `2963`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
