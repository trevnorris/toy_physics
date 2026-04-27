---
id: EQ_2PN_ADM_EQUALITY
title: 2PN ADM equality
type: equation
layer: equation_anchor
status: full_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Conservative 2PN Lagrangian Legendre-transforms exactly to generic-frame ADM target inside declared hierarchy.
future_paper_needed: false
source_files:
- research/4d_2pn/paper/4d_2pn.tex
- 4d_2pn_summary.md
- 4d_2pn.tex
legacy_sources:
- 4d_2pn_summary.md
- 4d_2pn.tex
source_links:
- '[[FILE_2PN_FULL]]'
tex_anchor:
  file: research/4d_2pn/paper/4d_2pn.tex
  line: 1215
  heading_level: section
  heading: Comparable-mass 2PN ADM lift
  nearest_label:
    name: sec:adm-lift
    line: 1216
  nearby_labels:
  - name: sec:adm-lift
    line: 1216
  match_basis: semantic_heading_match
  match_score: 0.57
  confidence: medium
equation_ids:
- EQ_1PN_EIH_EQUALITY
- EQ_P0_TARGET
claim_ids:
- CLAIM_2PN_ADM_WITHIN_HIERARCHY
source_ids:
- FILE_2PN_FULL
outgoing_edges:
- target: PN_2_ADM_FULL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_P0_TARGET
  relation: EXPOSES_SUPPORT_STRUCTURE_FOR
  status: derivation
  note: 2PN P0/P2 support structure supplies quadrupole representation content.
- target: CLAIM_2PN_ADM_WITHIN_HIERARCHY
  relation: SUPPORTS_CLAIM
  status: exact_assembly_within_closure
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_2PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_2PN_ADM_EQUALITY.
- source: FILE_2PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_2PN_ADM_EQUALITY.
- source: EQ_1PN_EIH_EQUALITY
  relation: FEEDS
  status: derivation
  note: 2PN imports the 1PN conservative ledger.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/full_assembly_within_closure
- topic/pn_chain
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 2PN ADM equality

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_2PN_ADM_EQUALITY`  
> **Status:** `full_assembly_within_closure`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Conservative 2PN Lagrangian Legendre-transforms exactly to generic-frame ADM target inside declared hierarchy.

## Physical Meaning

Conservative 2PN Lagrangian Legendre-transforms exactly to generic-frame ADM target inside declared hierarchy.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `full_assembly_within_closure`
- Parent node: `PN_2_ADM_FULL`

## Equation

$$
H_2 = H_2PN^ADM
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_1PN_EIH_EQUALITY]]
- [[EQ_P0_TARGET]]

### Related claims
- [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_2PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[PN_2_ADM_FULL]] | Equation anchor belongs to or formalizes this graph node. |
| `EXPOSES_SUPPORT_STRUCTURE_FOR` | [[EQ_P0_TARGET]] | 2PN P0/P2 support structure supplies quadrupole representation content. |
| `SUPPORTS_CLAIM` | [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_2PN_FULL]] | Paper backlink block references EQ_2PN_ADM_EQUALITY. |
| `CONTAINS_EQUATION` | [[FILE_2PN_FULL]] | Source artifact contains or supports EQ_2PN_ADM_EQUALITY. |
| `FEEDS` | [[EQ_1PN_EIH_EQUALITY]] | 2PN imports the 1PN conservative ledger. |

## Source Anchors

### Source anchor notes
- [[FILE_2PN_FULL]]

### Source files
- `research/4d_2pn/paper/4d_2pn.tex`
- `4d_2pn_summary.md`
- `4d_2pn.tex`

### TeX anchor
- File: `research/4d_2pn/paper/4d_2pn.tex`
- Line: `1215`
- Heading: Comparable-mass 2PN ADM lift
- Nearest label: `sec:adm-lift` at line `1216`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
