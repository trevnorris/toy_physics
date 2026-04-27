---
id: EQ_1PN_EIH_EQUALITY
title: 1PN EIH equality
type: equation
layer: equation_anchor
status: full_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Full conservative two-body 1PN assembly equals EIH inside declared closure hierarchy.
future_paper_needed: false
source_files:
- research/4d_1pn_full/paper/4d_1pn_full.tex
- 4d_1pn_full_summary.md
legacy_sources:
- 4d_1pn_full_summary.md
source_links:
- '[[FILE_1PN_FULL]]'
tex_anchor:
  file: research/4d_1pn_full/paper/4d_1pn_full.tex
  line: 2403
  heading_level: subsection
  heading: Standard EIH target and exact equality
  nearest_label:
    name: sec:eih-target
    line: 2404
  nearby_labels:
  - name: sec:eih-target
    line: 2404
  - name: eq:eih-target
    line: 2423
  - name: eq:eih-residual-zero
    line: 2431
  - name: tab:coefficient-provenance
    line: 2446
  match_basis: semantic_heading_match
  match_score: 0.581
  confidence: medium
equation_ids:
- EQ_2PN_ADM_EQUALITY
claim_ids:
- CLAIM_1PN_EIH_WITHIN_HIERARCHY
source_ids:
- FILE_1PN_FULL
outgoing_edges:
- target: PN_1_EIH_FULL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_2PN_ADM_EQUALITY
  relation: FEEDS
  status: derivation
  note: 2PN imports the 1PN conservative ledger.
- target: CLAIM_1PN_EIH_WITHIN_HIERARCHY
  relation: SUPPORTS_CLAIM
  status: exact_assembly_within_closure
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_1PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_1PN_EIH_EQUALITY.
- source: FILE_1PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_1PN_EIH_EQUALITY.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/full_assembly_within_closure
- topic/pn_chain
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 1PN EIH equality

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_1PN_EIH_EQUALITY`  
> **Status:** `full_assembly_within_closure`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Full conservative two-body 1PN assembly equals EIH inside declared closure hierarchy.

## Physical Meaning

Full conservative two-body 1PN assembly equals EIH inside declared closure hierarchy.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `full_assembly_within_closure`
- Parent node: `PN_1_EIH_FULL`

## Equation

$$
L_1PN^derived = L_1PN^EIH
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_2PN_ADM_EQUALITY]]

### Related claims
- [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_1PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[PN_1_EIH_FULL]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_2PN_ADM_EQUALITY]] | 2PN imports the 1PN conservative ledger. |
| `SUPPORTS_CLAIM` | [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_1PN_FULL]] | Paper backlink block references EQ_1PN_EIH_EQUALITY. |
| `CONTAINS_EQUATION` | [[FILE_1PN_FULL]] | Source artifact contains or supports EQ_1PN_EIH_EQUALITY. |

## Source Anchors

### Source anchor notes
- [[FILE_1PN_FULL]]

### Source files
- `research/4d_1pn_full/paper/4d_1pn_full.tex`
- `4d_1pn_full_summary.md`

### TeX anchor
- File: `research/4d_1pn_full/paper/4d_1pn_full.tex`
- Line: `2403`
- Heading: Standard EIH target and exact equality
- Nearest label: `sec:eih-target` at line `2404`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
