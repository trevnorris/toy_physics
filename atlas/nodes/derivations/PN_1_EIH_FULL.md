---
id: PN_1_EIH_FULL
title: Full conservative 1PN EIH assembly
type: PN_backbone
layer: derivation
status: full_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Conservative 1PN two-body EIH ledger closed within declared hierarchy, not a solved moving-throat PDE theorem.
future_paper_needed: false
source_files:
- research/4d_1pn_full/paper/4d_1pn_full.tex
- 4d_1pn_full_summary.md
legacy_sources:
- 4d_1pn_full_summary.md
source_links:
- '[[FILE_1PN_BRIDGE]]'
- '[[FILE_1PN_FULL]]'
tex_anchor:
  file: research/4d_1pn_full/paper/4d_1pn_full.tex
  line: 2920
  heading_level: paragraph
  heading: Full conservative two-body (1 )PN sector.
  nearest_label:
    name: eq:fixed-full-eih
    line: 2938
  nearby_labels:
  - name: eq:fixed-full-eih
    line: 2938
  - name: eq:fixed-full-equality
    line: 2943
  match_basis: semantic_heading_match
  match_score: 0.555
  confidence: medium
equation_ids:
- EQ_1PN_EIH_EQUALITY
claim_ids:
- CLAIM_1PN_EIH_WITHIN_HIERARCHY
source_ids:
- FILE_1PN_BRIDGE
- FILE_1PN_FULL
outgoing_edges:
- target: PN_2_ADM_FULL
  relation: FEEDS
  status: closure
  note: 1PN ledger carried into 2PN ADM assembly.
- target: PN_2_ADM_FULL
  relation: FEEDS
  status: derivation
  note: 2PN assembly imports 1PN closure hierarchy and EIH ledger.
incoming_edges:
- source: EQ_1PN_EIH_EQUALITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_1PN_BRIDGE
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: FILE_1PN_FULL
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PN_0_NEWTONIAN_HOOK
  relation: FEEDS
  status: closure
  note: Newtonian backbone plus closures feeds 1PN EIH assembly.
- source: CLAIM_1PN_EIH_WITHIN_HIERARCHY
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/full_assembly_within_closure
- topic/moving_throat
- topic/pn_chain
- type/pn_backbone
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Full conservative 1PN EIH assembly

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_1_EIH_FULL`
> **Status:** `full_assembly_within_closure`
> **Layer:** `derivation`
> **Type:** `PN_backbone`

## Summary

Conservative 1PN two-body EIH ledger closed within declared hierarchy, not a solved moving-throat PDE theorem.

## Physical Meaning

Conservative 1PN two-body EIH ledger closed within declared hierarchy, not a solved moving-throat PDE theorem.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_backbone`
- Status: `full_assembly_within_closure`

## Equation

$$
L_1PN^derived=L_1PN^EIH
$$

$$
beta_1PN=3
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

### Related claims
- [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_1PN_BRIDGE]]
- [[FILE_1PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[PN_2_ADM_FULL]] | 1PN ledger carried into 2PN ADM assembly. |
| `FEEDS` | [[PN_2_ADM_FULL]] | 2PN assembly imports 1PN closure hierarchy and EIH ledger. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_1PN_EIH_EQUALITY]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_1PN_BRIDGE]] | File anchor documents this node. |
| `DOCUMENTS` | [[FILE_1PN_FULL]] | File anchor documents this node. |
| `FEEDS` | [[PN_0_NEWTONIAN_HOOK]] | Newtonian backbone plus closures feeds 1PN EIH assembly. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_1PN_BRIDGE]]
- [[FILE_1PN_FULL]]

### Source files
- `research/4d_1pn_full/paper/4d_1pn_full.tex`
- `4d_1pn_full_summary.md`

### TeX anchor
- File: `research/4d_1pn_full/paper/4d_1pn_full.tex`
- Line: `2920`
- Heading: Full conservative two-body (1 )PN sector.
- Nearest label: `eq:fixed-full-eih` at line `2938`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
