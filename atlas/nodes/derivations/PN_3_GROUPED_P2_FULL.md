---
id: PN_3_GROUPED_P2_FULL
title: Conservative 3PN grouped P2 assembly
type: PN_backbone
layer: derivation
status: full_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: 3PN conservative sector assigned in fixed ADM chart; grouped real P2 is live conservative payload.
future_paper_needed: false
source_files:
- research/4d_3pn/paper/4d_3pn.tex
- 4d_3pn_summary.md
legacy_sources:
- 4d_3pn_summary.md
source_links:
- '[[FILE_3PN_FULL]]'
tex_anchor:
  file: research/4d_3pn/paper/4d_3pn.tex
  line: 3248
  heading_level: subsection
  heading: Why 3PN is the first conservative grouped-P 2 gate
  nearest_label:
    name: sec:interface-why
    line: 3249
  nearby_labels:
  - name: sec:interface-why
    line: 3249
  - name: eq:interface-com-split-raw
    line: 3269
  match_basis: semantic_heading_match
  match_score: 0.603
  confidence: medium
equation_ids:
- EQ_3PN_SPLIT
claim_ids:
- CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
source_ids:
- FILE_3PN_FULL
outgoing_edges:
- target: PN_4_LOCAL_TAIL
  relation: FEEDS
  status: derivation
  note: 3PN grouped P2 closure feeds 4PN local/tail interface.
- target: PN_4_LOCAL_TAIL
  relation: FEEDS
  status: closure
  note: Lower conservative ledger carried into local 4PN assembly.
incoming_edges:
- source: EQ_3PN_SPLIT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_3PN_FULL
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PN_2_5_QUAD_NARROWING
  relation: FEEDS
  status: derivation
  note: 2.5PN identifies grouped real P2 as higher-order payload.
- source: PN_2_5_QUAD_NARROWING
  relation: FEEDS
  status: closure
  note: Grouped real P2 bundle becomes 3PN conservative payload.
- source: CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/full_assembly_within_closure
- topic/pn_chain
- type/pn_backbone
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Conservative 3PN grouped P2 assembly

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_3_GROUPED_P2_FULL`
> **Status:** `full_assembly_within_closure`
> **Layer:** `derivation`
> **Type:** `PN_backbone`

## Summary

3PN conservative sector assigned in fixed ADM chart; grouped real P2 is live conservative payload.

## Physical Meaning

3PN conservative sector assigned in fixed ADM chart; grouped real P2 is live conservative payload.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_backbone`
- Status: `full_assembly_within_closure`

## Equation

$$
P20⊕P21⊕P22 middle block
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_3PN_SPLIT]]

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
| `FEEDS` | [[PN_4_LOCAL_TAIL]] | 3PN grouped P2 closure feeds 4PN local/tail interface. |
| `FEEDS` | [[PN_4_LOCAL_TAIL]] | Lower conservative ledger carried into local 4PN assembly. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_3PN_SPLIT]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_3PN_FULL]] | File anchor documents this node. |
| `FEEDS` | [[PN_2_5_QUAD_NARROWING]] | 2.5PN identifies grouped real P2 as higher-order payload. |
| `FEEDS` | [[PN_2_5_QUAD_NARROWING]] | Grouped real P2 bundle becomes 3PN conservative payload. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_3PN_FULL]]

### Source files
- `research/4d_3pn/paper/4d_3pn.tex`
- `4d_3pn_summary.md`

### TeX anchor
- File: `research/4d_3pn/paper/4d_3pn.tex`
- Line: `3248`
- Heading: Why 3PN is the first conservative grouped-P 2 gate
- Nearest label: `sec:interface-why` at line `3249`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
