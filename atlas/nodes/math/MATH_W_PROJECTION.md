---
id: MATH_W_PROJECTION
title: Projection kernel W(w)
type: kernel
layer: math_object
status: exact_measurement_definition
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Defines what brane observer measures; distinct from localization profile Z.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_plasma/paper/4d_plasma.tex
- 4d_summary.md
- 4d_plasma_summary.md
legacy_sources:
- 4d_summary.md
- 4d_plasma_summary.md
source_links:
- '[[FILE_PLASMA]]'
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 1457
  heading_level: paragraph
  heading: Choice of projection weight W(w).
  nearest_label: null
  nearby_labels: []
  match_basis: semantic_heading_match
  match_score: 0.637
  confidence: medium
physical_ids:
- PHYS_BRANE_OBSERVER
- PHYS_REG_OBSERVER_PROJECTION
math_ids:
- MATH_PROJECTED_CONTINUITY
status_firewall_ids:
- FIREWALL_PROJECTION_NOT_REDUCTION
source_ids:
- FILE_PLASMA
outgoing_edges:
- target: MATH_PROJECTED_CONTINUITY
  relation: DERIVES
  status: exact
  note: Projection of exact bulk continuity gives leakage source.
incoming_edges:
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_W_PROJECTION.
- source: PHYS_BRANE_OBSERVER
  relation: DEFINED_BY
  status: exact
  note: Projection kernel defines brane observables.
- source: FILE_PLASMA
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PHYS_REG_OBSERVER_PROJECTION
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FIREWALL_PROJECTION_NOT_REDUCTION
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_PROJECTION_EQUALS_REDUCTION
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MATH_W_PROJECTION.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_measurement_definition
- topic/projection
- type/kernel
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Projection kernel W(w)

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_W_PROJECTION`
> **Status:** `exact_measurement_definition`
> **Layer:** `math_object`
> **Type:** `kernel`

## Summary

Defines what brane observer measures; distinct from localization profile Z.

## Physical Meaning

Defines what brane observer measures; distinct from localization profile Z.

## Mathematical Role

- Layer: `math_object`
- Type: `kernel`
- Status: `exact_measurement_definition`

## Equation

$$
⟨Q⟩_W=∫WQdw
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_OBSERVER]]
- [[PHYS_REG_OBSERVER_PROJECTION]]

### Related math nodes
- [[MATH_PROJECTED_CONTINUITY]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- [[FIREWALL_PROJECTION_NOT_REDUCTION]]

### Source anchors
- [[FILE_PLASMA]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DERIVES` | [[MATH_PROJECTED_CONTINUITY]] | Projection of exact bulk continuity gives leakage source. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references MATH_W_PROJECTION. |
| `DEFINED_BY` | [[PHYS_BRANE_OBSERVER]] | Projection kernel defines brane observables. |
| `DOCUMENTS` | [[FILE_PLASMA]] | File anchor documents this node. |
| `LINKS_TO` | [[PHYS_REG_OBSERVER_PROJECTION]] | Physical register entry links to graph object. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_PROJECTION_NOT_REDUCTION]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_PROJECTION_EQUALS_REDUCTION]] | Negative query starts from MATH_W_PROJECTION. |

## Source Anchors

### Source anchor notes
- [[FILE_PLASMA]]

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_plasma/paper/4d_plasma.tex`
- `4d_summary.md`
- `4d_plasma_summary.md`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `1457`
- Heading: Choice of projection weight W(w).
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
