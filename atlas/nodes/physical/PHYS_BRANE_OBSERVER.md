---
id: PHYS_BRANE_OBSERVER
title: Operational brane observer
type: observer
layer: physical_ontology
status: exact_definition
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: A 3D observer accesses projected quantities, not the full bulk state.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_plasma/paper/4d_plasma.tex
- 4d_summary.md
- 4d_plasma_summary.md
legacy_sources:
- 4d_summary.md
- 4d_plasma_summary.md
tex_anchor:
  file: research/4d_plasma/paper/4d_plasma.tex
  line: 279
  heading_level: paragraph
  heading: Operational brane projection W(w).
  nearest_label:
    name: eq:W_normalized
    line: 287
  nearby_labels:
  - name: eq:W_normalized
    line: 287
  - name: eq:brane_projection
    line: 292
  match_basis: semantic_heading_match
  match_score: 0.68
  confidence: medium
physical_ids:
- PHYS_REG_OBSERVER_PROJECTION
math_ids:
- MATH_W_PROJECTION
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
outgoing_edges:
- target: MATH_W_PROJECTION
  relation: DEFINED_BY
  status: exact
  note: Projection kernel defines brane observables.
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: GROUNDS_PHYSICAL_MEANING
  status: reduced_sector_consequence
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_projection_plus_controlled_hook
  note: Physical ontology object grounded by this claim.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: GROUNDS_PHYSICAL_MEANING
  status: controlled_reduction
  note: Physical ontology object grounded by this claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_BRANE_OBSERVER.
- source: BACKLINK_ATOM_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_BRANE_OBSERVER.
- source: PHYS_REG_OBSERVER_PROJECTION
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_definition
- topic/projection
- type/observer
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Operational brane observer

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_BRANE_OBSERVER`
> **Status:** `exact_definition`
> **Layer:** `physical_ontology`
> **Type:** `observer`

## Summary

A 3D observer accesses projected quantities, not the full bulk state.

## Physical Meaning

A 3D observer accesses projected quantities, not the full bulk state.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `observer`
- Status: `exact_definition`

## Equation

$$
projection kernel W(w)
$$

$$
P_W[Q]=int W Q dw
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_REG_OBSERVER_PROJECTION]]

### Related math nodes
- [[MATH_W_PROJECTION]]

### Related equations
- none

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DEFINED_BY` | [[MATH_W_PROJECTION]] | Projection kernel defines brane observables. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Physical ontology object grounded by this claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references PHYS_BRANE_OBSERVER. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_ATOM_WORK]] | Paper backlink block references PHYS_BRANE_OBSERVER. |
| `LINKS_TO` | [[PHYS_REG_OBSERVER_PROJECTION]] | Physical register entry links to graph object. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_plasma/paper/4d_plasma.tex`
- `4d_summary.md`
- `4d_plasma_summary.md`

### TeX anchor
- File: `research/4d_plasma/paper/4d_plasma.tex`
- Line: `279`
- Heading: Operational brane projection W(w).
- Nearest label: `eq:W_normalized` at line `287`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
