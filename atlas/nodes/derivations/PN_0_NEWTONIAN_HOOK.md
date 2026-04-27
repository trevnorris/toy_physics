---
id: PN_0_NEWTONIAN_HOOK
title: Newtonian / 0PN inverse-square hook
type: PN_backbone
layer: derivation
status: controlled_reduction
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Poisson hook and worldtube reduction produce Newtonian brane-facing scalar sector.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_1pn_full/paper/4d_1pn_full.tex
- 4d_summary.md
- 4d_1pn_full_summary.md
legacy_sources:
- 4d_summary.md
- 4d_1pn_full_summary.md
tex_anchor:
  file: research/4d_1pn_full/paper/4d_1pn_full.tex
  line: 575
  heading_level: subsection
  heading: The exact longitudinal identity and the Poisson hook
  nearest_label:
    name: eq:dictionary-inverse-square
    line: 575
  nearby_labels:
  - name: eq:dictionary-inverse-square
    line: 575
  match_basis: semantic_label_match
  match_score: 0.515
  confidence: medium
math_ids:
- MATH_POISSON_HOOK
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
outgoing_edges:
- target: PN_1_EIH_FULL
  relation: FEEDS
  status: closure
  note: Newtonian backbone plus closures feeds 1PN EIH assembly.
incoming_edges:
- source: MATH_POISSON_HOOK
  relation: FEEDS
  status: controlled
  note: Poisson hook feeds Newtonian point-particle limit.
- source: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_06_07_NEWTONIAN_UNIVERSALITY
  relation: VALIDATES_BACKBONE
  status: controlled
  note: Projected continuity and Poisson hook support same-G Newtonian ledger.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/controlled_reduction
- topic/pn_chain
- topic/projection
- type/pn_backbone
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Newtonian / 0PN inverse-square hook

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_0_NEWTONIAN_HOOK`  
> **Status:** `controlled_reduction`  
> **Layer:** `derivation`  
> **Type:** `PN_backbone`

## Summary

Poisson hook and worldtube reduction produce Newtonian brane-facing scalar sector.

## Physical Meaning

Poisson hook and worldtube reduction produce Newtonian brane-facing scalar sector.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_backbone`
- Status: `controlled_reduction`

## Equation

$$
inverse-square longitudinal field
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_POISSON_HOOK]]

### Related equations
- none

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[PN_1_EIH_FULL]] | Newtonian backbone plus closures feeds 1PN EIH assembly. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MATH_POISSON_HOOK]] | Poisson hook feeds Newtonian point-particle limit. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Claim feeds this downstream object, output, or open gate. |
| `VALIDATES_BACKBONE` | [[MT_V2_06_07_NEWTONIAN_UNIVERSALITY]] | Projected continuity and Poisson hook support same-G Newtonian ledger. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_1pn_full/paper/4d_1pn_full.tex`
- `4d_summary.md`
- `4d_1pn_full_summary.md`

### TeX anchor
- File: `research/4d_1pn_full/paper/4d_1pn_full.tex`
- Line: `575`
- Heading: The exact longitudinal identity and the Poisson hook
- Nearest label: `eq:dictionary-inverse-square` at line `575`
- Match basis: `semantic_label_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
