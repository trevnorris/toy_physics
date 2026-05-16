---
id: MATH_GNLS_PSI
title: Gauged 4D GNLS field
type: field_equation
layer: math_object
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Bulk matter/order parameter sector with exact current and continuity.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- 4d_summary.md
legacy_sources:
- 4d_summary.md
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 365
  heading_level: subsection
  heading: 'Matter sector: gauged 4D GNLS with confinement'
  nearest_label:
    name: sec:action-matter
    line: 366
  nearby_labels:
  - name: sec:action-matter
    line: 366
  - name: eq:Lpsi
    line: 384
  match_basis: semantic_heading_match
  match_score: 0.636
  confidence: medium
physical_ids:
- PHYS_SUPERFLUID_MEDIUM
math_ids:
- MATH_PARENT_ACTION_CURRENT
equation_ids:
- EQ_GNLS_BULK
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
incoming_edges:
- source: EQ_GNLS_BULK
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: MATH_PARENT_ACTION_CURRENT
  relation: DERIVES
  status: exact
  note: Variation gives GNLS matter equation and continuity.
- source: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: FEEDS_OR_STATUS_OF
  status: exact_parent_current
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_SUPERFLUID_MEDIUM
  relation: REPRESENTED_BY
  status: exact
  note: Superfluid medium represented by gauged GNLS field and density/current.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact
- type/field_equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Gauged 4D GNLS field

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_GNLS_PSI`
> **Status:** `exact`
> **Layer:** `math_object`
> **Type:** `field_equation`

## Summary

Bulk matter/order parameter sector with exact current and continuity.

## Physical Meaning

Bulk matter/order parameter sector with exact current and continuity.

## Mathematical Role

- Layer: `math_object`
- Type: `field_equation`
- Status: `exact`

## Equation

$$
iℏD_tψ=[-ℏ²D_iD_i/(2m)+V_conf+h(ρ)]ψ
$$

$$
∂_tρ+∂_i j^i=0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_SUPERFLUID_MEDIUM]]

### Related math nodes
- [[MATH_PARENT_ACTION_CURRENT]]

### Related equations
- [[EQ_GNLS_BULK]]

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_GNLS_BULK]] | Equation anchor belongs to or formalizes this graph node. |
| `DERIVES` | [[MATH_PARENT_ACTION_CURRENT]] | Variation gives GNLS matter equation and continuity. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Claim feeds this downstream object, output, or open gate. |
| `REPRESENTED_BY` | [[PHYS_SUPERFLUID_MEDIUM]] | Superfluid medium represented by gauged GNLS field and density/current. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `4d_summary.md`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `365`
- Heading: Matter sector: gauged 4D GNLS with confinement
- Nearest label: `sec:action-matter` at line `366`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
