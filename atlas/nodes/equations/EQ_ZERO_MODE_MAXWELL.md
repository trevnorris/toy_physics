---
id: EQ_ZERO_MODE_MAXWELL
title: Zero-mode brane Maxwell reduction
type: equation
layer: equation_anchor
status: controlled_reduction
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Controlled far-field brane reduction of localized 4+1 Maxwell.
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
- research/4d/paper/4d.tex
- 4d_em_fields_summary.md
- 4d_summary.md
legacy_sources:
- 4d_em_fields_summary.md
- 4d_summary.md
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_EM_FIELDS]]'
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 1995
  heading_level: paragraph
  heading: EM localization reduction to a brane Maxwell theory.
  nearest_label: null
  nearby_labels: []
  match_basis: semantic_heading_match
  match_score: 0.528
  confidence: medium
math_ids:
- MATH_ZERO_MODE_MAXWELL
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
- CLAIM_PARENT_ACTION_CURRENT_EXACT
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
source_ids:
- FILE_4D_PARENT
- FILE_EM_FIELDS
outgoing_edges:
- target: MATH_ZERO_MODE_MAXWELL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: ATOM_HYDROGEN_ZERO_MODE
  relation: FEEDS
  status: reduced
  note: Provides the Coulomb sector used in first-pass atomic reduction.
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: SUPPORTS_CLAIM
  status: reduced_sector_consequence
  note: Equation anchor supports this named claim.
- target: CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
  relation: SUPPORTS_CLAIM
  status: safe_interpretation_or_structural_patch
  note: Equation anchor supports this named claim.
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: SUPPORTS_CLAIM
  status: controlled_reduction
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_ZERO_MODE_MAXWELL.
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_ZERO_MODE_MAXWELL.
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_ZERO_MODE_MAXWELL.
- source: FILE_EM_FIELDS
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_ZERO_MODE_MAXWELL.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/controlled_reduction
- topic/maxwell
- topic/projection
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Zero-mode brane Maxwell reduction

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_ZERO_MODE_MAXWELL`
> **Status:** `controlled_reduction`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Controlled far-field brane reduction of localized 4+1 Maxwell.

## Physical Meaning

Controlled far-field brane reduction of localized 4+1 Maxwell.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `controlled_reduction`
- Parent node: [[MATH_ZERO_MODE_MAXWELL]]

## Equation

$$
∂_μ F^{μν} = μ0_eff J_eff^ν,    μ0_eff = μ0/Z_int
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_ZERO_MODE_MAXWELL]]

### Related equations
- none

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]]
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_ZERO_MODE_MAXWELL]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[ATOM_HYDROGEN_ZERO_MODE]] | Provides the Coulomb sector used in first-pass atomic reduction. |
| `SUPPORTS_CLAIM` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references EQ_ZERO_MODE_MAXWELL. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references EQ_ZERO_MODE_MAXWELL. |
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_ZERO_MODE_MAXWELL. |
| `CONTAINS_EQUATION` | [[FILE_EM_FIELDS]] | Source artifact contains or supports EQ_ZERO_MODE_MAXWELL. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`
- `research/4d/paper/4d.tex`
- `4d_em_fields_summary.md`
- `4d_summary.md`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `1995`
- Heading: EM localization reduction to a brane Maxwell theory.
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
