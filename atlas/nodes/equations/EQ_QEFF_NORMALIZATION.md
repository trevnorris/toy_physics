---
id: EQ_QEFF_NORMALIZATION
title: Thickness-controlled brane charge
type: equation
layer: equation_anchor
status: exact_within_zero_mode_normalization
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Observable brane charge strength after canonical zero-mode normalization.
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
  file: research/4d_em_fields/paper/4d_em_fields.tex
  line: 822
  heading_level: subsection
  heading: Canonical normalization and thickness-controlled brane charge
  nearest_label:
    name: sec:canonical_charge
    line: 823
  nearby_labels:
  - name: sec:canonical_charge
    line: 823
  - name: eq:qeff_canonical
    line: 856
  match_basis: semantic_heading_match
  match_score: 0.977
  confidence: high
math_ids:
- MATH_QSTAR_QEFF
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
status_firewall_ids:
- FIREWALL_QEFF_THICKNESS_NOT_BREATHING
source_ids:
- FILE_4D_PARENT
- FILE_EM_FIELDS
outgoing_edges:
- target: MATH_QSTAR_QEFF
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: ATOM_HYDROGEN_ZERO_MODE
  relation: FEEDS
  status: reduced
  note: Sets observed charge scale in the reduced atomic sector.
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: SUPPORTS_CLAIM
  status: reduced_sector_consequence
  note: Equation anchor supports this named claim.
- target: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: SUPPORTS_CLAIM
  status: exact_bookkeeping_firewall
  note: Equation anchor supports this named claim.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: SUPPORTS_CLAIM
  status: controlled_reduction
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_QEFF_NORMALIZATION.
- source: BACKLINK_ATOM_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_QEFF_NORMALIZATION.
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_QEFF_NORMALIZATION.
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_QEFF_NORMALIZATION.
- source: FILE_EM_FIELDS
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_QEFF_NORMALIZATION.
- source: FIREWALL_QEFF_THICKNESS_NOT_BREATHING
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_zero_mode_normalization
- topic/charge
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Thickness-controlled brane charge

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_QEFF_NORMALIZATION`  
> **Status:** `exact_within_zero_mode_normalization`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Observable brane charge strength after canonical zero-mode normalization.

## Physical Meaning

Observable brane charge strength after canonical zero-mode normalization.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_zero_mode_normalization`
- Parent node: [[MATH_QSTAR_QEFF]]

## Equation

$$
q_eff = q_*/sqrt(Z_int) = η_Q e_eff
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_QSTAR_QEFF]]

### Related equations
- none

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]]

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_QSTAR_QEFF]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[ATOM_HYDROGEN_ZERO_MODE]] | Sets observed charge scale in the reduced atomic sector. |
| `SUPPORTS_CLAIM` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references EQ_QEFF_NORMALIZATION. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_ATOM_WORK]] | Paper backlink block references EQ_QEFF_NORMALIZATION. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references EQ_QEFF_NORMALIZATION. |
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_QEFF_NORMALIZATION. |
| `CONTAINS_EQUATION` | [[FILE_EM_FIELDS]] | Source artifact contains or supports EQ_QEFF_NORMALIZATION. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]] | Firewall preserves this correct status boundary. |

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
- File: `research/4d_em_fields/paper/4d_em_fields.tex`
- Line: `822`
- Heading: Canonical normalization and thickness-controlled brane charge
- Nearest label: `sec:canonical_charge` at line `823`
- Match basis: `semantic_heading_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
