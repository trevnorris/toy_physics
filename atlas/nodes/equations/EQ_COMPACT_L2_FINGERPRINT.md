---
id: EQ_COMPACT_L2_FINGERPRINT
title: Compact outgoing l=2 fingerprint
type: equation
layer: equation_anchor
status: exact_reduced
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Universal outgoing quadrupole low-frequency fingerprint in normalized response language.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_004.tex
- notes/pde_audit_full.md
- moving_throat_pde_stage004_maxwell_mixed_response.md
- pde_audit_full.md
legacy_sources:
- moving_throat_pde_stage004_maxwell_mixed_response.md
- pde_audit_full.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
tex_anchor:
  file: research/pde_ledger/paper/stages/stage_004.tex
  line: 113
  heading_level: paragraph
  heading: 5. Compact outgoing (l=2 ) fingerprint.
  nearest_label:
    name: eq:app-stage004-outgoing-fingerprint
    line: 116
  nearby_labels:
  - name: eq:app-stage004-outgoing-fingerprint
    line: 116
  - name: eq:app-stage004-gamma-port
    line: 125
  - name: eq:app-stage004-wall-odd
    line: 130
  match_basis: semantic_heading_match
  match_score: 1.0
  confidence: high
math_ids:
- MATH_COMPACT_L2_OUTGOING_FINGERPRINT
equation_ids:
- EQ_MAXWELL_MIXED_TRANSFER
- EQ_P0_TARGET
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
- CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_COMPACT_L2_OUTGOING_FINGERPRINT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_P0_TARGET
  relation: FEEDS
  status: conditional
  note: The iω^5 coefficient sets the universal quadrupole target.
- target: CLAIM_25PN_QUAD_NARROWING
  relation: SUPPORTS_CLAIM
  status: conditional_theorem_open_normalization
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_25PN
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_COMPACT_L2_FINGERPRINT.
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_COMPACT_L2_FINGERPRINT.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_COMPACT_L2_FINGERPRINT.
- source: EQ_MAXWELL_MIXED_TRANSFER
  relation: MULTIPLIES
  status: reduced
  note: Mixed transfer carries the outgoing fingerprint to the wall branch.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_reduced
- topic/maxwell
- topic/moving_throat
- topic/quadrupole
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Compact outgoing l=2 fingerprint

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_COMPACT_L2_FINGERPRINT`  
> **Status:** `exact_reduced`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Universal outgoing quadrupole low-frequency fingerprint in normalized response language.

## Physical Meaning

Universal outgoing quadrupole low-frequency fingerprint in normalized response language.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_reduced`
- Parent node: [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]]

## Equation

$$
Ŷ_2^out=1+a²ω²/(9c_s²)+4a⁴ω⁴/(81c_s⁴)+i a⁵ω⁵/(27c_s⁵)+O(ω⁶)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]]

### Related equations
- [[EQ_MAXWELL_MIXED_TRANSFER]]
- [[EQ_P0_TARGET]]

### Related claims
- [[CLAIM_25PN_QUAD_NARROWING]]
- [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_P0_TARGET]] | The iω^5 coefficient sets the universal quadrupole target. |
| `SUPPORTS_CLAIM` | [[CLAIM_25PN_QUAD_NARROWING]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_25PN]] | Paper backlink block references EQ_COMPACT_L2_FINGERPRINT. |
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_COMPACT_L2_FINGERPRINT. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_COMPACT_L2_FINGERPRINT. |
| `MULTIPLIES` | [[EQ_MAXWELL_MIXED_TRANSFER]] | Mixed transfer carries the outgoing fingerprint to the wall branch. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/pde_ledger/paper/stages/stage_004.tex`
- `notes/pde_audit_full.md`
- `moving_throat_pde_stage004_maxwell_mixed_response.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/pde_ledger/paper/stages/stage_004.tex`
- Line: `113`
- Heading: 5. Compact outgoing (l=2 ) fingerprint.
- Nearest label: `eq:app-stage004-outgoing-fingerprint` at line `116`
- Match basis: `semantic_heading_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
