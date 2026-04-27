---
id: EQ_P0_TARGET
title: Invariant quadrupole normalization target
type: equation
layer: equation_anchor
status: open_actual_branch_data
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Central remaining 2.5PN/4PN moving-throat normalization target.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_005.tex
- research/4d_2_5pn/paper/4d_2_5pn.tex
- research/4d_4pn/paper/4d_4pn.tex
- notes/pde_audit_full.md
- moving_throat_pde_stage005_grouped_p2_normalization_bridge.md
- 4d_2_5pn_summary.md
- 4d_4pn_summary.md
- pde_audit_full.md
legacy_sources:
- moving_throat_pde_stage005_grouped_p2_normalization_bridge.md
- 4d_2_5pn_summary.md
- 4d_4pn_summary.md
- pde_audit_full.md
source_links:
- '[[FILE_2_5PN]]'
- '[[FILE_4PN_FULL]]'
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
tex_anchor:
  file: research/4d_2_5pn/paper/4d_2_5pn.tex
  line: 4995
  heading_level: section
  heading: Quadrupole normalization package
  nearest_label:
    name: app:quadrupole-normalization-package
    line: 4996
  nearby_labels:
  - name: app:quadrupole-normalization-package
    line: 4996
  match_basis: semantic_heading_match
  match_score: 0.655
  confidence: medium
equation_ids:
- EQ_2PN_ADM_EQUALITY
- EQ_4PN_TAIL_BRIDGE
- EQ_COMPACT_L2_FINGERPRINT
- EQ_FULL_BUNDLE_TARGET_SURFACE
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
- CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE5_P0_TARGET
- CLAIM_STAGE8_TO_14_SELECTED_BRANCH
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
source_ids:
- FILE_2_5PN
- FILE_4PN_FULL
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: OPEN_QUAD_NORMALIZATION
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_FULL_BUNDLE_TARGET_SURFACE
  relation: IS_COMPONENT_OF
  status: conditional
  note: Universal quadrupole target is one part of the full-bundle surface.
- target: EQ_4PN_TAIL_BRIDGE
  relation: REQUIRED_BY
  status: conditional
  note: The 4PN tail bridge uses the same quadrupole normalization.
- target: CLAIM_25PN_QUAD_NARROWING
  relation: SUPPORTS_CLAIM
  status: conditional_theorem_open_normalization
  note: Equation anchor supports this named claim.
- target: CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
  relation: SUPPORTS_CLAIM
  status: local_closed_tail_conditional
  note: Equation anchor supports this named claim.
- target: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: SUPPORTS_CLAIM
  status: paper_facing_ontology_discipline
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE5_P0_TARGET
  relation: SUPPORTS_CLAIM
  status: exact_within_grouped_bridge
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE8_TO_14_SELECTED_BRANCH
  relation: SUPPORTS_CLAIM
  status: exact_within_selected_reduced_branch
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_25PN
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_P0_TARGET.
- source: FILE_2_5PN
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_P0_TARGET.
- source: FILE_4PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_P0_TARGET.
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_P0_TARGET.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_P0_TARGET.
- source: EQ_2PN_ADM_EQUALITY
  relation: EXPOSES_SUPPORT_STRUCTURE_FOR
  status: derivation
  note: 2PN P0/P2 support structure supplies quadrupole representation content.
- source: EQ_COMPACT_L2_FINGERPRINT
  relation: FEEDS
  status: conditional
  note: The iω^5 coefficient sets the universal quadrupole target.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/open_actual_branch_data
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Invariant quadrupole normalization target

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_P0_TARGET`  
> **Status:** `open_actual_branch_data`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Central remaining 2.5PN/4PN moving-throat normalization target.

## Physical Meaning

Central remaining 2.5PN/4PN moving-throat normalization target.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `open_actual_branch_data`
- Parent node: `OPEN_QUAD_NORMALIZATION`

## Equation

$$
mhat_0² P_0 = 54 G c_s⁵/(5 a⁵ c⁵)
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
- [[EQ_4PN_TAIL_BRIDGE]]
- [[EQ_COMPACT_L2_FINGERPRINT]]
- [[EQ_FULL_BUNDLE_TARGET_SURFACE]]

### Related claims
- [[CLAIM_25PN_QUAD_NARROWING]]
- [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE5_P0_TARGET]]
- [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- none

### Source anchors
- [[FILE_2_5PN]]
- [[FILE_4PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[OPEN_QUAD_NORMALIZATION]] | Equation anchor belongs to or formalizes this graph node. |
| `IS_COMPONENT_OF` | [[EQ_FULL_BUNDLE_TARGET_SURFACE]] | Universal quadrupole target is one part of the full-bundle surface. |
| `REQUIRED_BY` | [[EQ_4PN_TAIL_BRIDGE]] | The 4PN tail bridge uses the same quadrupole normalization. |
| `SUPPORTS_CLAIM` | [[CLAIM_25PN_QUAD_NARROWING]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE5_P0_TARGET]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_25PN]] | Paper backlink block references EQ_P0_TARGET. |
| `CONTAINS_EQUATION` | [[FILE_2_5PN]] | Source artifact contains or supports EQ_P0_TARGET. |
| `CONTAINS_EQUATION` | [[FILE_4PN_FULL]] | Source artifact contains or supports EQ_P0_TARGET. |
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_P0_TARGET. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_P0_TARGET. |
| `EXPOSES_SUPPORT_STRUCTURE_FOR` | [[EQ_2PN_ADM_EQUALITY]] | 2PN P0/P2 support structure supplies quadrupole representation content. |
| `FEEDS` | [[EQ_COMPACT_L2_FINGERPRINT]] | The iω^5 coefficient sets the universal quadrupole target. |

## Source Anchors

### Source anchor notes
- [[FILE_2_5PN]]
- [[FILE_4PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/pde_ledger/paper/stages/stage_005.tex`
- `research/4d_2_5pn/paper/4d_2_5pn.tex`
- `research/4d_4pn/paper/4d_4pn.tex`
- `notes/pde_audit_full.md`
- `moving_throat_pde_stage005_grouped_p2_normalization_bridge.md`
- `4d_2_5pn_summary.md`
- `4d_4pn_summary.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/4d_2_5pn/paper/4d_2_5pn.tex`
- Line: `4995`
- Heading: Quadrupole normalization package
- Nearest label: `app:quadrupole-normalization-package` at line `4996`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
