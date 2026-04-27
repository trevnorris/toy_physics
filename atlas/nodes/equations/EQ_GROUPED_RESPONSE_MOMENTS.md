---
id: EQ_GROUPED_RESPONSE_MOMENTS
title: Grouped response moment conversion
type: equation
layer: equation_anchor
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Bridge from conservative operator moments to normalized grouped response moments.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_005.tex
- notes/pde_audit_full.md
- moving_throat_pde_stage005_grouped_p2_normalization_bridge.md
- pde_audit_full.md
legacy_sources:
- moving_throat_pde_stage005_grouped_p2_normalization_bridge.md
- pde_audit_full.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
tex_anchor:
  file: research/pde_ledger/paper/stages/stage_005.tex
  line: 15
  heading_level: paragraph
  heading: 1. Conservative operator and normalized response.
  nearest_label:
    name: eq:app-stage005-d-operator
    line: 18
  nearby_labels:
  - name: eq:app-stage005-d-operator
    line: 18
  - name: eq:app-stage005-y-response
    line: 23
  - name: eq:app-stage005-u2u4
    line: 29
  match_basis: semantic_heading_match
  match_score: 0.589
  confidence: medium
equation_ids:
- EQ_BDG_SCHUR_KERNEL
- EQ_FULL_BUNDLE_TARGET_SURFACE
claim_ids:
- CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
- CLAIM_STAGE5_P0_TARGET
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MT_STAGE5_GROUPED_P2_BRIDGE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_FULL_BUNDLE_TARGET_SURFACE
  relation: FEEDS
  status: reduced
  note: Response moments enter one-pole and full-bundle target surface.
- target: CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
  relation: SUPPORTS_CLAIM
  status: exact_assembly_within_closure
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE5_P0_TARGET
  relation: SUPPORTS_CLAIM
  status: exact_within_grouped_bridge
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_GROUPED_RESPONSE_MOMENTS.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_GROUPED_RESPONSE_MOMENTS.
- source: EQ_BDG_SCHUR_KERNEL
  relation: FEEDS
  status: reduced
  note: BdG-dressed operator moments enter normalized response.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_reduced_bundle
- topic/moving_throat
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Grouped response moment conversion

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_GROUPED_RESPONSE_MOMENTS`  
> **Status:** `exact_within_reduced_bundle`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Bridge from conservative operator moments to normalized grouped response moments.

## Physical Meaning

Bridge from conservative operator moments to normalized grouped response moments.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_reduced_bundle`
- Parent node: `MT_STAGE5_GROUPED_P2_BRIDGE`

## Equation

$$
u_2^(A)=-D_A2/D_A0,   u_4^(A)=(D_A2²-D_A0D_A4)/D_A0²
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_BDG_SCHUR_KERNEL]]
- [[EQ_FULL_BUNDLE_TARGET_SURFACE]]

### Related claims
- [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]]
- [[CLAIM_STAGE5_P0_TARGET]]

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
| `ANCHORS` | [[MT_STAGE5_GROUPED_P2_BRIDGE]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_FULL_BUNDLE_TARGET_SURFACE]] | Response moments enter one-pole and full-bundle target surface. |
| `SUPPORTS_CLAIM` | [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE5_P0_TARGET]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_GROUPED_RESPONSE_MOMENTS. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_GROUPED_RESPONSE_MOMENTS. |
| `FEEDS` | [[EQ_BDG_SCHUR_KERNEL]] | BdG-dressed operator moments enter normalized response. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/pde_ledger/paper/stages/stage_005.tex`
- `notes/pde_audit_full.md`
- `moving_throat_pde_stage005_grouped_p2_normalization_bridge.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/pde_ledger/paper/stages/stage_005.tex`
- Line: `15`
- Heading: 1. Conservative operator and normalized response.
- Nearest label: `eq:app-stage005-d-operator` at line `18`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
