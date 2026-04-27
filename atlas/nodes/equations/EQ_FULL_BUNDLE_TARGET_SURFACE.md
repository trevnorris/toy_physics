---
id: EQ_FULL_BUNDLE_TARGET_SURFACE
title: Isotropic full-bundle target surface
type: equation
layer: equation_anchor
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Combined one-pole and quadrupole-normalization target surface.
future_paper_needed: false
source_files:
- notes/5pn/5pn_notes_full.md
- notes/pde_audit_full.md
- 5pn_notes_full.md
- pde_audit_full.md
legacy_sources:
- 5pn_notes_full.md
- pde_audit_full.md
source_links:
- '[[FILE_5PN_FULL]]'
- '[[FILE_PDE_AUDIT]]'
math_ids:
- MATH_FULL_BUNDLE_TARGET_SURFACE
equation_ids:
- EQ_GROUPED_RESPONSE_MOMENTS
- EQ_P0_TARGET
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_BRANCH_EXPORTER_REQUIRED
- CLAIM_PACKET_A_PACKET_B_SPLIT
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE6_FULL_BUNDLE_RATIO
- CLAIM_STAGE8_TO_14_SELECTED_BRANCH
source_ids:
- FILE_5PN_FULL
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_FULL_BUNDLE_TARGET_SURFACE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Equation anchor supports this named claim.
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: SUPPORTS_CLAIM
  status: open_actual_branch
  note: Equation anchor supports this named claim.
- target: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: SUPPORTS_CLAIM
  status: open_branch_packets
  note: Equation anchor supports this named claim.
- target: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: SUPPORTS_CLAIM
  status: paper_facing_ontology_discipline
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE6_FULL_BUNDLE_RATIO
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE8_TO_14_SELECTED_BRANCH
  relation: SUPPORTS_CLAIM
  status: exact_within_selected_reduced_branch
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_5PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_FULL_BUNDLE_TARGET_SURFACE.
- source: FILE_5PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_FULL_BUNDLE_TARGET_SURFACE.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_FULL_BUNDLE_TARGET_SURFACE.
- source: EQ_GROUPED_RESPONSE_MOMENTS
  relation: FEEDS
  status: reduced
  note: Response moments enter one-pole and full-bundle target surface.
- source: EQ_P0_TARGET
  relation: IS_COMPONENT_OF
  status: conditional
  note: Universal quadrupole target is one part of the full-bundle surface.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_reduced_bundle
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Isotropic full-bundle target surface

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_FULL_BUNDLE_TARGET_SURFACE`  
> **Status:** `exact_within_reduced_bundle`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Combined one-pole and quadrupole-normalization target surface.

## Physical Meaning

Combined one-pole and quadrupole-normalization target surface.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_reduced_bundle`
- Parent node: [[MATH_FULL_BUNDLE_TARGET_SURFACE]]

## Equation

$$
D0(B4+Z4)=3(M+B2+Z2)²;  mhat0²N0/D0=54Gc_s⁵/(5a⁵c⁵)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_FULL_BUNDLE_TARGET_SURFACE]]

### Related equations
- [[EQ_GROUPED_RESPONSE_MOMENTS]]
- [[EQ_P0_TARGET]]

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]]
- [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_5PN_FULL]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_FULL_BUNDLE_TARGET_SURFACE]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_5PN_FULL]] | Paper backlink block references EQ_FULL_BUNDLE_TARGET_SURFACE. |
| `CONTAINS_EQUATION` | [[FILE_5PN_FULL]] | Source artifact contains or supports EQ_FULL_BUNDLE_TARGET_SURFACE. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_FULL_BUNDLE_TARGET_SURFACE. |
| `FEEDS` | [[EQ_GROUPED_RESPONSE_MOMENTS]] | Response moments enter one-pole and full-bundle target surface. |
| `IS_COMPONENT_OF` | [[EQ_P0_TARGET]] | Universal quadrupole target is one part of the full-bundle surface. |

## Source Anchors

### Source anchor notes
- [[FILE_5PN_FULL]]
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/5pn/5pn_notes_full.md`
- `notes/pde_audit_full.md`
- `5pn_notes_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
