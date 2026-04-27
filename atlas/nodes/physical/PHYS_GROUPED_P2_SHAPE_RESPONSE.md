---
id: PHYS_GROUPED_P2_SHAPE_RESPONSE
title: Grouped real P2 shape/support response
type: response_sector
layer: physical_ontology
status: reduced_controlled
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The first serious finite-mouth/throat quadrupole support-response bundle.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- research/4d_2_5pn/paper/4d_2_5pn.tex
- research/4d_3pn/paper/4d_3pn.tex
- moving_throat_output_full.md
- 4d_2_5pn_summary.md
- 4d_3pn_summary.md
legacy_sources:
- moving_throat_output_full.md
- 4d_2_5pn_summary.md
- 4d_3pn_summary.md
physical_ids:
- PHYS_FINITE_MOUTH_SHAPE
- PHYS_REG_P2_QUAD
math_ids:
- MATH_WALL_MODAL_PDE
equation_ids:
- EQ_WALL_P2_STIFFNESS
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
- CLAIM_2PN_ADM_WITHIN_HIERARCHY
- CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
- CLAIM_STAGE1_GEOMETRY_LIFT
- CLAIM_STAGE2_AL_RECOVERY
- CLAIM_STAGE3_BDG_SCHUR
- CLAIM_STAGE6_FULL_BUNDLE_RATIO
- CLAIM_STAGE7_O3_ISOTROPY
outgoing_edges:
- target: CLAIM_25PN_QUAD_NARROWING
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_theorem_open_normalization
  note: Physical ontology object grounded by this claim.
- target: CLAIM_2PN_ADM_WITHIN_HIERARCHY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_assembly_within_closure
  note: Physical ontology object grounded by this claim.
- target: CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_assembly_within_closure
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: GROUNDS_PHYSICAL_MEANING
  status: effective_geometry_lift
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE2_AL_RECOVERY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_wall_action
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE3_BDG_SCHUR
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_stable_modes
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE6_FULL_BUNDLE_RATIO
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_bundle
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE7_O3_ISOTROPY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_angular_reduced
  note: Physical ontology object grounded by this claim.
incoming_edges:
- source: EQ_WALL_P2_STIFFNESS
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_2PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_GROUPED_P2_SHAPE_RESPONSE.
- source: BACKLINK_3PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_GROUPED_P2_SHAPE_RESPONSE.
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: FEEDS
  status: reduced
  note: Finite-mouth Hessian/tidal loading drives P0/P2 response.
- source: CLAIM_2PN_ADM_WITHIN_HIERARCHY
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_STAGE7_OVERLAP_ISOTROPY
  relation: GIVES_PHYSICAL_DIAGNOSTIC
  status: reduced
  note: O(3) isotropy and weak-axisymmetric b=3a pattern classify shape response.
- source: PHYS_REG_P2_QUAD
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: MATH_WALL_MODAL_PDE
  relation: SUPPORTS
  status: exact_if_closure
  note: l=2 modes are literal wall/support modes.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/reduced_controlled
- topic/pn_chain
- topic/quadrupole
- type/response_sector
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Grouped real P2 shape/support response

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_GROUPED_P2_SHAPE_RESPONSE`  
> **Status:** `reduced_controlled`  
> **Layer:** `physical_ontology`  
> **Type:** `response_sector`

## Summary

The first serious finite-mouth/throat quadrupole support-response bundle.

## Physical Meaning

The first serious finite-mouth/throat quadrupole support-response bundle.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `response_sector`
- Status: `reduced_controlled`

## Equation

$$
20,21,22 grouped lanes
$$

$$
real STF l=2
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_FINITE_MOUTH_SHAPE]]
- [[PHYS_REG_P2_QUAD]]

### Related math nodes
- [[MATH_WALL_MODAL_PDE]]

### Related equations
- [[EQ_WALL_P2_STIFFNESS]]

### Related claims
- [[CLAIM_25PN_QUAD_NARROWING]]
- [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]]
- [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]
- [[CLAIM_STAGE2_AL_RECOVERY]]
- [[CLAIM_STAGE3_BDG_SCHUR]]
- [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]]
- [[CLAIM_STAGE7_O3_ISOTROPY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_25PN_QUAD_NARROWING]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE2_AL_RECOVERY]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE3_BDG_SCHUR]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Physical ontology object grounded by this claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_WALL_P2_STIFFNESS]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_2PN_FULL]] | Paper backlink block references PHYS_GROUPED_P2_SHAPE_RESPONSE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_3PN_FULL]] | Paper backlink block references PHYS_GROUPED_P2_SHAPE_RESPONSE. |
| `FEEDS` | [[PHYS_FINITE_MOUTH_SHAPE]] | Finite-mouth Hessian/tidal loading drives P0/P2 response. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]] | Claim feeds this downstream object, output, or open gate. |
| `GIVES_PHYSICAL_DIAGNOSTIC` | [[MT_STAGE7_OVERLAP_ISOTROPY]] | O(3) isotropy and weak-axisymmetric b=3a pattern classify shape response. |
| `LINKS_TO` | [[PHYS_REG_P2_QUAD]] | Physical register entry links to graph object. |
| `SUPPORTS` | [[MATH_WALL_MODAL_PDE]] | l=2 modes are literal wall/support modes. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/moving_throat_notes_full.md`
- `research/4d_2_5pn/paper/4d_2_5pn.tex`
- `research/4d_3pn/paper/4d_3pn.tex`
- `moving_throat_output_full.md`
- `4d_2_5pn_summary.md`
- `4d_3pn_summary.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
