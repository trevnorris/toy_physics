---
id: CLAIM_STAGE2_AL_RECOVERY
title: Old a,L closure recovered as lowest wall truncation
type: claim
layer: claim_theorem
status: exact_within_wall_action
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The axisymmetric two-mode wall ansatz reduces the distributed wall theory back to an (a,L)-type matrix closure, while P2 is the next harmonic family of the same action.
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_WALL_VARIATION]]'
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
- PHYS_MOUTH_CROSS_SECTION
math_ids:
- MATH_WALL_ACTION_S_ETA
equation_ids:
- EQ_WALL_P2_STIFFNESS
claim_ids:
- CLAIM_STAGE1_GEOMETRY_LIFT
- CLAIM_STAGE3_BDG_SCHUR
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_PDE_WALL_VARIATION
outgoing_edges:
- target: MATH_WALL_ACTION_S_ETA
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_wall_action
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_STAGE2_BREATHING_REDUCTION
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_wall_action
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE3_BDG_SCHUR
  relation: PREPARES
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_PDE_WALL_VARIATION
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: S_eta gives wall PDE.
- source: BACKLINK_2PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE2_AL_RECOVERY.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_within_wall_action
  note: 'Claim class: exact_within_closure'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_wall_action
  note: Physical ontology object grounded by this claim.
- source: PHYS_MOUTH_CROSS_SECTION
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_wall_action
  note: Physical ontology object grounded by this claim.
- source: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: LIFTS_THEN_REDUCES_TO
  status: active
  note: Claim-level dependency added in v0.4.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_wall_action
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_wall_action
  note: Source artifact anchors this claim.
- source: EQ_WALL_P2_STIFFNESS
  relation: SUPPORTS_CLAIM
  status: exact_within_wall_action
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_within_wall_action
- topic/moving_throat
- topic/projection
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Old a,L closure recovered as lowest wall truncation

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_STAGE2_AL_RECOVERY`  
> **Status:** `exact_within_wall_action`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

The axisymmetric two-mode wall ansatz reduces the distributed wall theory back to an (a,L)-type matrix closure, while P2 is the next harmonic family of the same action.

## Claim

The axisymmetric two-mode wall ansatz reduces the distributed wall theory back to an (a,L)-type matrix closure, while P2 is the next harmonic family of the same action.

## Physical Meaning

The axisymmetric two-mode wall ansatz reduces the distributed wall theory back to an (a,L)-type matrix closure, while P2 is the next harmonic family of the same action.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_within_wall_action`
- Outputs: `MT_STAGE2_BREATHING_REDUCTION`, `MATH_WALL_ACTION_S_ETA`

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]
- [[PHYS_MOUTH_CROSS_SECTION]]

### Related math nodes
- [[MATH_WALL_ACTION_S_ETA]]

### Related equations
- [[EQ_WALL_P2_STIFFNESS]]

### Related claims
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]
- [[CLAIM_STAGE3_BDG_SCHUR]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_WALL_VARIATION]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_WALL_ACTION_S_ETA]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_STAGE2_BREATHING_REDUCTION]] | Claim feeds this downstream object, output, or open gate. |
| `PREPARES` | [[CLAIM_STAGE3_BDG_SCHUR]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_WALL_VARIATION]] | S_eta gives wall PDE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_2PN_FULL]] | Paper backlink block references CLAIM_STAGE2_AL_RECOVERY. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MOUTH_CROSS_SECTION]] | Physical ontology object grounded by this claim. |
| `LIFTS_THEN_REDUCES_TO` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Claim-level dependency added in v0.4. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_WALL_P2_STIFFNESS]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_WALL_VARIATION]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
