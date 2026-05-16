---
id: CLAIM_STAGE1_GEOMETRY_LIFT
title: Distributed Sigma/R geometry lift
type: claim
layer: claim_theorem
status: effective_geometry_lift
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The throat is represented as Sigma=r-R(Omega,w,t), so a(t), L(t) become collective moments and scalar/P2/geometry lanes become explicit wall/support modes.
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_MT_PARENT_FIELDS]]'
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
- PHYS_FINITE_MOUTH_SHAPE
- PHYS_GROUPED_P2_SHAPE_RESPONSE
- PHYS_REG_PARTICLE_THROAT
math_ids:
- MATH_SIGMA_R_FIELD
equation_ids:
- EQ_WALL_MODAL_PDE
claim_ids:
- CLAIM_LEPTON_MASS_REDUCED_LEDGER
- CLAIM_PARENT_ACTION_CURRENT_EXACT
- CLAIM_PARENT_WALL_STATUS_SPLIT
- CLAIM_STAGE2_AL_RECOVERY
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_MT_PARENT_FIELDS
outgoing_edges:
- target: MATH_SIGMA_R_FIELD
  relation: FEEDS_OR_STATUS_OF
  status: effective_geometry_lift
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_STAGE1_GEOMETRY_LIFT
  relation: FEEDS_OR_STATUS_OF
  status: effective_geometry_lift
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE2_AL_RECOVERY
  relation: LIFTS_THEN_REDUCES_TO
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_MT_PARENT_FIELDS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent/moving-throat fields.
- source: BACKLINK_2PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE1_GEOMETRY_LIFT.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE1_GEOMETRY_LIFT.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: effective_geometry_lift
  note: 'Claim class: effective_closure'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: GROUNDS_PHYSICAL_MEANING
  status: effective_geometry_lift
  note: Physical ontology object grounded by this claim.
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: GROUNDS_PHYSICAL_MEANING
  status: effective_geometry_lift
  note: Physical ontology object grounded by this claim.
- source: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: GROUNDS_PHYSICAL_MEANING
  status: effective_geometry_lift
  note: Physical ontology object grounded by this claim.
- source: PHYS_REG_PARTICLE_THROAT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: effective_geometry_lift
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: effective_geometry_lift
  note: Source artifact anchors this claim.
- source: CLAIM_LEPTON_MASS_REDUCED_LEDGER
  relation: PHYSICAL_THROAT_CONTEXT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: QUERY_P2_BECOMES_PHYSICAL
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
- source: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: STATUS_LIMITS
  status: active
  note: Claim-level dependency added in v0.4.
- source: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPLIES_PARENT_CONTEXT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_WALL_MODAL_PDE
  relation: SUPPORTS_CLAIM
  status: effective_geometry_lift
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/effective_geometry_lift
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Distributed Sigma/R geometry lift

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_STAGE1_GEOMETRY_LIFT`
> **Status:** `effective_geometry_lift`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The throat is represented as Sigma=r-R(Omega,w,t), so a(t), L(t) become collective moments and scalar/P2/geometry lanes become explicit wall/support modes.

## Claim

The throat is represented as Sigma=r-R(Omega,w,t), so a(t), L(t) become collective moments and scalar/P2/geometry lanes become explicit wall/support modes.

## Physical Meaning

The throat is represented as Sigma=r-R(Omega,w,t), so a(t), L(t) become collective moments and scalar/P2/geometry lanes become explicit wall/support modes.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `effective_geometry_lift`
- Outputs: `MT_STAGE1_GEOMETRY_LIFT`, `MATH_SIGMA_R_FIELD`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]
- [[PHYS_FINITE_MOUTH_SHAPE]]
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]
- [[PHYS_REG_PARTICLE_THROAT]]

### Related math nodes
- [[MATH_SIGMA_R_FIELD]]

### Related equations
- [[EQ_WALL_MODAL_PDE]]

### Related claims
- [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]]
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]
- [[CLAIM_STAGE2_AL_RECOVERY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_MT_PARENT_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_SIGMA_R_FIELD]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_STAGE1_GEOMETRY_LIFT]] | Claim feeds this downstream object, output, or open gate. |
| `LIFTS_THEN_REDUCES_TO` | [[CLAIM_STAGE2_AL_RECOVERY]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_PARENT_FIELDS]] | Parent/moving-throat fields. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_2PN_FULL]] | Paper backlink block references CLAIM_STAGE1_GEOMETRY_LIFT. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block references CLAIM_STAGE1_GEOMETRY_LIFT. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: effective_closure |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_FINITE_MOUTH_SHAPE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Physical ontology object grounded by this claim. |
| `LINKS_TO` | [[PHYS_REG_PARTICLE_THROAT]] | Physical register entry links to graph object. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `PHYSICAL_THROAT_CONTEXT_FOR` | [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]] | Claim-level dependency added in v0.4. |
| `STARTS_AT` | [[QUERY_P2_BECOMES_PHYSICAL]] | Query validation start node. |
| `STATUS_LIMITS` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Claim-level dependency added in v0.4. |
| `SUPPLIES_PARENT_CONTEXT_FOR` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_WALL_MODAL_PDE]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_MT_PARENT_FIELDS]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
