---
id: CLAIM_STAGE3_BDG_SCHUR
title: BdG-wall Schur complement gives even support ladder
type: claim
layer: claim_theorem
status: exact_within_reduced_stable_modes
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Stable BdG support modes integrate out by Schur complement, lowering static stiffness, increasing inertia, and generating the even low-frequency ladder for scalar and grouped P2...
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_BDG_AUDIT]]'
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
- PHYS_INTERIOR_SUPPORT
math_ids:
- MATH_BDG_SCHUR_COMPLEMENT
equation_ids:
- EQ_BDG_SCHUR_KERNEL
claim_ids:
- CLAIM_STAGE2_AL_RECOVERY
- CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_PDE_BDG_AUDIT
outgoing_edges:
- target: MATH_BDG_SCHUR_COMPLEMENT
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_stable_modes
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_STAGE3_BDG_COUPLING
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_stable_modes
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
  relation: PREPARES
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_PDE_BDG_AUDIT
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: BdG-wall Schur complement and softening audit.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE3_BDG_SCHUR.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_within_reduced_stable_modes
  note: 'Claim class: exact_within_closure'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_stable_modes
  note: Physical ontology object grounded by this claim.
- source: PHYS_INTERIOR_SUPPORT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_stable_modes
  note: Physical ontology object grounded by this claim.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_stable_modes
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_stable_modes
  note: Source artifact anchors this claim.
- source: CLAIM_STAGE2_AL_RECOVERY
  relation: PREPARES
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_BDG_SCHUR_KERNEL
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_stable_modes
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_within_reduced_stable_modes
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# BdG-wall Schur complement gives even support ladder

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_STAGE3_BDG_SCHUR`  
> **Status:** `exact_within_reduced_stable_modes`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

Stable BdG support modes integrate out by Schur complement, lowering static stiffness, increasing inertia, and generating the even low-frequency ladder for scalar and grouped P2 lanes.

## Claim

Stable BdG support modes integrate out by Schur complement, lowering static stiffness, increasing inertia, and generating the even low-frequency ladder for scalar and grouped P2 lanes.

## Physical Meaning

Stable BdG support modes integrate out by Schur complement, lowering static stiffness, increasing inertia, and generating the even low-frequency ladder for scalar and grouped P2 lanes.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_within_reduced_stable_modes`
- Outputs: `MT_STAGE3_BDG_COUPLING`, `MATH_BDG_SCHUR_COMPLEMENT`

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]
- [[PHYS_INTERIOR_SUPPORT]]

### Related math nodes
- [[MATH_BDG_SCHUR_COMPLEMENT]]

### Related equations
- [[EQ_BDG_SCHUR_KERNEL]]

### Related claims
- [[CLAIM_STAGE2_AL_RECOVERY]]
- [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_BDG_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_BDG_SCHUR_COMPLEMENT]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_STAGE3_BDG_COUPLING]] | Claim feeds this downstream object, output, or open gate. |
| `PREPARES` | [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_BDG_AUDIT]] | BdG-wall Schur complement and softening audit. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_STAGE3_BDG_SCHUR. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_INTERIOR_SUPPORT]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `PREPARES` | [[CLAIM_STAGE2_AL_RECOVERY]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_BDG_SCHUR_KERNEL]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_BDG_AUDIT]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
