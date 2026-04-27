---
id: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
title: Localized Maxwell/mixed block carries first honest outgoing bridge
type: claim
layer: claim_theorem
status: exact_within_reduced_mixed_kernel
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The Maxwell/mixed block contributes conservative self-energy and transfers a passive/outgoing l=2 odd fingerprint to the wall through a nonnegative static transfer factor.
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_MIXED_KERNEL]]'
physical_ids:
- PHYS_MIXED_EM_CORE
- PHYS_OUTGOING_QUADRUPOLE_PORT
math_ids:
- MATH_COMPACT_L2_OUTGOING_FINGERPRINT
- MATH_MAXWELL_MIXED_KERNEL
equation_ids:
- EQ_COMPACT_L2_FINGERPRINT
- EQ_MAXWELL_MIXED_TRANSFER
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_STAGE3_BDG_SCHUR
- CLAIM_STAGE5_P0_TARGET
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_PDE_MIXED_KERNEL
outgoing_edges:
- target: CLAIM_STAGE5_P0_TARGET
  relation: FEEDS
  status: active
  note: Claim-level dependency added in v0.4.
- target: MATH_COMPACT_L2_OUTGOING_FINGERPRINT
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_mixed_kernel
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_MAXWELL_MIXED_KERNEL
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_mixed_kernel
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_STAGE4_MAXWELL_MIXED
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_mixed_kernel
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_PDE_MIXED_KERNEL
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mixed kernel audit.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_within_reduced_mixed_kernel
  note: 'Claim class: exact_within_closure'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: ENABLES
  status: active
  note: Claim-level dependency added in v0.4.
- source: QUERY_MIXED_SECTOR
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- source: PHYS_MIXED_EM_CORE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_mixed_kernel
  note: Physical ontology object grounded by this claim.
- source: PHYS_OUTGOING_QUADRUPOLE_PORT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_mixed_kernel
  note: Physical ontology object grounded by this claim.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Source artifact anchors this claim.
- source: CLAIM_STAGE3_BDG_SCHUR
  relation: PREPARES
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_COMPACT_L2_FINGERPRINT
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Equation anchor supports this named claim.
- source: EQ_MAXWELL_MIXED_TRANSFER
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_within_reduced_mixed_kernel
- topic/maxwell
- topic/moving_throat
- topic/quadrupole
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Localized Maxwell/mixed block carries first honest outgoing bridge

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER`  
> **Status:** `exact_within_reduced_mixed_kernel`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

The Maxwell/mixed block contributes conservative self-energy and transfers a passive/outgoing l=2 odd fingerprint to the wall through a nonnegative static transfer factor.

## Claim

The Maxwell/mixed block contributes conservative self-energy and transfers a passive/outgoing l=2 odd fingerprint to the wall through a nonnegative static transfer factor.

## Physical Meaning

The Maxwell/mixed block contributes conservative self-energy and transfers a passive/outgoing l=2 odd fingerprint to the wall through a nonnegative static transfer factor.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_within_reduced_mixed_kernel`
- Outputs: `MT_STAGE4_MAXWELL_MIXED`, `MATH_MAXWELL_MIXED_KERNEL`, `MATH_COMPACT_L2_OUTGOING_FINGERPRINT`

## Atlas Links

### Related physical nodes
- [[PHYS_MIXED_EM_CORE]]
- [[PHYS_OUTGOING_QUADRUPOLE_PORT]]

### Related math nodes
- [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]]
- [[MATH_MAXWELL_MIXED_KERNEL]]

### Related equations
- [[EQ_COMPACT_L2_FINGERPRINT]]
- [[EQ_MAXWELL_MIXED_TRANSFER]]

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_STAGE3_BDG_SCHUR]]
- [[CLAIM_STAGE5_P0_TARGET]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_MIXED_KERNEL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[CLAIM_STAGE5_P0_TARGET]] | Claim-level dependency added in v0.4. |
| `FEEDS_OR_STATUS_OF` | [[MATH_COMPACT_L2_OUTGOING_FINGERPRINT]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_MAXWELL_MIXED_KERNEL]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_STAGE4_MAXWELL_MIXED]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_MIXED_KERNEL]] | Mixed kernel audit. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block references CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `ENABLES` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim-level dependency added in v0.4. |
| `EXPECTS_TARGET` | [[QUERY_MIXED_SECTOR]] | Query validation expected target node. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MIXED_EM_CORE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_OUTGOING_QUADRUPOLE_PORT]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `PREPARES` | [[CLAIM_STAGE3_BDG_SCHUR]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_COMPACT_L2_FINGERPRINT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_MAXWELL_MIXED_TRANSFER]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_MIXED_KERNEL]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
