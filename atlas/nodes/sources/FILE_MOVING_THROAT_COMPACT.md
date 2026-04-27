---
id: FILE_MOVING_THROAT_COMPACT
title: moving_throat_pde_program_compact.md
type: source_file
layer: file_anchor
status: master_anchor
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Compact master ledger for moving-throat PDE program and current active bottlenecks.
source_kind: future_paper_note
future_paper_needed: true
source_files:
- notes/moving_throat_pde_program_compact.md
- moving_throat_pde_program_compact.md
legacy_sources:
- moving_throat_pde_program_compact.md
source_links:
- '[[SEC_MT_MIXED_INVARIANTS]]'
- '[[SEC_MT_PARENT_FIELDS]]'
- '[[SEC_MT_PROJECTION_ZERO_MODE]]'
- '[[SEC_MT_READING_RULES]]'
- '[[SEC_MT_THEOREM_STATUS]]'
equation_ids:
- EQ_BDG_SCHUR_KERNEL
- EQ_BULK_CONTINUITY
- EQ_COMPACT_L2_FINGERPRINT
- EQ_GNLS_BULK
- EQ_GROUPED_RESPONSE_MOMENTS
- EQ_MAXWELL_MIXED_TRANSFER
- EQ_P0_TARGET
- EQ_WALL_MODAL_PDE
- EQ_WALL_P2_STIFFNESS
- EQ_WEAK_AXISYM_SIGNATURE
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_BRANCH_EXPORTER_REQUIRED
- CLAIM_MATERIAL_CLOSURE_GAP
- CLAIM_MIXED_RECIRCULATION_OPEN
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PACKET_A_PACKET_B_SPLIT
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE1_GEOMETRY_LIFT
- CLAIM_STAGE2_AL_RECOVERY
- CLAIM_STAGE3_BDG_SCHUR
- CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
- CLAIM_STAGE5_P0_TARGET
- CLAIM_STAGE6_FULL_BUNDLE_RATIO
- CLAIM_STAGE7_O3_ISOTROPY
- CLAIM_STAGE8_TO_14_SELECTED_BRANCH
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
source_ids:
- SEC_MT_MIXED_INVARIANTS
- SEC_MT_PARENT_FIELDS
- SEC_MT_PROJECTION_ZERO_MODE
- SEC_MT_READING_RULES
- SEC_MT_THEOREM_STATUS
outgoing_edges:
- target: EQ_BDG_SCHUR_KERNEL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_BDG_SCHUR_KERNEL.
- target: EQ_BULK_CONTINUITY
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_BULK_CONTINUITY.
- target: EQ_COMPACT_L2_FINGERPRINT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_COMPACT_L2_FINGERPRINT.
- target: EQ_GNLS_BULK
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_GNLS_BULK.
- target: EQ_GROUPED_RESPONSE_MOMENTS
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_GROUPED_RESPONSE_MOMENTS.
- target: EQ_MAXWELL_MIXED_TRANSFER
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MAXWELL_MIXED_TRANSFER.
- target: EQ_P0_TARGET
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_P0_TARGET.
- target: EQ_WALL_MODAL_PDE
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WALL_MODAL_PDE.
- target: EQ_WALL_P2_STIFFNESS
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WALL_P2_STIFFNESS.
- target: EQ_WEAK_AXISYM_SIGNATURE
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE.
- target: MT_STAGE1_GEOMETRY_LIFT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- target: BACKLINK_MOVING_THROAT_COMPACT
  relation: HAS_PAPER_BACKLINK_BLOCK
  status: v06
  note: Paste-ready atlas backlink block prepared in v0.6.
- target: SEC_MT_MIXED_INVARIANTS
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: moving_throat_pde_program_compact.md:483 — 2.5.6 Exact mixed-sector gauge invariants
- target: SEC_MT_PARENT_FIELDS
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: moving_throat_pde_program_compact.md:180 — 2. Parent Theory and Exact Bulk Equations
- target: SEC_MT_PROJECTION_ZERO_MODE
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: moving_throat_pde_program_compact.md:507 — 2.5.7 Cold-start projection and zero-mode reduction hooks
- target: SEC_MT_READING_RULES
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: moving_throat_pde_program_compact.md:49 — 1. Reading Rules and Status Firewall
- target: SEC_MT_THEOREM_STATUS
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: moving_throat_pde_program_compact.md:84 — 1.3 Present theorem-status summary
- target: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Source artifact anchors this claim.
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open_actual_branch
  note: Source artifact anchors this claim.
- target: CLAIM_MATERIAL_CLOSURE_GAP
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open
  note: Source artifact anchors this claim.
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open
  note: Source artifact anchors this claim.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Source artifact anchors this claim.
- target: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open_branch_packets
  note: Source artifact anchors this claim.
- target: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: OWNS_OR_ANCHORS_CLAIM
  status: paper_facing_ontology_discipline
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: effective_geometry_lift
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE2_AL_RECOVERY
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_wall_action
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE3_BDG_SCHUR
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_stable_modes
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE5_P0_TARGET
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_grouped_bridge
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE6_FULL_BUNDLE_RATIO
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_bundle
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE7_O3_ISOTROPY
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_angular_reduced
  note: Source artifact anchors this claim.
- target: CLAIM_STAGE8_TO_14_SELECTED_BRANCH
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_selected_reduced_branch
  note: Source artifact anchors this claim.
tags:
- atlas/node
- atlas/sources
- layer/file_anchor
- status/master_anchor
- topic/moving_throat
- type/source_file
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# moving_throat_pde_program_compact.md

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `FILE_MOVING_THROAT_COMPACT`  
> **Status:** `master_anchor`  
> **Layer:** `file_anchor`  
> **Type:** `source_file`

## Summary

Compact master ledger for moving-throat PDE program and current active bottlenecks.

> [!note] Future paper needed
> This node is intentionally anchored to a notes/derivation source until a maintained paper exists.

## Physical Meaning

Compact master ledger for moving-throat PDE program and current active bottlenecks.

## Mathematical Role

- Layer: `file_anchor`
- Type: `source_file`
- Status: `master_anchor`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_BDG_SCHUR_KERNEL]]
- [[EQ_BULK_CONTINUITY]]
- [[EQ_COMPACT_L2_FINGERPRINT]]
- [[EQ_GNLS_BULK]]
- [[EQ_GROUPED_RESPONSE_MOMENTS]]
- [[EQ_MAXWELL_MIXED_TRANSFER]]
- [[EQ_P0_TARGET]]
- [[EQ_WALL_MODAL_PDE]]
- [[EQ_WALL_P2_STIFFNESS]]
- [[EQ_WEAK_AXISYM_SIGNATURE]]

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]
- [[CLAIM_MATERIAL_CLOSURE_GAP]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]
- [[CLAIM_STAGE2_AL_RECOVERY]]
- [[CLAIM_STAGE3_BDG_SCHUR]]
- [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]]
- [[CLAIM_STAGE5_P0_TARGET]]
- [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]]
- [[CLAIM_STAGE7_O3_ISOTROPY]]
- [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- [[SEC_MT_MIXED_INVARIANTS]]
- [[SEC_MT_PARENT_FIELDS]]
- [[SEC_MT_PROJECTION_ZERO_MODE]]
- [[SEC_MT_READING_RULES]]
- [[SEC_MT_THEOREM_STATUS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[EQ_BDG_SCHUR_KERNEL]] | Source artifact contains or supports EQ_BDG_SCHUR_KERNEL. |
| `CONTAINS_EQUATION` | [[EQ_BULK_CONTINUITY]] | Source artifact contains or supports EQ_BULK_CONTINUITY. |
| `CONTAINS_EQUATION` | [[EQ_COMPACT_L2_FINGERPRINT]] | Source artifact contains or supports EQ_COMPACT_L2_FINGERPRINT. |
| `CONTAINS_EQUATION` | [[EQ_GNLS_BULK]] | Source artifact contains or supports EQ_GNLS_BULK. |
| `CONTAINS_EQUATION` | [[EQ_GROUPED_RESPONSE_MOMENTS]] | Source artifact contains or supports EQ_GROUPED_RESPONSE_MOMENTS. |
| `CONTAINS_EQUATION` | [[EQ_MAXWELL_MIXED_TRANSFER]] | Source artifact contains or supports EQ_MAXWELL_MIXED_TRANSFER. |
| `CONTAINS_EQUATION` | [[EQ_P0_TARGET]] | Source artifact contains or supports EQ_P0_TARGET. |
| `CONTAINS_EQUATION` | [[EQ_WALL_MODAL_PDE]] | Source artifact contains or supports EQ_WALL_MODAL_PDE. |
| `CONTAINS_EQUATION` | [[EQ_WALL_P2_STIFFNESS]] | Source artifact contains or supports EQ_WALL_P2_STIFFNESS. |
| `CONTAINS_EQUATION` | [[EQ_WEAK_AXISYM_SIGNATURE]] | Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE. |
| `DOCUMENTS` | [[MT_STAGE1_GEOMETRY_LIFT]] | File anchor documents this node. |
| `DOCUMENTS` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | File anchor documents this node. |
| `HAS_PAPER_BACKLINK_BLOCK` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paste-ready atlas backlink block prepared in v0.6. |
| `HAS_SECTION_ANCHOR` | [[SEC_MT_MIXED_INVARIANTS]] | moving_throat_pde_program_compact.md:483 — 2.5.6 Exact mixed-sector gauge invariants |
| `HAS_SECTION_ANCHOR` | [[SEC_MT_PARENT_FIELDS]] | moving_throat_pde_program_compact.md:180 — 2. Parent Theory and Exact Bulk Equations |
| `HAS_SECTION_ANCHOR` | [[SEC_MT_PROJECTION_ZERO_MODE]] | moving_throat_pde_program_compact.md:507 — 2.5.7 Cold-start projection and zero-mode reduction hooks |
| `HAS_SECTION_ANCHOR` | [[SEC_MT_READING_RULES]] | moving_throat_pde_program_compact.md:49 — 1. Reading Rules and Status Firewall |
| `HAS_SECTION_ANCHOR` | [[SEC_MT_THEOREM_STATUS]] | moving_throat_pde_program_compact.md:84 — 1.3 Present theorem-status summary |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_MATERIAL_CLOSURE_GAP]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE2_AL_RECOVERY]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE3_BDG_SCHUR]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE5_P0_TARGET]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]] | Source artifact anchors this claim. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- [[SEC_MT_MIXED_INVARIANTS]]
- [[SEC_MT_PARENT_FIELDS]]
- [[SEC_MT_PROJECTION_ZERO_MODE]]
- [[SEC_MT_READING_RULES]]
- [[SEC_MT_THEOREM_STATUS]]

### Source files
- `notes/moving_throat_pde_program_compact.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
