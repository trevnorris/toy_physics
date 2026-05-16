---
id: CLAIM_BRANCH_EXPORTER_REQUIRED
title: Actual target-blind branch exporter required
type: claim
layer: claim_theorem
status: open_actual_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The program must produce target-blind branch data before comparison; executable fixtures and smoke tests do not close the theorem without a physical branch exporter.
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_5PN_STAGE24_26_PLACEMENT]]'
- '[[SEC_5PN_SUMMARY_BOTTOM_LINE]]'
- '[[SEC_5PN_SUMMARY_SCOPE]]'
- '[[SEC_MT_READING_RULES]]'
- '[[SEC_MT_THEOREM_STATUS]]'
- '[[SEC_PDE_ACTUAL_BRANCH_PACKETS]]'
- '[[SEC_PDE_BRANCH_FIXTURE]]'
- '[[SEC_PDE_BRANCH_FREEZE]]'
- '[[SEC_PDE_PHYSICAL_ONTOLOGY]]'
- '[[SEC_PDE_WEAK_FORM_EXPORTER]]'
physical_ids:
- PHYS_BRANCH_EXPORTER
- PHYS_REG_BRANCH_EXPORT
- PHYS_TARGET_BLIND_BRANCH_SELECTION
equation_ids:
- EQ_FULL_BUNDLE_TARGET_SURFACE
claim_ids:
- CLAIM_MATERIAL_CLOSURE_GAP
- CLAIM_PACKET_A_PACKET_B_SPLIT
- CLAIM_PARENT_WALL_STATUS_SPLIT
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
- OPEN_EXECUTABLE_BRANCH_SOLVER
- TARGET_PACKET_A
- TARGET_PACKET_B
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_5PN_STAGE24_26_PLACEMENT
- SEC_5PN_SUMMARY_BOTTOM_LINE
- SEC_5PN_SUMMARY_SCOPE
- SEC_MT_READING_RULES
- SEC_MT_THEOREM_STATUS
- SEC_PDE_ACTUAL_BRANCH_PACKETS
- SEC_PDE_BRANCH_FIXTURE
- SEC_PDE_BRANCH_FREEZE
- SEC_PDE_PHYSICAL_ONTOLOGY
- SEC_PDE_WEAK_FORM_EXPORTER
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
- target: OPEN_EXECUTABLE_BRANCH_SOLVER
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
- target: TARGET_PACKET_A
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
- target: TARGET_PACKET_B
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_5PN_STAGE24_26_PLACEMENT
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Continuum-kernel extraction and placement map.
- source: SEC_5PN_SUMMARY_BOTTOM_LINE
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_BRANCH_EXPORTER_REQUIRED
- source: SEC_5PN_SUMMARY_SCOPE
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_BRANCH_EXPORTER_REQUIRED
- source: SEC_MT_READING_RULES
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Status firewall.
- source: SEC_MT_THEOREM_STATUS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Current bottleneck summary.
- source: SEC_PDE_ACTUAL_BRANCH_PACKETS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Packet A/B actual branch protocol.
- source: SEC_PDE_BRANCH_FIXTURE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Observable packet fixture and extraction formulas.
- source: SEC_PDE_BRANCH_FREEZE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: No-refit branch freeze protocol.
- source: SEC_PDE_PHYSICAL_ONTOLOGY
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Physical picture and ontology checklist.
- source: SEC_PDE_WEAK_FORM_EXPORTER
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Weak-form branch extraction preparation.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_BRANCH_EXPORTER_REQUIRED.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_BRANCH_EXPORTER_REQUIRED.
- source: CLAIM_MATERIAL_CLOSURE_GAP
  relation: BLOCKS_FULL_BRANCH_REALIZATION_IF_UNSOLVED
  status: active
  note: Claim-level dependency added in v0.4.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: open_actual_branch
  note: 'Claim class: open_gate'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_BRANCH_EXPORTER
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_actual_branch
  note: Physical ontology object grounded by this claim.
- source: PHYS_TARGET_BLIND_BRANCH_SELECTION
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_actual_branch
  note: Physical ontology object grounded by this claim.
- source: PHYS_REG_BRANCH_EXPORT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: OPENS_REQUIRED_WORK
  status: active
  note: Claim-level dependency added in v0.4.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open_actual_branch
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open_actual_branch
  note: Source artifact anchors this claim.
- source: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: REQUIRES
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_FULL_BUNDLE_TARGET_SURFACE
  relation: SUPPORTS_CLAIM
  status: open_actual_branch
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/open_actual_branch
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Actual target-blind branch exporter required

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_BRANCH_EXPORTER_REQUIRED`
> **Status:** `open_actual_branch`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The program must produce target-blind branch data before comparison; executable fixtures and smoke tests do not close the theorem without a physical branch exporter.

## Claim

The program must produce target-blind branch data before comparison; executable fixtures and smoke tests do not close the theorem without a physical branch exporter.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

The program must produce target-blind branch data before comparison; executable fixtures and smoke tests do not close the theorem without a physical branch exporter.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `open_actual_branch`
- Outputs: `OPEN_ACTUAL_BRANCH_EXPORTER`, `OPEN_EXECUTABLE_BRANCH_SOLVER`, `TARGET_PACKET_A`, `TARGET_PACKET_B`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANCH_EXPORTER]]
- [[PHYS_REG_BRANCH_EXPORT]]
- [[PHYS_TARGET_BLIND_BRANCH_SELECTION]]

### Related math nodes
- none

### Related equations
- [[EQ_FULL_BUNDLE_TARGET_SURFACE]]

### Related claims
- [[CLAIM_MATERIAL_CLOSURE_GAP]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]
- [[OPEN_EXECUTABLE_BRANCH_SOLVER]]
- [[TARGET_PACKET_A]]
- [[TARGET_PACKET_B]]

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_5PN_STAGE24_26_PLACEMENT]]
- [[SEC_5PN_SUMMARY_BOTTOM_LINE]]
- [[SEC_5PN_SUMMARY_SCOPE]]
- [[SEC_MT_READING_RULES]]
- [[SEC_MT_THEOREM_STATUS]]
- [[SEC_PDE_ACTUAL_BRANCH_PACKETS]]
- [[SEC_PDE_BRANCH_FIXTURE]]
- [[SEC_PDE_BRANCH_FREEZE]]
- [[SEC_PDE_PHYSICAL_ONTOLOGY]]
- [[SEC_PDE_WEAK_FORM_EXPORTER]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_EXECUTABLE_BRANCH_SOLVER]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[TARGET_PACKET_A]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[TARGET_PACKET_B]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_STAGE24_26_PLACEMENT]] | Continuum-kernel extraction and placement map. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_BOTTOM_LINE]] | v0.6 section anchor for CLAIM_BRANCH_EXPORTER_REQUIRED |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_SCOPE]] | v0.6 section anchor for CLAIM_BRANCH_EXPORTER_REQUIRED |
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_READING_RULES]] | Status firewall. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_THEOREM_STATUS]] | Current bottleneck summary. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_ACTUAL_BRANCH_PACKETS]] | Packet A/B actual branch protocol. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_BRANCH_FIXTURE]] | Observable packet fixture and extraction formulas. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_BRANCH_FREEZE]] | No-refit branch freeze protocol. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_PHYSICAL_ONTOLOGY]] | Physical picture and ontology checklist. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_WEAK_FORM_EXPORTER]] | Weak-form branch extraction preparation. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block references CLAIM_BRANCH_EXPORTER_REQUIRED. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_BRANCH_EXPORTER_REQUIRED. |
| `BLOCKS_FULL_BRANCH_REALIZATION_IF_UNSOLVED` | [[CLAIM_MATERIAL_CLOSURE_GAP]] | Claim-level dependency added in v0.4. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: open_gate |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_BRANCH_EXPORTER]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_TARGET_BLIND_BRANCH_SELECTION]] | Physical ontology object grounded by this claim. |
| `LINKS_TO` | [[PHYS_REG_BRANCH_EXPORT]] | Physical register entry links to graph object. |
| `OPENS_REQUIRED_WORK` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Claim-level dependency added in v0.4. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `REQUIRES` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_FULL_BUNDLE_TARGET_SURFACE]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_5PN_STAGE24_26_PLACEMENT]]
- [[SEC_5PN_SUMMARY_BOTTOM_LINE]]
- [[SEC_5PN_SUMMARY_SCOPE]]
- [[SEC_MT_READING_RULES]]
- [[SEC_MT_THEOREM_STATUS]]
- [[SEC_PDE_ACTUAL_BRANCH_PACKETS]]
- [[SEC_PDE_BRANCH_FIXTURE]]
- [[SEC_PDE_BRANCH_FREEZE]]
- [[SEC_PDE_PHYSICAL_ONTOLOGY]]
- [[SEC_PDE_WEAK_FORM_EXPORTER]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
