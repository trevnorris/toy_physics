---
id: CLAIM_PACKET_A_PACKET_B_SPLIT
title: Packet A / Packet B target split
type: claim
layer: claim_theorem
status: open_branch_packets
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Current moving-throat finish line separates conservative/outgoing grouped target data (Packet A) from orbit-lock/weak-axisymmetric data (Packet B), both requiring actual branch ...
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_MT_THEOREM_STATUS]]'
- '[[SEC_PDE_ACTUAL_BRANCH_PACKETS]]'
- '[[SEC_PDE_BRANCH_FREEZE]]'
physical_ids:
- PHYS_FINITE_MOUTH_SHAPE
- PHYS_RESPONSE_READOUTS
equation_ids:
- EQ_FULL_BUNDLE_TARGET_SURFACE
- EQ_XI1_PREF_SLOPE
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
- CLAIM_G2_COMMON_QUOTIENT
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE023_FULL_BUNDLE_RATIO
open_gate_ids:
- OPEN_WEAK_AXISYM_ORBIT_LOCK
- TARGET_PACKET_A
- TARGET_PACKET_B
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_MT_THEOREM_STATUS
- SEC_PDE_ACTUAL_BRANCH_PACKETS
- SEC_PDE_BRANCH_FREEZE
outgoing_edges:
- target: OPEN_WEAK_AXISYM_ORBIT_LOCK
  relation: FEEDS_OR_STATUS_OF
  status: open_branch_packets
  note: Claim feeds this downstream object, output, or open gate.
- target: TARGET_PACKET_A
  relation: FEEDS_OR_STATUS_OF
  status: open_branch_packets
  note: Claim feeds this downstream object, output, or open gate.
- target: TARGET_PACKET_B
  relation: FEEDS_OR_STATUS_OF
  status: open_branch_packets
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: REQUIRES
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_MT_THEOREM_STATUS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Current bottleneck summary.
- source: SEC_PDE_ACTUAL_BRANCH_PACKETS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Packet A/B actual branch protocol.
- source: SEC_PDE_BRANCH_FREEZE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: No-refit branch freeze protocol.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_PACKET_A_PACKET_B_SPLIT.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: open_branch_packets
  note: 'Claim class: open_gate'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: CLAIM_STAGE023_FULL_BUNDLE_RATIO
  relation: DEFINES_PACKET_A
  status: active
  note: Claim-level dependency added in v0.4.
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_branch_packets
  note: Physical ontology object grounded by this claim.
- source: PHYS_RESPONSE_READOUTS
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_branch_packets
  note: Physical ontology object grounded by this claim.
- source: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: ORGANIZES
  status: active
  note: Claim-level dependency added in v0.4.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open_branch_packets
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open_branch_packets
  note: Source artifact anchors this claim.
- source: EQ_FULL_BUNDLE_TARGET_SURFACE
  relation: SUPPORTS_CLAIM
  status: open_branch_packets
  note: Equation anchor supports this named claim.
- source: EQ_XI1_PREF_SLOPE
  relation: SUPPORTS_CLAIM
  status: open_branch_packets
  note: Equation anchor supports this named claim.
- source: CLAIM_G2_COMMON_QUOTIENT
  relation: USES_PACKET_B
  status: active
  note: Claim-level dependency added in v0.4.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/open_branch_packets
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Packet A / Packet B target split

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_PACKET_A_PACKET_B_SPLIT`
> **Status:** `open_branch_packets`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

Current moving-throat finish line separates conservative/outgoing grouped target data (Packet A) from orbit-lock/weak-axisymmetric data (Packet B), both requiring actual branch realization.

## Claim

Current moving-throat finish line separates conservative/outgoing grouped target data (Packet A) from orbit-lock/weak-axisymmetric data (Packet B), both requiring actual branch realization.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

Current moving-throat finish line separates conservative/outgoing grouped target data (Packet A) from orbit-lock/weak-axisymmetric data (Packet B), both requiring actual branch realization.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `open_branch_packets`
- Outputs: `TARGET_PACKET_A`, `TARGET_PACKET_B`, `OPEN_WEAK_AXISYM_ORBIT_LOCK`

## Atlas Links

### Related physical nodes
- [[PHYS_FINITE_MOUTH_SHAPE]]
- [[PHYS_RESPONSE_READOUTS]]

### Related math nodes
- none

### Related equations
- [[EQ_FULL_BUNDLE_TARGET_SURFACE]]
- [[EQ_XI1_PREF_SLOPE]]

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]
- [[CLAIM_G2_COMMON_QUOTIENT]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE023_FULL_BUNDLE_RATIO]]

### Open gates
- [[OPEN_WEAK_AXISYM_ORBIT_LOCK]]
- [[TARGET_PACKET_A]]
- [[TARGET_PACKET_B]]

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_MT_THEOREM_STATUS]]
- [[SEC_PDE_ACTUAL_BRANCH_PACKETS]]
- [[SEC_PDE_BRANCH_FREEZE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[OPEN_WEAK_AXISYM_ORBIT_LOCK]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[TARGET_PACKET_A]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[TARGET_PACKET_B]] | Claim feeds this downstream object, output, or open gate. |
| `REQUIRES` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_THEOREM_STATUS]] | Current bottleneck summary. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_ACTUAL_BRANCH_PACKETS]] | Packet A/B actual branch protocol. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_BRANCH_FREEZE]] | No-refit branch freeze protocol. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block references CLAIM_PACKET_A_PACKET_B_SPLIT. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: open_gate |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `DEFINES_PACKET_A` | [[CLAIM_STAGE023_FULL_BUNDLE_RATIO]] | Claim-level dependency added in v0.4. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_FINITE_MOUTH_SHAPE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_RESPONSE_READOUTS]] | Physical ontology object grounded by this claim. |
| `ORGANIZES` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Claim-level dependency added in v0.4. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_FULL_BUNDLE_TARGET_SURFACE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_XI1_PREF_SLOPE]] | Equation anchor supports this named claim. |
| `USES_PACKET_B` | [[CLAIM_G2_COMMON_QUOTIENT]] | Claim-level dependency added in v0.4. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_MT_THEOREM_STATUS]]
- [[SEC_PDE_ACTUAL_BRANCH_PACKETS]]
- [[SEC_PDE_BRANCH_FREEZE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
