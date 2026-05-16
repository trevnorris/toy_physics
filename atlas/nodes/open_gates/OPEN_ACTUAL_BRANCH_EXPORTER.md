---
id: OPEN_ACTUAL_BRANCH_EXPORTER
title: 'Open gate: actual moving-throat branch exporter'
type: branch_realization_gate
layer: open_gate
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Need target-blind physical branch that exports required response packets before comparison.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-26
- pde_audit_full.md:V2-25
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[SEC_5PN_SUMMARY_CLAIMS]]'
physical_ids:
- PHYS_REG_BRANCH_EXPORT
- PHYS_TARGET_BLIND_BRANCH_SELECTION
math_ids:
- MATH_FULL_BUNDLE_TARGET_SURFACE
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
open_gate_ids:
- OPEN_EXECUTABLE_BRANCH_SOLVER
- OPEN_MATERIAL_CLOSURE
- OPEN_N2_N4_MOMENT_SHAPE
- OPEN_SOURCE_PORT_NORMALIZATION
- OPEN_WEAK_AXISYM_ORBIT_LOCK
- TARGET_PACKET_A
- TARGET_PACKET_B
status_firewall_ids:
- FIREWALL_ATOM_REDUCED_SECTOR
- FIREWALL_READOUTS_NOT_THROAT
source_ids:
- FILE_MOVING_THROAT_COMPACT
- SEC_5PN_SUMMARY_CLAIMS
outgoing_edges:
- target: OPEN_MATERIAL_CLOSURE
  relation: BLOCKED_BY
  status: open
  note: Material closure gap may block actual moving-throat PDE branch.
- target: OPEN_N2_N4_MOMENT_SHAPE
  relation: DECOMPOSES_INTO
  status: open
  note: Exporter must derive outgoing moment-shape controls.
- target: OPEN_SOURCE_PORT_NORMALIZATION
  relation: DECOMPOSES_INTO
  status: open
  note: Exporter must derive source/port normalization law.
- target: OPEN_WEAK_AXISYM_ORBIT_LOCK
  relation: DECOMPOSES_INTO
  status: open
  note: Exporter must derive weak-axisymmetric tangent/orbit-lock packet.
- target: OPEN_MATERIAL_CLOSURE
  relation: DEPENDS_ON
  status: open
  note: Exporter may require material closure to select branch.
- target: TARGET_PACKET_A
  relation: MUST_REALIZE
  status: open
  note: Physical branch must output Packet A target-blind.
- target: TARGET_PACKET_B
  relation: MUST_REALIZE
  status: open
  note: Physical branch must output Packet B where relevant.
- target: OPEN_EXECUTABLE_BRANCH_SOLVER
  relation: NEEDS
  status: open
  note: Needs branch solver or analytic branch exporter.
incoming_edges:
- source: SEC_5PN_SUMMARY_CLAIMS
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for OPEN_ACTUAL_BRANCH_EXPORTER
- source: WORKFLOW_STEP_4_OPEN_GATES
  relation: CHECKS
  status: active
  note: Major gate for moving-throat extensions.
- source: MT_V2_26_STATUS
  relation: DECLARES_NEXT_ARTIFACT
  status: open
  note: Next required artifact is actual moving-throat branch packet.
- source: FILE_MOVING_THROAT_COMPACT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: MT_V2_23_OPEN_BRANCH_SOLVER
  relation: DOES_NOT_CLOSE
  status: open
  note: First target-blind simulated branch missed target; exporter still incomplete.
- source: QUERY_READOUT_DISCIPLINE
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- source: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: BACKLINK_2PN_FULL
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER.
- source: BACKLINK_3PN_FULL
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER.
- source: BACKLINK_5PN_FULL
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER.
- source: BACKLINK_LEPTON_MASS
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER.
- source: BACKLINK_PDE_AUDIT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER.
- source: ATLAS_CURRENT_READINESS_V05
  relation: FLAGS_REMAINING_GATE
  status: v05
  note: Still open after atlas organization; atlas tracks it but does not solve it.
- source: ATLAS_CURRENT_READINESS_V06
  relation: FLAGS_REMAINING_GATE
  status: v06
  note: Still open after v0.6 organization pass.
- source: MT_V2_16_NO_REFIT
  relation: GOVERNS
  status: mandatory
  note: Branch exporter must freeze protocol before target comparison.
- source: PHYS_TARGET_BLIND_BRANCH_SELECTION
  relation: GOVERNS
  status: mandatory
  note: Exporter must be frozen before comparing to target packets.
- source: PHYS_REG_BRANCH_EXPORT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: MT_V2_20_WEAK_FORM_EXTRACTOR
  relation: PREPARES
  status: prepared
  note: Weak-form extraction protocol is a preparation for actual branch exporter.
- source: FIREWALL_ATOM_REDUCED_SECTOR
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_READOUTS_NOT_THROAT
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: MT_V2_25_ACTUAL_BRANCH_PROTOCOL
  relation: REFINES
  status: open
  note: Defines Packet A/B and no-refit actual branch protocol.
- source: MATH_FULL_BUNDLE_TARGET_SURFACE
  relation: REQUIRES
  status: open
  note: Target surface can only be tested after branch exporter supplies frozen data.
- source: MT_V2_26_STATUS
  relation: SUMMARIZES
  status: open
  note: 'Current status: real audit work yes; completed one-PDE derivation no.'
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open
- topic/moving_throat
- type/branch_realization_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Open gate: actual moving-throat branch exporter

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_ACTUAL_BRANCH_EXPORTER`
> **Status:** `open`
> **Layer:** `open_gate`
> **Type:** `branch_realization_gate`

## Summary

Need target-blind physical branch that exports required response packets before comparison.

## What Remains Open

Need target-blind physical branch that exports required response packets before comparison.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Need target-blind physical branch that exports required response packets before comparison.

## Mathematical Role

- Layer: `open_gate`
- Type: `branch_realization_gate`
- Status: `open`

## Equation

$$
solve branch -> extract K,M,B_n,Z_n,N_n -> freeze -> compare
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_REG_BRANCH_EXPORT]]
- [[PHYS_TARGET_BLIND_BRANCH_SELECTION]]

### Related math nodes
- [[MATH_FULL_BUNDLE_TARGET_SURFACE]]

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]

### Open gates
- [[OPEN_EXECUTABLE_BRANCH_SOLVER]]
- [[OPEN_MATERIAL_CLOSURE]]
- [[OPEN_N2_N4_MOMENT_SHAPE]]
- [[OPEN_SOURCE_PORT_NORMALIZATION]]
- [[OPEN_WEAK_AXISYM_ORBIT_LOCK]]
- [[TARGET_PACKET_A]]
- [[TARGET_PACKET_B]]

### Status firewalls
- [[FIREWALL_ATOM_REDUCED_SECTOR]]
- [[FIREWALL_READOUTS_NOT_THROAT]]

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[SEC_5PN_SUMMARY_CLAIMS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `BLOCKED_BY` | [[OPEN_MATERIAL_CLOSURE]] | Material closure gap may block actual moving-throat PDE branch. |
| `DECOMPOSES_INTO` | [[OPEN_N2_N4_MOMENT_SHAPE]] | Exporter must derive outgoing moment-shape controls. |
| `DECOMPOSES_INTO` | [[OPEN_SOURCE_PORT_NORMALIZATION]] | Exporter must derive source/port normalization law. |
| `DECOMPOSES_INTO` | [[OPEN_WEAK_AXISYM_ORBIT_LOCK]] | Exporter must derive weak-axisymmetric tangent/orbit-lock packet. |
| `DEPENDS_ON` | [[OPEN_MATERIAL_CLOSURE]] | Exporter may require material closure to select branch. |
| `MUST_REALIZE` | [[TARGET_PACKET_A]] | Physical branch must output Packet A target-blind. |
| `MUST_REALIZE` | [[TARGET_PACKET_B]] | Physical branch must output Packet B where relevant. |
| `NEEDS` | [[OPEN_EXECUTABLE_BRANCH_SOLVER]] | Needs branch solver or analytic branch exporter. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_CLAIMS]] | v0.6 section anchor for OPEN_ACTUAL_BRANCH_EXPORTER |
| `CHECKS` | [[WORKFLOW_STEP_4_OPEN_GATES]] | Major gate for moving-throat extensions. |
| `DECLARES_NEXT_ARTIFACT` | [[MT_V2_26_STATUS]] | Next required artifact is actual moving-throat branch packet. |
| `DOCUMENTS` | [[FILE_MOVING_THROAT_COMPACT]] | File anchor documents this node. |
| `DOES_NOT_CLOSE` | [[MT_V2_23_OPEN_BRANCH_SOLVER]] | First target-blind simulated branch missed target; exporter still incomplete. |
| `EXPECTS_TARGET` | [[QUERY_READOUT_DISCIPLINE]] | Query validation expected target node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Claim feeds this downstream object, output, or open gate. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_2PN_FULL]] | Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_3PN_FULL]] | Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_5PN_FULL]] | Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_LEPTON_MASS]] | Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block flags open gate OPEN_ACTUAL_BRANCH_EXPORTER. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V05]] | Still open after atlas organization; atlas tracks it but does not solve it. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V06]] | Still open after v0.6 organization pass. |
| `GOVERNS` | [[MT_V2_16_NO_REFIT]] | Branch exporter must freeze protocol before target comparison. |
| `GOVERNS` | [[PHYS_TARGET_BLIND_BRANCH_SELECTION]] | Exporter must be frozen before comparing to target packets. |
| `LINKS_TO` | [[PHYS_REG_BRANCH_EXPORT]] | Physical register entry links to graph object. |
| `PREPARES` | [[MT_V2_20_WEAK_FORM_EXTRACTOR]] | Weak-form extraction protocol is a preparation for actual branch exporter. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ATOM_REDUCED_SECTOR]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_READOUTS_NOT_THROAT]] | Firewall preserves this correct status boundary. |
| `REFINES` | [[MT_V2_25_ACTUAL_BRANCH_PROTOCOL]] | Defines Packet A/B and no-refit actual branch protocol. |
| `REQUIRES` | [[MATH_FULL_BUNDLE_TARGET_SURFACE]] | Target surface can only be tested after branch exporter supplies frozen data. |
| `SUMMARIZES` | [[MT_V2_26_STATUS]] | Current status: real audit work yes; completed one-PDE derivation no. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[SEC_5PN_SUMMARY_CLAIMS]]

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
