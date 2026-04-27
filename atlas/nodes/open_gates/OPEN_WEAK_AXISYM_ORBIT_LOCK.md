---
id: OPEN_WEAK_AXISYM_ORBIT_LOCK
title: Weak-axisymmetric tangent and orbit-lock packet
type: orbit_lock_gate
layer: open_gate
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Derive weak-axisymmetric tangent, Xi1, q_tr/q_nt/q_eta packet, and orbit-lock conditions on the actual branch.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
source_links:
- '[[SEC_5PN_STAGE11_MONOMIALS]]'
- '[[SEC_5PN_STAGE20_EVEN_GATES]]'
- '[[SEC_5PN_SUMMARY_CLAIMS]]'
- '[[SEC_G2_PREF_SLOPE]]'
math_ids:
- MATH_WEAK_AXISYM_SPLITTING
- MATH_XI1_PREF_SLOPE
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_G2_COMMON_QUOTIENT
- CLAIM_PACKET_A_PACKET_B_SPLIT
- CLAIM_STAGE7_O3_ISOTROPY
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
source_ids:
- SEC_5PN_STAGE11_MONOMIALS
- SEC_5PN_STAGE20_EVEN_GATES
- SEC_5PN_SUMMARY_CLAIMS
- SEC_G2_PREF_SLOPE
incoming_edges:
- source: SEC_5PN_STAGE11_MONOMIALS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Monomial invariants and similarity-orbit closure.
- source: SEC_5PN_STAGE20_EVEN_GATES
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Full even-gate solve with reinstated Z-sector.
- source: SEC_5PN_SUMMARY_CLAIMS
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for OPEN_WEAK_AXISYM_ORBIT_LOCK
- source: SEC_G2_PREF_SLOPE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Xi1/P1/P0 bridge.
- source: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: DECOMPOSES_INTO
  status: open
  note: Exporter must derive weak-axisymmetric tangent/orbit-lock packet.
- source: MATH_WEAK_AXISYM_SPLITTING
  relation: FEEDS
  status: open
  note: Weak-axisymmetric tangent/orbit-lock packet depends on this signature.
- source: MATH_XI1_PREF_SLOPE
  relation: FEEDS
  status: open
  note: Actual branch must export Xi1 or show orbit-lock.
- source: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle_open_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_G2_COMMON_QUOTIENT
  relation: FEEDS_OR_STATUS_OF
  status: conditional_reduced_residual
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: FEEDS_OR_STATUS_OF
  status: open_branch_packets
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_STAGE7_O3_ISOTROPY
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- source: BACKLINK_5PN_FULL
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_WEAK_AXISYM_ORBIT_LOCK.
- source: BACKLINK_G2_OUTPUT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_WEAK_AXISYM_ORBIT_LOCK.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_WEAK_AXISYM_ORBIT_LOCK.
- source: ATLAS_CURRENT_READINESS_V05
  relation: FLAGS_REMAINING_GATE
  status: v05
  note: Still open after atlas organization; atlas tracks it but does not solve it.
- source: ATLAS_CURRENT_READINESS_V06
  relation: FLAGS_REMAINING_GATE
  status: v06
  note: Still open after v0.6 organization pass.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open
- topic/moving_throat
- type/orbit_lock_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weak-axisymmetric tangent and orbit-lock packet

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_WEAK_AXISYM_ORBIT_LOCK`  
> **Status:** `open`  
> **Layer:** `open_gate`  
> **Type:** `orbit_lock_gate`

## Summary

Derive weak-axisymmetric tangent, Xi1, q_tr/q_nt/q_eta packet, and orbit-lock conditions on the actual branch.

## What Remains Open

Derive weak-axisymmetric tangent, Xi1, q_tr/q_nt/q_eta packet, and orbit-lock conditions on the actual branch.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Derive weak-axisymmetric tangent, Xi1, q_tr/q_nt/q_eta packet, and orbit-lock conditions on the actual branch.

## Mathematical Role

- Layer: `open_gate`
- Type: `orbit_lock_gate`
- Status: `open`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WEAK_AXISYM_SPLITTING]]
- [[MATH_XI1_PREF_SLOPE]]

### Related equations
- none

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_G2_COMMON_QUOTIENT]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]
- [[CLAIM_STAGE7_O3_ISOTROPY]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- [[SEC_5PN_STAGE11_MONOMIALS]]
- [[SEC_5PN_STAGE20_EVEN_GATES]]
- [[SEC_5PN_SUMMARY_CLAIMS]]
- [[SEC_G2_PREF_SLOPE]]

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_STAGE11_MONOMIALS]] | Monomial invariants and similarity-orbit closure. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_STAGE20_EVEN_GATES]] | Full even-gate solve with reinstated Z-sector. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_CLAIMS]] | v0.6 section anchor for OPEN_WEAK_AXISYM_ORBIT_LOCK |
| `ANCHORS_CLAIM_SECTION` | [[SEC_G2_PREF_SLOPE]] | Xi1/P1/P0 bridge. |
| `DECOMPOSES_INTO` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Exporter must derive weak-axisymmetric tangent/orbit-lock packet. |
| `FEEDS` | [[MATH_WEAK_AXISYM_SPLITTING]] | Weak-axisymmetric tangent/orbit-lock packet depends on this signature. |
| `FEEDS` | [[MATH_XI1_PREF_SLOPE]] | Actual branch must export Xi1 or show orbit-lock. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_G2_COMMON_QUOTIENT]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Claim feeds this downstream object, output, or open gate. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_5PN_FULL]] | Paper backlink block flags open gate OPEN_WEAK_AXISYM_ORBIT_LOCK. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_G2_OUTPUT]] | Paper backlink block flags open gate OPEN_WEAK_AXISYM_ORBIT_LOCK. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block flags open gate OPEN_WEAK_AXISYM_ORBIT_LOCK. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V05]] | Still open after atlas organization; atlas tracks it but does not solve it. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V06]] | Still open after v0.6 organization pass. |

## Source Anchors

### Source anchor notes
- [[SEC_5PN_STAGE11_MONOMIALS]]
- [[SEC_5PN_STAGE20_EVEN_GATES]]
- [[SEC_5PN_SUMMARY_CLAIMS]]
- [[SEC_G2_PREF_SLOPE]]

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_pde_program_compact.md`
- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
