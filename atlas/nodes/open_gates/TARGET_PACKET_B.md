---
id: TARGET_PACKET_B
title: 'Packet B: orbit-lock / weak-axisymmetric packet'
type: target_packet
layer: open_gate
status: open_branch_realization
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Actual branch must satisfy orbit-lock and weak-axisymmetric prefactor conditions where relevant.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/5pn/5pn_notes_full.md
- pde_audit_full.md
- 5pn_notes_full.md
legacy_sources:
- pde_audit_full.md:V2-26
- 5pn_notes_full.md
physical_ids:
- PHYS_FINITE_MOUTH_SHAPE
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
- CLAIM_PACKET_A_PACKET_B_SPLIT
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
incoming_edges:
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: AFFECTS
  status: reduced/open
  note: P22/orbit-lock/mouth bracing enters weak-axisymmetric packet.
- source: MT_V2_26_STATUS
  relation: DEFINES
  status: open
  note: Defines orbit-lock/weak-axisymmetric branch-output packet.
- source: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: FEEDS_OR_STATUS_OF
  status: open_branch_packets
  note: Claim feeds this downstream object, output, or open gate.
- source: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: MUST_REALIZE
  status: open
  note: Physical branch must output Packet B where relevant.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open_branch_realization
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/target_packet
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Packet B: orbit-lock / weak-axisymmetric packet

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `TARGET_PACKET_B`
> **Status:** `open_branch_realization`
> **Layer:** `open_gate`
> **Type:** `target_packet`

## Summary

Actual branch must satisfy orbit-lock and weak-axisymmetric prefactor conditions where relevant.

## What Remains Open

Actual branch must satisfy orbit-lock and weak-axisymmetric prefactor conditions where relevant.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Actual branch must satisfy orbit-lock and weak-axisymmetric prefactor conditions where relevant.

## Mathematical Role

- Layer: `open_gate`
- Type: `target_packet`
- Status: `open_branch_realization`

## Equation

$$
dln_R_tr=0
$$

$$
dln_R_target=0
$$

$$
dln_epsilon_eta=0
$$

$$
N_Q=1
$$

$$
Xi1=P1/P0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_FINITE_MOUTH_SHAPE]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `AFFECTS` | [[PHYS_FINITE_MOUTH_SHAPE]] | P22/orbit-lock/mouth bracing enters weak-axisymmetric packet. |
| `DEFINES` | [[MT_V2_26_STATUS]] | Defines orbit-lock/weak-axisymmetric branch-output packet. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Claim feeds this downstream object, output, or open gate. |
| `MUST_REALIZE` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Physical branch must output Packet B where relevant. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/5pn/5pn_notes_full.md`
- `pde_audit_full.md`
- `5pn_notes_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
