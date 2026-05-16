---
id: MT_V2_26_STATUS
title: V2-26 post-audit program status
type: status_summary
layer: status_audit
status: boundary_set
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Reduced/effective framework is internally consistent and executable, but actual branch realization is unsolved.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-26
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
- TARGET_PACKET_A
- TARGET_PACKET_B
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: DECLARES_NEXT_ARTIFACT
  status: open
  note: Next required artifact is actual moving-throat branch packet.
- target: TARGET_PACKET_A
  relation: DEFINES
  status: open
  note: Defines conservative/outgoing branch-output packet.
- target: TARGET_PACKET_B
  relation: DEFINES
  status: open
  note: Defines orbit-lock/weak-axisymmetric branch-output packet.
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: SUMMARIZES
  status: open
  note: 'Current status: real audit work yes; completed one-PDE derivation no.'
incoming_edges:
- source: MT_V2_23_OPEN_BRANCH_SOLVER
  relation: FEEDS_STATUS
  status: audit
  note: Target-blind misses inform honest status boundary.
- source: MT_V2_27_EXECUTABLE_AUDIT
  relation: SUPPORTS
  status: reproducibility
  note: Executable harness supports audit claims but does not prove branch existence.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/boundary_set
- topic/moving_throat
- topic/quadrupole
- type/status_summary
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-26 post-audit program status

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_26_STATUS`
> **Status:** `boundary_set`
> **Layer:** `status_audit`
> **Type:** `status_summary`

## Summary

Reduced/effective framework is internally consistent and executable, but actual branch realization is unsolved.

## Physical Meaning

Reduced/effective framework is internally consistent and executable, but actual branch realization is unsolved.

## Mathematical Role

- Layer: `status_audit`
- Type: `status_summary`
- Status: `boundary_set`

## Equation

$$
D0*C/(3A²)=1
$$

$$
P0=N0/D0
$$

$$
P2=0
$$

$$
P4=0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- none

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]
- [[TARGET_PACKET_A]]
- [[TARGET_PACKET_B]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DECLARES_NEXT_ARTIFACT` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Next required artifact is actual moving-throat branch packet. |
| `DEFINES` | [[TARGET_PACKET_A]] | Defines conservative/outgoing branch-output packet. |
| `DEFINES` | [[TARGET_PACKET_B]] | Defines orbit-lock/weak-axisymmetric branch-output packet. |
| `SUMMARIZES` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Current status: real audit work yes; completed one-PDE derivation no. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_STATUS` | [[MT_V2_23_OPEN_BRANCH_SOLVER]] | Target-blind misses inform honest status boundary. |
| `SUPPORTS` | [[MT_V2_27_EXECUTABLE_AUDIT]] | Executable harness supports audit claims but does not prove branch existence. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
