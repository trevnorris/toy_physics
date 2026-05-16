---
id: MT_V2_10_HAMILTONIAN_STABILITY
title: V2-10 Hamiltonian/stability audit
type: audit_gate
layer: status_audit
status: mandatory_gates
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Defines Schur-complement stability, outgoing passivity, and failure modes for the one-lane conservative bundle.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
math_ids:
- MATH_MAXWELL_MIXED_KERNEL
outgoing_edges:
- target: MATH_MAXWELL_MIXED_KERNEL
  relation: GATES
  status: mandatory
  note: Adds internal-block, wall-softening, ghost/Krein, and dark-port failure modes.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/mandatory_gates
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-10 Hamiltonian/stability audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_10_HAMILTONIAN_STABILITY`
> **Status:** `mandatory_gates`
> **Layer:** `status_audit`
> **Type:** `audit_gate`

## Summary

Defines Schur-complement stability, outgoing passivity, and failure modes for the one-lane conservative bundle.

## Physical Meaning

Defines Schur-complement stability, outgoing passivity, and failure modes for the one-lane conservative bundle.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `mandatory_gates`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MAXWELL_MIXED_KERNEL]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GATES` | [[MATH_MAXWELL_MIXED_KERNEL]] | Adds internal-block, wall-softening, ghost/Krein, and dark-port failure modes. |

## Incoming Edges

- none

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
