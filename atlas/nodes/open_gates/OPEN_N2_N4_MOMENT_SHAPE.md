---
id: OPEN_N2_N4_MOMENT_SHAPE
title: Outgoing moment-shape control
type: moment_shape_gate
layer: open_gate
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Derive independent outgoing N2 and N4 moment-shape controls on the same frozen branch.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
incoming_edges:
- source: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: DECOMPOSES_INTO
  status: open
  note: Exporter must derive outgoing moment-shape controls.
- source: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle_open_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: PN_5_FULL_BUNDLE_SURFACE
  relation: REQUIRES
  status: open
  note: Higher-order full-bundle surface needs outgoing moment-shape control.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open
- topic/moving_throat
- type/moment_shape_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Outgoing moment-shape control

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_N2_N4_MOMENT_SHAPE`
> **Status:** `open`
> **Layer:** `open_gate`
> **Type:** `moment_shape_gate`

## Summary

Derive independent outgoing N2 and N4 moment-shape controls on the same frozen branch.

## What Remains Open

Derive independent outgoing N2 and N4 moment-shape controls on the same frozen branch.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Derive independent outgoing N2 and N4 moment-shape controls on the same frozen branch.

## Mathematical Role

- Layer: `open_gate`
- Type: `moment_shape_gate`
- Status: `open`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]

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
| `DECOMPOSES_INTO` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Exporter must derive outgoing moment-shape controls. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Claim feeds this downstream object, output, or open gate. |
| `REQUIRES` | [[PN_5_FULL_BUNDLE_SURFACE]] | Higher-order full-bundle surface needs outgoing moment-shape control. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_pde_program_compact.md`
- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
