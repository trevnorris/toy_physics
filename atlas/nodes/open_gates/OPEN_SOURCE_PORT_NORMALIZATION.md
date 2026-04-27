---
id: OPEN_SOURCE_PORT_NORMALIZATION
title: Source/port normalization law
type: normalization_gate
layer: open_gate
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Derive mhat and P0 as branch data rather than fitted scales.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
math_ids:
- MATH_STF_SOURCE_MAP
claim_ids:
- CLAIM_STAGE8_TO_14_SELECTED_BRANCH
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
status_firewall_ids:
- FIREWALL_ANGULAR_NOT_RADIAL
incoming_edges:
- source: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: DECOMPOSES_INTO
  status: open
  note: Exporter must derive source/port normalization law.
- source: CLAIM_STAGE8_TO_14_SELECTED_BRANCH
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_selected_reduced_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: BACKLINK_25PN
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_SOURCE_PORT_NORMALIZATION.
- source: MATH_STF_SOURCE_MAP
  relation: NARROWS
  status: open
  note: Leaves radial/axial source/port normalization only.
- source: FIREWALL_ANGULAR_NOT_RADIAL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open
- topic/moving_throat
- topic/quadrupole
- type/normalization_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Source/port normalization law

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_SOURCE_PORT_NORMALIZATION`  
> **Status:** `open`  
> **Layer:** `open_gate`  
> **Type:** `normalization_gate`

## Summary

Derive mhat and P0 as branch data rather than fitted scales.

## What Remains Open

Derive mhat and P0 as branch data rather than fitted scales.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Derive mhat and P0 as branch data rather than fitted scales.

## Mathematical Role

- Layer: `open_gate`
- Type: `normalization_gate`
- Status: `open`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_STF_SOURCE_MAP]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- [[FIREWALL_ANGULAR_NOT_RADIAL]]

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DECOMPOSES_INTO` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Exporter must derive source/port normalization law. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]] | Claim feeds this downstream object, output, or open gate. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_25PN]] | Paper backlink block flags open gate OPEN_SOURCE_PORT_NORMALIZATION. |
| `NARROWS` | [[MATH_STF_SOURCE_MAP]] | Leaves radial/axial source/port normalization only. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ANGULAR_NOT_RADIAL]] | Firewall preserves this correct status boundary. |

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
