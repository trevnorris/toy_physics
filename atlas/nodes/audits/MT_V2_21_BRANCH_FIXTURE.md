---
id: MT_V2_21_BRANCH_FIXTURE
title: V2-21 branch-extraction fixture
type: fixture_gate
layer: status_audit
status: implemented
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Defines primitive/direct coefficient inputs, open-throat gate, extraction formulas, grouped decomposition, residuals, and stability gates.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
open_gate_ids:
- OPEN_EXECUTABLE_BRANCH_SOLVER
outgoing_edges:
- target: OPEN_EXECUTABLE_BRANCH_SOLVER
  relation: TESTS
  status: implemented
  note: Fixture defines observable packets and built-in test branches.
incoming_edges:
- source: MT_V2_22A_PROFILE_ADAPTER
  relation: FEEDS
  status: implemented
  note: Profile adapter supplies coefficients to fixture extraction.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/implemented
- topic/moving_throat
- type/fixture_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-21 branch-extraction fixture

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_21_BRANCH_FIXTURE`  
> **Status:** `implemented`  
> **Layer:** `status_audit`  
> **Type:** `fixture_gate`

## Summary

Defines primitive/direct coefficient inputs, open-throat gate, extraction formulas, grouped decomposition, residuals, and stability gates.

## Physical Meaning

Defines primitive/direct coefficient inputs, open-throat gate, extraction formulas, grouped decomposition, residuals, and stability gates.

## Mathematical Role

- Layer: `status_audit`
- Type: `fixture_gate`
- Status: `implemented`

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
- [[OPEN_EXECUTABLE_BRANCH_SOLVER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `TESTS` | [[OPEN_EXECUTABLE_BRANCH_SOLVER]] | Fixture defines observable packets and built-in test branches. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MT_V2_22A_PROFILE_ADAPTER]] | Profile adapter supplies coefficients to fixture extraction. |

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
