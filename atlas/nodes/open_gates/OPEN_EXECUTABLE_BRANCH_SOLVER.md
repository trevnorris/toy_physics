---
id: OPEN_EXECUTABLE_BRANCH_SOLVER
title: Executable physical branch solver
type: solver_gate
layer: open_gate
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Numerically or analytically solve the actual open-throat branch and export frozen profile/response data.
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
- '[[SEC_PDE_WEAK_FORM_EXPORTER]]'
math_ids:
- MATH_WEAK_FORM_BRANCH_EXTRACTOR
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
source_ids:
- SEC_PDE_WEAK_FORM_EXPORTER
incoming_edges:
- source: SEC_PDE_WEAK_FORM_EXPORTER
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Weak-form branch extraction preparation.
- source: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: FEEDS_OR_STATUS_OF
  status: open_actual_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: BACKLINK_5PN_FULL
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_EXECUTABLE_BRANCH_SOLVER.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_EXECUTABLE_BRANCH_SOLVER.
- source: BACKLINK_PDE_AUDIT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_EXECUTABLE_BRANCH_SOLVER.
- source: ATLAS_CURRENT_READINESS_V06
  relation: FLAGS_REMAINING_GATE
  status: v06
  note: Still open after v0.6 organization pass.
- source: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: NEEDS
  status: open
  note: Needs branch solver or analytic branch exporter.
- source: MATH_WEAK_FORM_BRANCH_EXTRACTOR
  relation: PREPARES
  status: open
  note: Schema is needed by executable branch solver.
- source: MT_V2_21_BRANCH_FIXTURE
  relation: TESTS
  status: implemented
  note: Fixture defines observable packets and built-in test branches.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open
- topic/moving_throat
- type/solver_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Executable physical branch solver

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_EXECUTABLE_BRANCH_SOLVER`  
> **Status:** `open`  
> **Layer:** `open_gate`  
> **Type:** `solver_gate`

## Summary

Numerically or analytically solve the actual open-throat branch and export frozen profile/response data.

## What Remains Open

Numerically or analytically solve the actual open-throat branch and export frozen profile/response data.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Numerically or analytically solve the actual open-throat branch and export frozen profile/response data.

## Mathematical Role

- Layer: `open_gate`
- Type: `solver_gate`
- Status: `open`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WEAK_FORM_BRANCH_EXTRACTOR]]

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- [[SEC_PDE_WEAK_FORM_EXPORTER]]

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_WEAK_FORM_EXPORTER]] | Weak-form branch extraction preparation. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Claim feeds this downstream object, output, or open gate. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_5PN_FULL]] | Paper backlink block flags open gate OPEN_EXECUTABLE_BRANCH_SOLVER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block flags open gate OPEN_EXECUTABLE_BRANCH_SOLVER. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block flags open gate OPEN_EXECUTABLE_BRANCH_SOLVER. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V06]] | Still open after v0.6 organization pass. |
| `NEEDS` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Needs branch solver or analytic branch exporter. |
| `PREPARES` | [[MATH_WEAK_FORM_BRANCH_EXTRACTOR]] | Schema is needed by executable branch solver. |
| `TESTS` | [[MT_V2_21_BRANCH_FIXTURE]] | Fixture defines observable packets and built-in test branches. |

## Source Anchors

### Source anchor notes
- [[SEC_PDE_WEAK_FORM_EXPORTER]]

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_pde_program_compact.md`
- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
