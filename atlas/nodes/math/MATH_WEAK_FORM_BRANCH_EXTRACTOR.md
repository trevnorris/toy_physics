---
id: MATH_WEAK_FORM_BRANCH_EXTRACTOR
title: Weak-form branch extractor
type: solver_schema
layer: math_object
status: prepared_open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Galerkin/weak-form extraction protocol for branch data K,M,B_n,Z_n,N_n, grouped projectors, residual packet, and stability gates.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
open_gate_ids:
- OPEN_EXECUTABLE_BRANCH_SOLVER
outgoing_edges:
- target: OPEN_EXECUTABLE_BRANCH_SOLVER
  relation: PREPARES
  status: open
  note: Schema is needed by executable branch solver.
incoming_edges:
- source: MT_V2_20_WEAK_FORM_EXTRACTOR
  relation: PREPARES
  status: prepared
  note: Weak-form/Galerkin extraction schema.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/prepared_open
- topic/moving_throat
- type/solver_schema
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weak-form branch extractor

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_WEAK_FORM_BRANCH_EXTRACTOR`
> **Status:** `prepared_open`
> **Layer:** `math_object`
> **Type:** `solver_schema`

## Summary

Galerkin/weak-form extraction protocol for branch data K,M,B_n,Z_n,N_n, grouped projectors, residual packet, and stability gates.

## Physical Meaning

Galerkin/weak-form extraction protocol for branch data K,M,B_n,Z_n,N_n, grouped projectors, residual packet, and stability gates.

## Mathematical Role

- Layer: `math_object`
- Type: `solver_schema`
- Status: `prepared_open`

## Equation

$$
profiles -> coefficients -> residual packet
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
- [[OPEN_EXECUTABLE_BRANCH_SOLVER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `PREPARES` | [[OPEN_EXECUTABLE_BRANCH_SOLVER]] | Schema is needed by executable branch solver. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `PREPARES` | [[MT_V2_20_WEAK_FORM_EXTRACTOR]] | Weak-form/Galerkin extraction schema. |

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
