---
id: PHYS_OPEN_FINITE_EXIT
title: Open finite-radius exit
type: geometry_boundary
layer: physical_ontology
status: effective_branch_condition
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The physically preferred throat branch has a finite exit into the bulk/reservoir side, R(L)>0, not a hard cap R(L)=0.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
source_links:
- '[[FILE_PDE_AUDIT]]'
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
source_ids:
- FILE_PDE_AUDIT
incoming_edges:
- source: FILE_PDE_AUDIT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: HAS_BOUNDARY_FEATURE
  status: physical
  note: Open finite exit is part of the updated physical picture.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/effective_branch_condition
- topic/moving_throat
- type/geometry_boundary
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Open finite-radius exit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_OPEN_FINITE_EXIT`  
> **Status:** `effective_branch_condition`  
> **Layer:** `physical_ontology`  
> **Type:** `geometry_boundary`

## Summary

The physically preferred throat branch has a finite exit into the bulk/reservoir side, R(L)>0, not a hard cap R(L)=0.

## Physical Meaning

The physically preferred throat branch has a finite exit into the bulk/reservoir side, R(L)>0, not a hard cap R(L)=0.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `geometry_boundary`
- Status: `effective_branch_condition`

## Equation

$$
R(0)=a
$$

$$
R(L)>0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_PDE_AUDIT]] | File anchor documents this node. |
| `HAS_BOUNDARY_FEATURE` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Open finite exit is part of the updated physical picture. |

## Source Anchors

### Source anchor notes
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
