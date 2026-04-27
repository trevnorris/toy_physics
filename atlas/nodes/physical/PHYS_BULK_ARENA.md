---
id: PHYS_BULK_ARENA
title: 4+1 bulk arena
type: arena
layer: physical_ontology
status: exact_declared_parent
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Fundamental spacetime arena with brane coordinates plus transverse/bulk coordinate w.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/moving_throat_pde_program_compact.md
- 4d_summary.md
- moving_throat_pde_program_compact.md
legacy_sources:
- 4d_summary.md
- moving_throat_pde_program_compact.md
math_ids:
- MATH_PARENT_ACTION_CURRENT
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
outgoing_edges:
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_parent_current
  note: Physical ontology object grounded by this claim.
- target: MATH_PARENT_ACTION_CURRENT
  relation: HOSTS
  status: exact
  note: Parent action lives on the 4+1 bulk arena.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_declared_parent
- topic/moving_throat
- type/arena
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 4+1 bulk arena

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_BULK_ARENA`  
> **Status:** `exact_declared_parent`  
> **Layer:** `physical_ontology`  
> **Type:** `arena`

## Summary

Fundamental spacetime arena with brane coordinates plus transverse/bulk coordinate w.

## Physical Meaning

Fundamental spacetime arena with brane coordinates plus transverse/bulk coordinate w.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `arena`
- Status: `exact_declared_parent`

## Equation

$$
x^M=(t,x,y,z,w)
$$

$$
X=(x,y,z,w)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_PARENT_ACTION_CURRENT]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Physical ontology object grounded by this claim. |
| `HOSTS` | [[MATH_PARENT_ACTION_CURRENT]] | Parent action lives on the 4+1 bulk arena. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `notes/moving_throat_pde_program_compact.md`
- `4d_summary.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
