---
id: EQ_BDG_SCHUR_KERNEL
title: BdG Schur-complement kernel
type: equation
layer: equation_anchor
status: exact_within_stable_mode_reduction
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Stable BdG support renormalizes wall stiffness/inertia and generates even moments.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_003.tex
- notes/pde_audit_full.md
- moving_throat_pde_stage003_bdg_coupling.md
- pde_audit_full.md
legacy_sources:
- moving_throat_pde_stage003_bdg_coupling.md
- pde_audit_full.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
math_ids:
- MATH_BDG_SCHUR_COMPLEMENT
equation_ids:
- EQ_GROUPED_RESPONSE_MOMENTS
- EQ_WALL_P2_STIFFNESS
claim_ids:
- CLAIM_STAGE3_BDG_SCHUR
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_BDG_SCHUR_COMPLEMENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_GROUPED_RESPONSE_MOMENTS
  relation: FEEDS
  status: reduced
  note: BdG-dressed operator moments enter normalized response.
- target: CLAIM_STAGE3_BDG_SCHUR
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_stable_modes
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_BDG_SCHUR_KERNEL.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_BDG_SCHUR_KERNEL.
- source: EQ_WALL_P2_STIFFNESS
  relation: FEEDS
  status: reduced
  note: Bare P2 wall operator is then dressed by BdG support.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_stable_mode_reduction
- topic/moving_throat
- topic/projection
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# BdG Schur-complement kernel

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_BDG_SCHUR_KERNEL`  
> **Status:** `exact_within_stable_mode_reduction`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Stable BdG support renormalizes wall stiffness/inertia and generates even moments.

## Physical Meaning

Stable BdG support renormalizes wall stiffness/inertia and generates even moments.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_stable_mode_reduction`
- Parent node: [[MATH_BDG_SCHUR_COMPLEMENT]]

## Equation

$$
D_eff(ω)=K-ω²M-C(Ω²-ω²I)^(-1)Cᵀ
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_BDG_SCHUR_COMPLEMENT]]

### Related equations
- [[EQ_GROUPED_RESPONSE_MOMENTS]]
- [[EQ_WALL_P2_STIFFNESS]]

### Related claims
- [[CLAIM_STAGE3_BDG_SCHUR]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_BDG_SCHUR_COMPLEMENT]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_GROUPED_RESPONSE_MOMENTS]] | BdG-dressed operator moments enter normalized response. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE3_BDG_SCHUR]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_BDG_SCHUR_KERNEL. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_BDG_SCHUR_KERNEL. |
| `FEEDS` | [[EQ_WALL_P2_STIFFNESS]] | Bare P2 wall operator is then dressed by BdG support. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/pde_ledger/paper/stages/stage_003.tex`
- `notes/pde_audit_full.md`
- `moving_throat_pde_stage003_bdg_coupling.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
