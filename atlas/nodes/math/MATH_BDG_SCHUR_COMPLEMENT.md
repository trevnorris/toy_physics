---
id: MATH_BDG_SCHUR_COMPLEMENT
title: BdG-wall Schur complement
type: reduced_kernel
layer: math_object
status: exact_within_stable_mode_reduction
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Stable BdG modes lower static stiffness, increase inertia, and generate even response moments by Schur complement.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
equation_ids:
- EQ_BDG_SCHUR_KERNEL
claim_ids:
- CLAIM_STAGE3_BDG_SCHUR
outgoing_edges:
- target: MT_STAGE3_BDG_COUPLING
  relation: FORMALIZES
  status: reduced
  note: Matches Stage-3 wall/BdG coupling kernel.
incoming_edges:
- source: EQ_BDG_SCHUR_KERNEL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_BDG_SCHUR_COMPLEMENT.
- source: CLAIM_STAGE3_BDG_SCHUR
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_stable_modes
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_08_BDG_SCHUR
  relation: VALIDATES
  status: reduced
  note: Stable BdG modes produce Schur-complement support moments.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_within_stable_mode_reduction
- topic/moving_throat
- topic/projection
- type/reduced_kernel
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# BdG-wall Schur complement

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_BDG_SCHUR_COMPLEMENT`
> **Status:** `exact_within_stable_mode_reduction`
> **Layer:** `math_object`
> **Type:** `reduced_kernel`

## Summary

Stable BdG modes lower static stiffness, increase inertia, and generate even response moments by Schur complement.

## Physical Meaning

Stable BdG modes lower static stiffness, increase inertia, and generate even response moments by Schur complement.

## Mathematical Role

- Layer: `math_object`
- Type: `reduced_kernel`
- Status: `exact_within_stable_mode_reduction`

## Equation

$$
D_eff=K-omega^2 M-C(Omega^2-omega^2 I)^-1 C^T
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_BDG_SCHUR_KERNEL]]

### Related claims
- [[CLAIM_STAGE3_BDG_SCHUR]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FORMALIZES` | [[MT_STAGE3_BDG_COUPLING]] | Matches Stage-3 wall/BdG coupling kernel. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_BDG_SCHUR_KERNEL]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references MATH_BDG_SCHUR_COMPLEMENT. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE3_BDG_SCHUR]] | Claim feeds this downstream object, output, or open gate. |
| `VALIDATES` | [[MT_V2_08_BDG_SCHUR]] | Stable BdG modes produce Schur-complement support moments. |

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
