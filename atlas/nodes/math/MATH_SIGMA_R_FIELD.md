---
id: MATH_SIGMA_R_FIELD
title: Moving surface field Sigma/R
type: geometry_field
layer: math_object
status: effective_closure_unless_Sigma_action_promoted
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Distributed throat shape variable replacing pure a,L closure.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- notes/pde_audit_full.md
- moving_throat_pde_stage025_minimal_isotropic_normalization.md
- pde_audit_full.md
legacy_sources:
- moving_throat_pde_stage025_minimal_isotropic_normalization.md
- pde_audit_full.md:V2-28
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
- PHYS_REG_PARTICLE_THROAT
math_ids:
- MATH_PARENT_ACTION_CURRENT
claim_ids:
- CLAIM_STAGE1_GEOMETRY_LIFT
outgoing_edges:
- target: MT_STAGE1_GEOMETRY_LIFT
  relation: ORGANIZES
  status: effective
  note: Stage 1 uses Sigma/R to create distributed geometry field.
incoming_edges:
- source: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: FEEDS_OR_STATUS_OF
  status: effective_geometry_lift
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_REG_PARTICLE_THROAT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: REPRESENTED_BY
  status: effective closure
  note: Finite throat geometry represented by Sigma=r-R.
- source: MATH_PARENT_ACTION_CURRENT
  relation: USES_AS_COUPLING_ARGUMENT
  status: exact coupling argument
  note: Geometry currently enters through V_conf(Sigma).
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/effective_closure_unless_sigma_action_promoted
- topic/moving_throat
- type/geometry_field
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Moving surface field Sigma/R

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_SIGMA_R_FIELD`
> **Status:** `effective_closure_unless_Sigma_action_promoted`
> **Layer:** `math_object`
> **Type:** `geometry_field`

## Summary

Distributed throat shape variable replacing pure a,L closure.

## Physical Meaning

Distributed throat shape variable replacing pure a,L closure.

## Mathematical Role

- Layer: `math_object`
- Type: `geometry_field`
- Status: `effective_closure_unless_Sigma_action_promoted`

## Equation

$$
Sigma=r-R(Omega,w,t)
$$

$$
eta=R-R0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]
- [[PHYS_REG_PARTICLE_THROAT]]

### Related math nodes
- [[MATH_PARENT_ACTION_CURRENT]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ORGANIZES` | [[MT_STAGE1_GEOMETRY_LIFT]] | Stage 1 uses Sigma/R to create distributed geometry field. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Claim feeds this downstream object, output, or open gate. |
| `LINKS_TO` | [[PHYS_REG_PARTICLE_THROAT]] | Physical register entry links to graph object. |
| `REPRESENTED_BY` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Finite throat geometry represented by Sigma=r-R. |
| `USES_AS_COUPLING_ARGUMENT` | [[MATH_PARENT_ACTION_CURRENT]] | Geometry currently enters through V_conf(Sigma). |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/moving_throat_notes_full.md`
- `notes/pde_audit_full.md`
- `moving_throat_pde_stage025_minimal_isotropic_normalization.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
