---
id: MATH_GROUPED_PROJECTORS_GGRP
title: Grouped P2 weighted projectors
type: projector_calculus
layer: math_object
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Weighted grouped-real P2 trace/anomaly decomposition with metric diag(1,2,2).
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
- CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
outgoing_edges:
- target: MT_STAGE023_FULL_GROUPED_BUNDLE
  relation: ORGANIZES
  status: exact
  note: Organizes the full grouped P2 bundle and anisotropy extraction.
incoming_edges:
- source: BACKLINK_3PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_GROUPED_PROJECTORS_GGRP.
- source: CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_11_GROUPED_PROJECTORS
  relation: FREEZES
  status: exact
  note: Freezes grouped P2 trace/anomaly projectors.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact
- topic/g2
- topic/moving_throat
- type/projector_calculus
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Grouped P2 weighted projectors

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_GROUPED_PROJECTORS_GGRP`
> **Status:** `exact`
> **Layer:** `math_object`
> **Type:** `projector_calculus`

## Summary

Weighted grouped-real P2 trace/anomaly decomposition with metric diag(1,2,2).

## Physical Meaning

Weighted grouped-real P2 trace/anomaly decomposition with metric diag(1,2,2).

## Mathematical Role

- Layer: `math_object`
- Type: `projector_calculus`
- Status: `exact`

## Equation

$$
Ggrp=diag(1,2,2)
$$

$$
P_bar,P_a,P_b
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
- [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ORGANIZES` | [[MT_STAGE023_FULL_GROUPED_BUNDLE]] | Organizes the full grouped P2 bundle and anisotropy extraction. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_3PN_FULL]] | Paper backlink block references MATH_GROUPED_PROJECTORS_GGRP. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_3PN_GROUPED_P2_WITHIN_HIERARCHY]] | Claim feeds this downstream object, output, or open gate. |
| `FREEZES` | [[MT_V2_11_GROUPED_PROJECTORS]] | Freezes grouped P2 trace/anomaly projectors. |

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
