---
id: MT_STAGE7_OVERLAP_ISOTROPY
title: Stage 7 overlap isotropy theorem
type: moving_throat_stage
layer: derivation
status: exact_within_O3_reduced_kernel
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Proves angular source map identity and O(3) isotropy; first weak-axisymmetric splitting law.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- notes/pde_audit_full.md
- moving_throat_output_full.md
- pde_audit_full.md
legacy_sources:
- moving_throat_output_full.md
- pde_audit_full.md:V2-12/V2-17
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
claim_ids:
- CLAIM_STAGE7_O3_ISOTROPY
outgoing_edges:
- target: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: GIVES_PHYSICAL_DIAGNOSTIC
  status: reduced
  note: O(3) isotropy and weak-axisymmetric b=3a pattern classify shape response.
- target: MT_STAGE8_MINIMAL_NORMALIZATION
  relation: NARROWS_TO
  status: reduced
  note: Angular side closed; radial/axial normalization remains.
incoming_edges:
- source: CLAIM_STAGE7_O3_ISOTROPY
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_STAGE6_FULL_GROUPED_BUNDLE
  relation: SPECIALIZES_TO
  status: exact within O3 reduced kernel
  note: O(3)-invariant overlap model collapses grouped lanes.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/exact_within_o3_reduced_kernel
- topic/moving_throat
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 7 overlap isotropy theorem

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE7_OVERLAP_ISOTROPY`  
> **Status:** `exact_within_O3_reduced_kernel`  
> **Layer:** `derivation`  
> **Type:** `moving_throat_stage`

## Summary

Proves angular source map identity and O(3) isotropy; first weak-axisymmetric splitting law.

## Physical Meaning

Proves angular source map identity and O(3) isotropy; first weak-axisymmetric splitting law.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `exact_within_O3_reduced_kernel`

## Equation

$$
mhat_ang=1
$$

$$
(20,21,22)~(1,1/2,-1)
$$

$$
b=3a
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_STAGE7_O3_ISOTROPY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GIVES_PHYSICAL_DIAGNOSTIC` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | O(3) isotropy and weak-axisymmetric b=3a pattern classify shape response. |
| `NARROWS_TO` | [[MT_STAGE8_MINIMAL_NORMALIZATION]] | Angular side closed; radial/axial normalization remains. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Claim feeds this downstream object, output, or open gate. |
| `SPECIALIZES_TO` | [[MT_STAGE6_FULL_GROUPED_BUNDLE]] | O(3)-invariant overlap model collapses grouped lanes. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/moving_throat_notes_full.md`
- `notes/pde_audit_full.md`
- `moving_throat_output_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
