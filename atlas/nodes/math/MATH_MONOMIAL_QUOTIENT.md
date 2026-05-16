---
id: MATH_MONOMIAL_QUOTIENT
title: Monomial quotient coordinates
type: quotient_map
layer: math_object
status: exact_within_coherent_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Three exact quotient coordinates separate physical branch motion from the five-dimensional similarity orbit.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
math_ids:
- MATH_XI1_PREF_SLOPE
equation_ids:
- EQ_MONOMIAL_QUOTIENT
outgoing_edges:
- target: MATH_XI1_PREF_SLOPE
  relation: ORGANIZES
  status: exact_within_branch
  note: Xi1 is a direct defect coordinate/prefactor slope within the quotient packet.
incoming_edges:
- source: EQ_MONOMIAL_QUOTIENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_5PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_MONOMIAL_QUOTIENT.
- source: BACKLINK_G2_OUTPUT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_MONOMIAL_QUOTIENT.
- source: MT_V2_18_MONOMIAL_ORBIT
  relation: DERIVES
  status: exact_within_branch
  note: Derives quotient coordinates and similarity-orbit split.
- source: ANOMALY_G2_COMMON_QUOTIENT
  relation: IMPORTS
  status: conditional
  note: g-2 common layer uses the same quotient coordinates.
- source: PN_5_FULL_BUNDLE_SURFACE
  relation: USES
  status: reduced
  note: 5PN weak-axisymmetric transport uses monomial quotient/orbit-lock machinery.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_within_coherent_branch
- topic/moving_throat
- type/quotient_map
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Monomial quotient coordinates

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_MONOMIAL_QUOTIENT`
> **Status:** `exact_within_coherent_branch`
> **Layer:** `math_object`
> **Type:** `quotient_map`

## Summary

Three exact quotient coordinates separate physical branch motion from the five-dimensional similarity orbit.

## Physical Meaning

Three exact quotient coordinates separate physical branch motion from the five-dimensional similarity orbit.

## Mathematical Role

- Layer: `math_object`
- Type: `quotient_map`
- Status: `exact_within_coherent_branch`

## Equation

$$
q_tr
$$

$$
q_nt
$$

$$
q_eta
$$

$$
rank(M*)=3
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_XI1_PREF_SLOPE]]

### Related equations
- [[EQ_MONOMIAL_QUOTIENT]]

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ORGANIZES` | [[MATH_XI1_PREF_SLOPE]] | Xi1 is a direct defect coordinate/prefactor slope within the quotient packet. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_MONOMIAL_QUOTIENT]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_5PN_FULL]] | Paper backlink block references MATH_MONOMIAL_QUOTIENT. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_G2_OUTPUT]] | Paper backlink block references MATH_MONOMIAL_QUOTIENT. |
| `DERIVES` | [[MT_V2_18_MONOMIAL_ORBIT]] | Derives quotient coordinates and similarity-orbit split. |
| `IMPORTS` | [[ANOMALY_G2_COMMON_QUOTIENT]] | g-2 common layer uses the same quotient coordinates. |
| `USES` | [[PN_5_FULL_BUNDLE_SURFACE]] | 5PN weak-axisymmetric transport uses monomial quotient/orbit-lock machinery. |

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
