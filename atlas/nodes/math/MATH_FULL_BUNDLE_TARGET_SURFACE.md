---
id: MATH_FULL_BUNDLE_TARGET_SURFACE
title: Isotropic full-bundle target surface
type: target_surface
layer: math_object
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Combined surface for one-pole conservative response, universal quadrupole normalization, and constant-prefactor branch.
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
- EQ_FULL_BUNDLE_TARGET_SURFACE
claim_ids:
- CLAIM_STAGE6_FULL_BUNDLE_RATIO
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: REQUIRES
  status: open
  note: Target surface can only be tested after branch exporter supplies frozen data.
incoming_edges:
- source: EQ_FULL_BUNDLE_TARGET_SURFACE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: MT_V2_13_GROUPED_NORMALIZATION
  relation: DEFINES
  status: exact_within_bundle
  note: Defines normalized response, P0, P2, P4, target surface.
- source: CLAIM_STAGE6_FULL_BUNDLE_RATIO
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_19_TARGET_SURFACE
  relation: FREEZES
  status: exact_within_bundle
  note: Freezes the target packet for actual-branch export.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_within_reduced_bundle
- topic/moving_throat
- topic/quadrupole
- type/target_surface
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Isotropic full-bundle target surface

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_FULL_BUNDLE_TARGET_SURFACE`  
> **Status:** `exact_within_reduced_bundle`  
> **Layer:** `math_object`  
> **Type:** `target_surface`

## Summary

Combined surface for one-pole conservative response, universal quadrupole normalization, and constant-prefactor branch.

## Physical Meaning

Combined surface for one-pole conservative response, universal quadrupole normalization, and constant-prefactor branch.

## Mathematical Role

- Layer: `math_object`
- Type: `target_surface`
- Status: `exact_within_reduced_bundle`

## Equation

$$
D0(B4+Z4)=3(M+B2+Z2)^2
$$

$$
mhat0^2 N0/D0=54 G c_s^5/(5a^5c^5)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_FULL_BUNDLE_TARGET_SURFACE]]

### Related claims
- [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REQUIRES` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Target surface can only be tested after branch exporter supplies frozen data. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_FULL_BUNDLE_TARGET_SURFACE]] | Equation anchor belongs to or formalizes this graph node. |
| `DEFINES` | [[MT_V2_13_GROUPED_NORMALIZATION]] | Defines normalized response, P0, P2, P4, target surface. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]] | Claim feeds this downstream object, output, or open gate. |
| `FREEZES` | [[MT_V2_19_TARGET_SURFACE]] | Freezes the target packet for actual-branch export. |

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
