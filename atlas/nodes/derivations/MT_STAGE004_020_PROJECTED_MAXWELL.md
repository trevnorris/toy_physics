---
id: MT_STAGE004_020_PROJECTED_MAXWELL
title: Stages 004--020 projected Maxwell packet
type: moving_throat_stage
layer: derivation
status: exact_projection_derivation
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Imports the derivation-only projected Maxwell chain from notes/em_projected through step 18. This is the canonical EM derivation block; reduction is no longer the primary EM der...
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_004.tex
- research/pde_ledger/paper/stages/stage_020.tex
- notes/moving_throat_pde_program_compact.md
- notes/em_projected/step_01_projected_maxwell_readme.md
- notes/em_projected/step_18_parent_throat_action_weak_axisym_packet_notes.md
legacy_sources:
- notes/em_projected/step_01_projected_maxwell_readme.md
- notes/em_projected/step_18_parent_throat_action_weak_axisym_packet_notes.md
math_ids:
- MATH_MAXWELL_MIXED_KERNEL
- MATH_MIXED_FIELDS_EW_CA
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
outgoing_edges:
- target: MT_STAGE021_REDUCED_ONE_PORT
  relation: CLOSED_BY
  status: retained_reduced_adapter
  note: Projection-first EM derivation is closed for downstream response work by the retained one-port normal form.
incoming_edges:
- source: MATH_MIXED_FIELDS_EW_CA
  relation: ENABLES
  status: exact/reduced
  note: Mixed fields are the microscopic place for outgoing bridge.
- source: MT_STAGE3_BDG_COUPLING
  relation: FEEDS
  status: exact_projection_derivation
  note: BdG support data feeds the projection-first EM packet before the retained one-port adapter.
- source: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_mixed_kernel
  note: Claim feeds this downstream object, output, or open gate.
- source: MATH_MAXWELL_MIXED_KERNEL
  relation: FORMALIZES
  status: reduced
  note: Matches Stage-4 port-active mixed-sector bridge.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/exact_projection_derivation
- topic/maxwell
- topic/moving_throat
- topic/projection
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stages 004--020 projected Maxwell packet

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE004_020_PROJECTED_MAXWELL`
> **Status:** `exact_projection_derivation`
> **Layer:** `derivation`
> **Type:** `moving_throat_stage`

## Summary

Imports the derivation-only projected Maxwell chain from notes/em_projected through step 18. This is the canonical EM derivation block; reduction is no longer the primary EM derivation.

## Physical Meaning

Imports the derivation-only projected Maxwell chain from notes/em_projected through step 18. This is the canonical EM derivation block; reduction is no longer the primary EM derivation.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `exact_projection_derivation`

## Equation

$$
projection-first Maxwell sector
$$

$$
derivation notes step 01 through step 18, with no step 06 source file
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MAXWELL_MIXED_KERNEL]]
- [[MATH_MIXED_FIELDS_EW_CA]]

### Related equations
- none

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `CLOSED_BY` | [[MT_STAGE021_REDUCED_ONE_PORT]] | Projection-first EM derivation is closed for downstream response work by the retained one-port normal form. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ENABLES` | [[MATH_MIXED_FIELDS_EW_CA]] | Mixed fields are the microscopic place for outgoing bridge. |
| `FEEDS` | [[MT_STAGE3_BDG_COUPLING]] | BdG support data feeds the projection-first EM packet before the retained one-port adapter. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Claim feeds this downstream object, output, or open gate. |
| `FORMALIZES` | [[MATH_MAXWELL_MIXED_KERNEL]] | Matches Stage-4 port-active mixed-sector bridge. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/pde_ledger/paper/stages/stage_004.tex`
- `research/pde_ledger/paper/stages/stage_020.tex`
- `notes/moving_throat_pde_program_compact.md`
- `notes/em_projected/step_01_projected_maxwell_readme.md`
- `notes/em_projected/step_18_parent_throat_action_weak_axisym_packet_notes.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
