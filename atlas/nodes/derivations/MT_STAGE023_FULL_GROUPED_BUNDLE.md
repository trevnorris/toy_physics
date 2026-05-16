---
id: MT_STAGE023_FULL_GROUPED_BUNDLE
title: Stage 023 full grouped bundle
type: moving_throat_stage
layer: derivation
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Writes full 20/21/22 bundle, projectors, anisotropy defects, and target ratio.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_023.tex
- notes/moving_throat_pde_program_compact.md
- moving_throat_pde_stage023_full_grouped_bundle.md
legacy_sources:
- moving_throat_pde_stage023_full_grouped_bundle.md
math_ids:
- MATH_GROUPED_PROJECTORS_GGRP
claim_ids:
- CLAIM_STAGE023_FULL_BUNDLE_RATIO
outgoing_edges:
- target: READOUT_D0_C_P0_N2_N4
  relation: COMPRESSES_TO
  status: target packet
  note: Full bundle exports compressed audit readouts.
- target: MT_STAGE024_OVERLAP_ISOTROPY
  relation: SPECIALIZES_TO
  status: exact within O3 reduced kernel
  note: O(3)-invariant overlap model collapses grouped lanes.
- target: PN_5_FULL_BUNDLE_SURFACE
  relation: SUPPORTS
  status: reduced
  note: Full bundle target surface matches later 5PN/orbit-lock compression.
incoming_edges:
- source: CLAIM_STAGE023_FULL_BUNDLE_RATIO
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_STAGE022_GROUPED_P2_BRIDGE
  relation: GENERALIZES_TO
  status: exact within bundle
  note: Full 20/21/22 bundle and projectors.
- source: MATH_GROUPED_PROJECTORS_GGRP
  relation: ORGANIZES
  status: exact
  note: Organizes the full grouped P2 bundle and anisotropy extraction.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/exact_within_reduced_bundle
- topic/moving_throat
- topic/quadrupole
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 023 full grouped bundle

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE023_FULL_GROUPED_BUNDLE`
> **Status:** `exact_within_reduced_bundle`
> **Layer:** `derivation`
> **Type:** `moving_throat_stage`

## Summary

Writes full 20/21/22 bundle, projectors, anisotropy defects, and target ratio.

## Physical Meaning

Writes full 20/21/22 bundle, projectors, anisotropy defects, and target ratio.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `exact_within_reduced_bundle`

## Equation

$$
D0=K-B0-Z0
$$

$$
P0=N0/D0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_GROUPED_PROJECTORS_GGRP]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE023_FULL_BUNDLE_RATIO]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `COMPRESSES_TO` | [[READOUT_D0_C_P0_N2_N4]] | Full bundle exports compressed audit readouts. |
| `SPECIALIZES_TO` | [[MT_STAGE024_OVERLAP_ISOTROPY]] | O(3)-invariant overlap model collapses grouped lanes. |
| `SUPPORTS` | [[PN_5_FULL_BUNDLE_SURFACE]] | Full bundle target surface matches later 5PN/orbit-lock compression. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE023_FULL_BUNDLE_RATIO]] | Claim feeds this downstream object, output, or open gate. |
| `GENERALIZES_TO` | [[MT_STAGE022_GROUPED_P2_BRIDGE]] | Full 20/21/22 bundle and projectors. |
| `ORGANIZES` | [[MATH_GROUPED_PROJECTORS_GGRP]] | Organizes the full grouped P2 bundle and anisotropy extraction. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/pde_ledger/paper/stages/stage_023.tex`
- `notes/moving_throat_pde_program_compact.md`
- `moving_throat_pde_stage023_full_grouped_bundle.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
