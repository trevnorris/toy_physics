---
id: MT_STAGE022_GROUPED_P2_BRIDGE
title: Stage 022 grouped P2 normalization bridge
type: moving_throat_stage
layer: derivation
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Converts conservative operator moments into normalized grouped response and isolates P0=N0/D0.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_022.tex
- notes/moving_throat_pde_program_compact.md
- moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
legacy_sources:
- moving_throat_pde_stage022_grouped_p2_normalization_bridge.md
equation_ids:
- EQ_GROUPED_RESPONSE_MOMENTS
claim_ids:
- CLAIM_STAGE022_P0_TARGET
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
outgoing_edges:
- target: MT_STAGE023_FULL_GROUPED_BUNDLE
  relation: GENERALIZES_TO
  status: exact within bundle
  note: Full 20/21/22 bundle and projectors.
- target: OPEN_QUAD_NORMALIZATION
  relation: ISOLATES
  status: open
  note: Invariant product mhat0²P0 is the remaining normalization target.
- target: PN_2_5_QUAD_NARROWING
  relation: SUPPORTS
  status: reduced
  note: Moving-throat bridge supplies exact reduced target mhat0²P0.
incoming_edges:
- source: EQ_GROUPED_RESPONSE_MOMENTS
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: MT_STAGE021_REDUCED_ONE_PORT
  relation: FEEDS
  status: exact_closure_adapter
  note: One-lane outgoing transfer lifted to grouped P2 normalization language.
- source: CLAIM_STAGE022_P0_TARGET
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_grouped_bridge
  note: Claim feeds this downstream object, output, or open gate.
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

# Stage 022 grouped P2 normalization bridge

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE022_GROUPED_P2_BRIDGE`
> **Status:** `exact_within_reduced_bundle`
> **Layer:** `derivation`
> **Type:** `moving_throat_stage`

## Summary

Converts conservative operator moments into normalized grouped response and isolates P0=N0/D0.

## Physical Meaning

Converts conservative operator moments into normalized grouped response and isolates P0=N0/D0.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `exact_within_reduced_bundle`

## Equation

$$
u2=-D2/D0
$$

$$
P0=N0/D0
$$

$$
mhat0²P0=54Gc_s^5/(5a^5c^5)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_GROUPED_RESPONSE_MOMENTS]]

### Related claims
- [[CLAIM_STAGE022_P0_TARGET]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GENERALIZES_TO` | [[MT_STAGE023_FULL_GROUPED_BUNDLE]] | Full 20/21/22 bundle and projectors. |
| `ISOLATES` | [[OPEN_QUAD_NORMALIZATION]] | Invariant product mhat0²P0 is the remaining normalization target. |
| `SUPPORTS` | [[PN_2_5_QUAD_NARROWING]] | Moving-throat bridge supplies exact reduced target mhat0²P0. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_GROUPED_RESPONSE_MOMENTS]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[MT_STAGE021_REDUCED_ONE_PORT]] | One-lane outgoing transfer lifted to grouped P2 normalization language. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE022_P0_TARGET]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/pde_ledger/paper/stages/stage_022.tex`
- `notes/moving_throat_pde_program_compact.md`
- `moving_throat_pde_stage022_grouped_p2_normalization_bridge.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
