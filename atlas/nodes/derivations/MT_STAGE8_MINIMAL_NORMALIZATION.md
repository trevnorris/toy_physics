---
id: MT_STAGE8_MINIMAL_NORMALIZATION
title: Stage 8 minimal isotropic normalization
type: moving_throat_stage
layer: derivation
status: reduced_model
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Collapses radial/axial normalization to one explicit minimal scalar expression.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- moving_throat_output_full.md
legacy_sources:
- moving_throat_output_full.md
claim_ids:
- CLAIM_STAGE8_TO_14_SELECTED_BRANCH
incoming_edges:
- source: CLAIM_STAGE8_TO_14_SELECTED_BRANCH
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_selected_reduced_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_STAGE7_OVERLAP_ISOTROPY
  relation: NARROWS_TO
  status: reduced
  note: Angular side closed; radial/axial normalization remains.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/reduced_model
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 8 minimal isotropic normalization

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE8_MINIMAL_NORMALIZATION`  
> **Status:** `reduced_model`  
> **Layer:** `derivation`  
> **Type:** `moving_throat_stage`

## Summary

Collapses radial/axial normalization to one explicit minimal scalar expression.

## Physical Meaning

Collapses radial/axial normalization to one explicit minimal scalar expression.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `reduced_model`

## Equation

$$
mhat_rad²P²/[Delta(KDelta-Delta C²/varpi²-Q)] = target
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
- [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]] | Claim feeds this downstream object, output, or open gate. |
| `NARROWS_TO` | [[MT_STAGE7_OVERLAP_ISOTROPY]] | Angular side closed; radial/axial normalization remains. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/moving_throat_notes_full.md`
- `moving_throat_output_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
