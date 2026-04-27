---
id: MT_STAGE4_MAXWELL_MIXED
title: Stage 4 Maxwell/mixed outgoing bridge
type: moving_throat_stage
layer: derivation
status: reduced_controlled
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Adds localized Maxwell/mixed U/W block and passive outgoing port transfer.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- notes/pde_audit_full.md
- moving_throat_output_full.md
- pde_audit_full.md
legacy_sources:
- moving_throat_output_full.md
- pde_audit_full.md:V2-09
math_ids:
- MATH_MAXWELL_MIXED_KERNEL
- MATH_MIXED_FIELDS_EW_CA
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
outgoing_edges:
- target: MT_STAGE5_GROUPED_P2_BRIDGE
  relation: FEEDS
  status: reduced
  note: One-lane outgoing transfer lifted to grouped P2 normalization language.
incoming_edges:
- source: MATH_MIXED_FIELDS_EW_CA
  relation: ENABLES
  status: exact/reduced
  note: Mixed fields are the microscopic place for outgoing bridge.
- source: MT_STAGE3_BDG_COUPLING
  relation: FEEDS
  status: reduced
  note: Conservative support self-energy enlarged by localized Maxwell/mixed block.
- source: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
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
- status/reduced_controlled
- topic/maxwell
- topic/moving_throat
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 4 Maxwell/mixed outgoing bridge

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE4_MAXWELL_MIXED`  
> **Status:** `reduced_controlled`  
> **Layer:** `derivation`  
> **Type:** `moving_throat_stage`

## Summary

Adds localized Maxwell/mixed U/W block and passive outgoing port transfer.

## Physical Meaning

Adds localized Maxwell/mixed U/W block and passive outgoing port transfer.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `reduced_controlled`

## Equation

$$
N_l(0)=[Ω_A²g_W+Rg_A]²/[Ω_A²Ω_W²-R²]²
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
- [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MT_STAGE5_GROUPED_P2_BRIDGE]] | One-lane outgoing transfer lifted to grouped P2 normalization language. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ENABLES` | [[MATH_MIXED_FIELDS_EW_CA]] | Mixed fields are the microscopic place for outgoing bridge. |
| `FEEDS` | [[MT_STAGE3_BDG_COUPLING]] | Conservative support self-energy enlarged by localized Maxwell/mixed block. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]] | Claim feeds this downstream object, output, or open gate. |
| `FORMALIZES` | [[MATH_MAXWELL_MIXED_KERNEL]] | Matches Stage-4 port-active mixed-sector bridge. |

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
