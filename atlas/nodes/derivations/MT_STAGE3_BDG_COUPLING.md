---
id: MT_STAGE3_BDG_COUPLING
title: Stage 3 BdG-wall coupling
type: moving_throat_stage
layer: derivation
status: reduced_controlled
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Integrates stable BdG support modes out by Schur complement, giving conservative self-energies and pole shifts.
future_paper_needed: false
source_files:
- notes/moving_throat_notes_full.md
- notes/pde_audit_full.md
- moving_throat_output_full.md
- pde_audit_full.md
legacy_sources:
- moving_throat_output_full.md
- pde_audit_full.md:V2-08
math_ids:
- MATH_BDG_SCHUR_COMPLEMENT
claim_ids:
- CLAIM_STAGE3_BDG_SCHUR
outgoing_edges:
- target: MT_STAGE4_MAXWELL_MIXED
  relation: FEEDS
  status: reduced
  note: Conservative support self-energy enlarged by localized Maxwell/mixed block.
incoming_edges:
- source: MT_STAGE2_BREATHING_REDUCTION
  relation: FEEDS
  status: reduced
  note: Wall modes then coupled to stable matter support modes.
- source: CLAIM_STAGE3_BDG_SCHUR
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_stable_modes
  note: Claim feeds this downstream object, output, or open gate.
- source: MATH_BDG_SCHUR_COMPLEMENT
  relation: FORMALIZES
  status: reduced
  note: Matches Stage-3 wall/BdG coupling kernel.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/reduced_controlled
- topic/moving_throat
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 3 BdG-wall coupling

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE3_BDG_COUPLING`  
> **Status:** `reduced_controlled`  
> **Layer:** `derivation`  
> **Type:** `moving_throat_stage`

## Summary

Integrates stable BdG support modes out by Schur complement, giving conservative self-energies and pole shifts.

## Physical Meaning

Integrates stable BdG support modes out by Schur complement, giving conservative self-energies and pole shifts.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `reduced_controlled`

## Equation

$$
D_A^eff(ω)=K_A-M_Aω²-Σ g²/(varpi²-ω²)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_BDG_SCHUR_COMPLEMENT]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE3_BDG_SCHUR]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MT_STAGE4_MAXWELL_MIXED]] | Conservative support self-energy enlarged by localized Maxwell/mixed block. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MT_STAGE2_BREATHING_REDUCTION]] | Wall modes then coupled to stable matter support modes. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE3_BDG_SCHUR]] | Claim feeds this downstream object, output, or open gate. |
| `FORMALIZES` | [[MATH_BDG_SCHUR_COMPLEMENT]] | Matches Stage-3 wall/BdG coupling kernel. |

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
