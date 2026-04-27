---
id: MATH_PARENT_ACTION_PROMOTED
title: Promoted parent action with S_Sigma
type: action_target
layer: math_object
status: open_patch_required
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Parent-complete formulation needed for autonomous moving-throat PDE claim.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-01
math_ids:
- MATH_WALL_ACTION_S_ETA
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
open_gate_ids:
- OPEN_PARENT_PROMOTION_S_SIGMA
incoming_edges:
- source: MATH_WALL_ACTION_S_ETA
  relation: APPROXIMATES
  status: effective-to-parent patch
  note: Quadratic wall action can be first approximation to promoted S_Sigma.
- source: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: FEEDS_OR_STATUS_OF
  status: strict_parent_fail_effective_wall_pass
  note: Claim feeds this downstream object, output, or open gate.
- source: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: REQUIRES
  status: open
  note: Parent-complete statement needs S_Sigma.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/open_patch_required
- topic/moving_throat
- type/action_target
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Promoted parent action with S_Sigma

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_PARENT_ACTION_PROMOTED`  
> **Status:** `open_patch_required`  
> **Layer:** `math_object`  
> **Type:** `action_target`

## Summary

Parent-complete formulation needed for autonomous moving-throat PDE claim.

## Physical Meaning

Parent-complete formulation needed for autonomous moving-throat PDE claim.

## Mathematical Role

- Layer: `math_object`
- Type: `action_target`
- Status: `open_patch_required`

## Equation

$$
S_total=S_psi+S_EM+S_Sigma
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WALL_ACTION_S_ETA]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- [[OPEN_PARENT_PROMOTION_S_SIGMA]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `APPROXIMATES` | [[MATH_WALL_ACTION_S_ETA]] | Quadratic wall action can be first approximation to promoted S_Sigma. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Claim feeds this downstream object, output, or open gate. |
| `REQUIRES` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Parent-complete statement needs S_Sigma. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
