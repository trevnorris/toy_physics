---
id: MATH_WALL_ACTION_S_ETA
title: Quadratic wall action S_eta^(2)
type: effective_action
layer: math_object
status: effective_linear_wall_closure_passed
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Consistent reduced wall action that generates scalar and P2 wall PDEs, but must be promoted for strict parent-level status.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_notes_full.md
- pde_audit_full.md
- moving_throat_output_full.md
legacy_sources:
- pde_audit_full.md:V2-01
- moving_throat_output_full.md
math_ids:
- MATH_PARENT_ACTION_PROMOTED
claim_ids:
- CLAIM_STAGE2_AL_RECOVERY
status_firewall_ids:
- FIREWALL_PARENT_WALL_NOT_STRICT
- FIREWALL_WALL_COEFFS_BRANCH_DATA
outgoing_edges:
- target: MATH_PARENT_ACTION_PROMOTED
  relation: APPROXIMATES
  status: effective-to-parent patch
  note: Quadratic wall action can be first approximation to promoted S_Sigma.
- target: MT_STAGE2_BREATHING_REDUCTION
  relation: REDUCES_TO
  status: exact within closure
  note: Axisymmetric two-mode truncation recovers old a,L closure.
incoming_edges:
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_WALL_ACTION_S_ETA.
- source: CLAIM_STAGE2_AL_RECOVERY
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_wall_action
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_STAGE1_GEOMETRY_LIFT
  relation: INTRODUCES
  status: effective
  note: Stage 1 introduces minimal quadratic wall action as new ansatz.
- source: FIREWALL_PARENT_WALL_NOT_STRICT
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_WALL_COEFFS_BRANCH_DATA
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_WALL_COEFFS_ARE_FIT_KNOBS
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MATH_WALL_ACTION_S_ETA.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/effective_linear_wall_closure_passed
- topic/moving_throat
- type/effective_action
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Quadratic wall action S_eta^(2)

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_WALL_ACTION_S_ETA`
> **Status:** `effective_linear_wall_closure_passed`
> **Layer:** `math_object`
> **Type:** `effective_action`

## Summary

Consistent reduced wall action that generates scalar and P2 wall PDEs, but must be promoted for strict parent-level status.

## Physical Meaning

Consistent reduced wall action that generates scalar and P2 wall PDEs, but must be promoted for strict parent-level status.

## Mathematical Role

- Layer: `math_object`
- Type: `effective_action`
- Status: `effective_linear_wall_closure_passed`

## Equation

$$
S_eta^(2)=1/2∫(mu eta_t^2-T_w eta_w^2-T_Omega eta(-Delta_S2)eta-K_eta eta^2)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_PARENT_ACTION_PROMOTED]]

### Related equations
- none

### Related claims
- [[CLAIM_STAGE2_AL_RECOVERY]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_PARENT_WALL_NOT_STRICT]]
- [[FIREWALL_WALL_COEFFS_BRANCH_DATA]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `APPROXIMATES` | [[MATH_PARENT_ACTION_PROMOTED]] | Quadratic wall action can be first approximation to promoted S_Sigma. |
| `REDUCES_TO` | [[MT_STAGE2_BREATHING_REDUCTION]] | Axisymmetric two-mode truncation recovers old a,L closure. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references MATH_WALL_ACTION_S_ETA. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE2_AL_RECOVERY]] | Claim feeds this downstream object, output, or open gate. |
| `INTRODUCES` | [[MT_STAGE1_GEOMETRY_LIFT]] | Stage 1 introduces minimal quadratic wall action as new ansatz. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_PARENT_WALL_NOT_STRICT]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_WALL_COEFFS_BRANCH_DATA]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_WALL_COEFFS_ARE_FIT_KNOBS]] | Negative query starts from MATH_WALL_ACTION_S_ETA. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_notes_full.md`
- `pde_audit_full.md`
- `moving_throat_output_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
