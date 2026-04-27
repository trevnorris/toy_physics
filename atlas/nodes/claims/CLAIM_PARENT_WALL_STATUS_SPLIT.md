---
id: CLAIM_PARENT_WALL_STATUS_SPLIT
title: 'Parent wall status split: force yes, autonomous PDE no'
type: claim
layer: claim_theorem
status: strict_parent_fail_effective_wall_pass
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The confinement-only parent supplies a wall force/source but no autonomous wall PDE; S_eta^(2) is a consistent linear effective closure and must be promoted to S_Sigma for paren...
future_paper_needed: false
source_links:
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_PARENT_WALL_EXEC]]'
- '[[SEC_PDE_WALL_VARIATION]]'
physical_ids:
- PHYS_FINITE_MOUTH_SHAPE
- PHYS_RESPONSE_READOUTS
math_ids:
- MATH_PARENT_ACTION_PROMOTED
equation_ids:
- EQ_PARENT_ACTION_CURRENT
- EQ_WALL_MODAL_PDE
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
- CLAIM_STAGE1_GEOMETRY_LIFT
open_gate_ids:
- OPEN_PARENT_PROMOTION_S_SIGMA
status_firewall_ids:
- FIREWALL_PARENT_WALL_NOT_STRICT
- FIREWALL_WALL_COEFFS_BRANCH_DATA
source_ids:
- FILE_PDE_AUDIT
- SEC_PDE_PARENT_WALL_EXEC
- SEC_PDE_WALL_VARIATION
outgoing_edges:
- target: MATH_PARENT_ACTION_PROMOTED
  relation: FEEDS_OR_STATUS_OF
  status: strict_parent_fail_effective_wall_pass
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_V2_01_PARENT_WALL_AUDIT
  relation: FEEDS_OR_STATUS_OF
  status: strict_parent_fail_effective_wall_pass
  note: Claim feeds this downstream object, output, or open gate.
- target: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: FEEDS_OR_STATUS_OF
  status: strict_parent_fail_effective_wall_pass
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: OPENS_REQUIRED_WORK
  status: active
  note: Claim-level dependency added in v0.4.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: STATUS_LIMITS
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_PDE_PARENT_WALL_EXEC
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent wall status split.
- source: SEC_PDE_WALL_VARIATION
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: S_eta gives wall PDE.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_PARENT_WALL_STATUS_SPLIT.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: strict_parent_fail_effective_wall_pass
  note: 'Claim class: patched_required'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: GROUNDS_PHYSICAL_MEANING
  status: strict_parent_fail_effective_wall_pass
  note: Physical ontology object grounded by this claim.
- source: PHYS_RESPONSE_READOUTS
  relation: GROUNDS_PHYSICAL_MEANING
  status: strict_parent_fail_effective_wall_pass
  note: Physical ontology object grounded by this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: strict_parent_fail_effective_wall_pass
  note: Source artifact anchors this claim.
- source: FIREWALL_PARENT_WALL_NOT_STRICT
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_WALL_COEFFS_BRANCH_DATA
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: QUERY_PARENT_WALL_STATUS
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
- source: EQ_PARENT_ACTION_CURRENT
  relation: SUPPORTS_CLAIM
  status: strict_parent_fail_effective_wall_pass
  note: Equation anchor supports this named claim.
- source: EQ_WALL_MODAL_PDE
  relation: SUPPORTS_CLAIM
  status: strict_parent_fail_effective_wall_pass
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/strict_parent_fail_effective_wall_pass
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Parent wall status split: force yes, autonomous PDE no

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_PARENT_WALL_STATUS_SPLIT`  
> **Status:** `strict_parent_fail_effective_wall_pass`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

The confinement-only parent supplies a wall force/source but no autonomous wall PDE; S_eta^(2) is a consistent linear effective closure and must be promoted to S_Sigma for parent-level dynamics.

## Claim

The confinement-only parent supplies a wall force/source but no autonomous wall PDE; S_eta^(2) is a consistent linear effective closure and must be promoted to S_Sigma for parent-level dynamics.

## Physical Meaning

The confinement-only parent supplies a wall force/source but no autonomous wall PDE; S_eta^(2) is a consistent linear effective closure and must be promoted to S_Sigma for parent-level dynamics.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `strict_parent_fail_effective_wall_pass`
- Outputs: `MT_V2_01_PARENT_WALL_AUDIT`, `OPEN_PARENT_PROMOTION_S_SIGMA`, `MATH_PARENT_ACTION_PROMOTED`

## Atlas Links

### Related physical nodes
- [[PHYS_FINITE_MOUTH_SHAPE]]
- [[PHYS_RESPONSE_READOUTS]]

### Related math nodes
- [[MATH_PARENT_ACTION_PROMOTED]]

### Related equations
- [[EQ_PARENT_ACTION_CURRENT]]
- [[EQ_WALL_MODAL_PDE]]

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- [[OPEN_PARENT_PROMOTION_S_SIGMA]]

### Status firewalls
- [[FIREWALL_PARENT_WALL_NOT_STRICT]]
- [[FIREWALL_WALL_COEFFS_BRANCH_DATA]]

### Source anchors
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_PARENT_WALL_EXEC]]
- [[SEC_PDE_WALL_VARIATION]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_PARENT_ACTION_PROMOTED]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_V2_01_PARENT_WALL_AUDIT]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Claim feeds this downstream object, output, or open gate. |
| `OPENS_REQUIRED_WORK` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Claim-level dependency added in v0.4. |
| `STATUS_LIMITS` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_PARENT_WALL_EXEC]] | Parent wall status split. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_WALL_VARIATION]] | S_eta gives wall PDE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_PARENT_WALL_STATUS_SPLIT. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: patched_required |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_FINITE_MOUTH_SHAPE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_RESPONSE_READOUTS]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_PARENT_WALL_NOT_STRICT]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_WALL_COEFFS_BRANCH_DATA]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[QUERY_PARENT_WALL_STATUS]] | Query validation start node. |
| `SUPPORTS_CLAIM` | [[EQ_PARENT_ACTION_CURRENT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_WALL_MODAL_PDE]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_PARENT_WALL_EXEC]]
- [[SEC_PDE_WALL_VARIATION]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
