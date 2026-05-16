---
id: OPEN_PARENT_PROMOTION_S_SIGMA
title: 'Open gate: promote S_Sigma'
type: parent_status_gate
layer: open_gate
status: open_patch_required
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Needed before saying moving throat is autonomous parent-level dynamical field.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-01
source_links:
- '[[SEC_PDE_PARENT_WALL_EXEC]]'
math_ids:
- MATH_PARENT_ACTION_PROMOTED
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
open_gate_ids:
- OPEN_NONLINEAR_S_SIGMA
status_firewall_ids:
- FIREWALL_PARENT_WALL_NOT_STRICT
source_ids:
- SEC_PDE_PARENT_WALL_EXEC
outgoing_edges:
- target: OPEN_NONLINEAR_S_SIGMA
  relation: REFINES_TO
  status: open
  note: Best patch is a nonlinear S_Sigma whose quadratic limit is S_eta^(2).
- target: MATH_PARENT_ACTION_PROMOTED
  relation: REQUIRES
  status: open
  note: Parent-complete statement needs S_Sigma.
incoming_edges:
- source: SEC_PDE_PARENT_WALL_EXEC
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent wall status split.
- source: QUERY_PARENT_WALL_STATUS
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- source: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: FEEDS_OR_STATUS_OF
  status: strict_parent_fail_effective_wall_pass
  note: Claim feeds this downstream object, output, or open gate.
- source: BACKLINK_PDE_AUDIT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_PARENT_PROMOTION_S_SIGMA.
- source: ATLAS_CURRENT_READINESS_V05
  relation: FLAGS_REMAINING_GATE
  status: v05
  note: Still open after atlas organization; atlas tracks it but does not solve it.
- source: ATLAS_CURRENT_READINESS_V06
  relation: FLAGS_REMAINING_GATE
  status: v06
  note: Still open after v0.6 organization pass.
- source: MT_V2_01_PARENT_WALL_AUDIT
  relation: OPENS_GATE
  status: patch_required
  note: Strict parent-level dynamic wall claim fails unless S_eta/S_Sigma is promoted.
- source: MT_V2_01_PARENT_WALL_AUDIT
  relation: OPENS_GATE
  status: open patch
  note: Strict parent-level wall dynamics fail unless S_eta/S_Sigma is added.
- source: FIREWALL_PARENT_WALL_NOT_STRICT
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open_patch_required
- topic/moving_throat
- type/parent_status_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Open gate: promote S_Sigma

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_PARENT_PROMOTION_S_SIGMA`
> **Status:** `open_patch_required`
> **Layer:** `open_gate`
> **Type:** `parent_status_gate`

## Summary

Needed before saying moving throat is autonomous parent-level dynamical field.

## What Remains Open

Needed before saying moving throat is autonomous parent-level dynamical field.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Needed before saying moving throat is autonomous parent-level dynamical field.

## Mathematical Role

- Layer: `open_gate`
- Type: `parent_status_gate`
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
- [[MATH_PARENT_ACTION_PROMOTED]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- [[OPEN_NONLINEAR_S_SIGMA]]

### Status firewalls
- [[FIREWALL_PARENT_WALL_NOT_STRICT]]

### Source anchors
- [[SEC_PDE_PARENT_WALL_EXEC]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REFINES_TO` | [[OPEN_NONLINEAR_S_SIGMA]] | Best patch is a nonlinear S_Sigma whose quadratic limit is S_eta^(2). |
| `REQUIRES` | [[MATH_PARENT_ACTION_PROMOTED]] | Parent-complete statement needs S_Sigma. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_PARENT_WALL_EXEC]] | Parent wall status split. |
| `EXPECTS_TARGET` | [[QUERY_PARENT_WALL_STATUS]] | Query validation expected target node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Claim feeds this downstream object, output, or open gate. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block flags open gate OPEN_PARENT_PROMOTION_S_SIGMA. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V05]] | Still open after atlas organization; atlas tracks it but does not solve it. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V06]] | Still open after v0.6 organization pass. |
| `OPENS_GATE` | [[MT_V2_01_PARENT_WALL_AUDIT]] | Strict parent-level dynamic wall claim fails unless S_eta/S_Sigma is promoted. |
| `OPENS_GATE` | [[MT_V2_01_PARENT_WALL_AUDIT]] | Strict parent-level wall dynamics fail unless S_eta/S_Sigma is added. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_PARENT_WALL_NOT_STRICT]] | Firewall preserves this correct status boundary. |

## Source Anchors

### Source anchor notes
- [[SEC_PDE_PARENT_WALL_EXEC]]

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
