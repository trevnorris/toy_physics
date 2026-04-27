---
id: MT_V2_01_PARENT_WALL_AUDIT
title: V2-01 parent/wall action audit
type: audit_gate
layer: status_audit
status: strict_parent_fail_effective_closure_pass
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Confinement-only parent gives wall force, not autonomous wall PDE; S_eta/S_Sigma needed for parent dynamics.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-01
source_links:
- '[[FILE_PDE_AUDIT]]'
math_ids:
- MATH_PARENT_ACTION_CURRENT
- MATH_WALL_MODAL_PDE
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
open_gate_ids:
- OPEN_PARENT_PROMOTION_S_SIGMA
source_ids:
- FILE_PDE_AUDIT
outgoing_edges:
- target: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: OPENS_GATE
  status: patch_required
  note: Strict parent-level dynamic wall claim fails unless S_eta/S_Sigma is promoted.
- target: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: OPENS_GATE
  status: open patch
  note: Strict parent-level wall dynamics fail unless S_eta/S_Sigma is added.
- target: MATH_WALL_MODAL_PDE
  relation: VALIDATES_IF_INCLUDED
  status: effective_closure_pass
  note: Quadratic wall action supplies modal wall PDE.
incoming_edges:
- source: MATH_PARENT_ACTION_CURRENT
  relation: AUDITED_BY
  status: audit
  note: Audit checks whether Sigma/R is autonomous dynamical field.
- source: FILE_PDE_AUDIT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: FEEDS_OR_STATUS_OF
  status: strict_parent_fail_effective_wall_pass
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/strict_parent_fail_effective_closure_pass
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-01 parent/wall action audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_01_PARENT_WALL_AUDIT`  
> **Status:** `strict_parent_fail_effective_closure_pass`  
> **Layer:** `status_audit`  
> **Type:** `audit_gate`

## Summary

Confinement-only parent gives wall force, not autonomous wall PDE; S_eta/S_Sigma needed for parent dynamics.

## Physical Meaning

Confinement-only parent gives wall force, not autonomous wall PDE; S_eta/S_Sigma needed for parent dynamics.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `strict_parent_fail_effective_closure_pass`

## Equation

$$
STRICT_PARENT_DYNAMIC_WALL: FAIL unless S_eta/S_Sigma included
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_PARENT_ACTION_CURRENT]]
- [[MATH_WALL_MODAL_PDE]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- [[OPEN_PARENT_PROMOTION_S_SIGMA]]

### Status firewalls
- none

### Source anchors
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `OPENS_GATE` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Strict parent-level dynamic wall claim fails unless S_eta/S_Sigma is promoted. |
| `OPENS_GATE` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Strict parent-level wall dynamics fail unless S_eta/S_Sigma is added. |
| `VALIDATES_IF_INCLUDED` | [[MATH_WALL_MODAL_PDE]] | Quadratic wall action supplies modal wall PDE. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `AUDITED_BY` | [[MATH_PARENT_ACTION_CURRENT]] | Audit checks whether Sigma/R is autonomous dynamical field. |
| `DOCUMENTS` | [[FILE_PDE_AUDIT]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
