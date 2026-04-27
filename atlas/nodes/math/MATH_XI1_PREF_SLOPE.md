---
id: MATH_XI1_PREF_SLOPE
title: Weak-axisymmetric prefactor slope Xi1
type: prefactor_slope
layer: math_object
status: exact_within_reduced_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The weak-axisymmetric normalization defect equals the outgoing-prefactor slope P1/P0 and raw transfer-shape slope.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
math_ids:
- MATH_MONOMIAL_QUOTIENT
equation_ids:
- EQ_XI1_PREF_SLOPE
open_gate_ids:
- OPEN_WEAK_AXISYM_ORBIT_LOCK
status_firewall_ids:
- FIREWALL_G2_COMMON_CONDITIONAL
outgoing_edges:
- target: OPEN_WEAK_AXISYM_ORBIT_LOCK
  relation: FEEDS
  status: open
  note: Actual branch must export Xi1 or show orbit-lock.
incoming_edges:
- source: EQ_XI1_PREF_SLOPE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_5PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_XI1_PREF_SLOPE.
- source: BACKLINK_G2_OUTPUT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_XI1_PREF_SLOPE.
- source: QUERY_G2_XI1
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- source: ANOMALY_G2_COMMON_QUOTIENT
  relation: IMPORTS
  status: conditional
  note: g-2 residual can be expressed as transfer-shape / outgoing-prefactor slope.
- source: MATH_MONOMIAL_QUOTIENT
  relation: ORGANIZES
  status: exact_within_branch
  note: Xi1 is a direct defect coordinate/prefactor slope within the quotient packet.
- source: FIREWALL_G2_COMMON_CONDITIONAL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_within_reduced_branch
- topic/moving_throat
- topic/quadrupole
- type/prefactor_slope
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weak-axisymmetric prefactor slope Xi1

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_XI1_PREF_SLOPE`  
> **Status:** `exact_within_reduced_branch`  
> **Layer:** `math_object`  
> **Type:** `prefactor_slope`

## Summary

The weak-axisymmetric normalization defect equals the outgoing-prefactor slope P1/P0 and raw transfer-shape slope.

## Physical Meaning

The weak-axisymmetric normalization defect equals the outgoing-prefactor slope P1/P0 and raw transfer-shape slope.

## Mathematical Role

- Layer: `math_object`
- Type: `prefactor_slope`
- Status: `exact_within_reduced_branch`

## Equation

$$
Xi1=P1/P0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MONOMIAL_QUOTIENT]]

### Related equations
- [[EQ_XI1_PREF_SLOPE]]

### Related claims
- none

### Open gates
- [[OPEN_WEAK_AXISYM_ORBIT_LOCK]]

### Status firewalls
- [[FIREWALL_G2_COMMON_CONDITIONAL]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[OPEN_WEAK_AXISYM_ORBIT_LOCK]] | Actual branch must export Xi1 or show orbit-lock. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_XI1_PREF_SLOPE]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_5PN_FULL]] | Paper backlink block references MATH_XI1_PREF_SLOPE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_G2_OUTPUT]] | Paper backlink block references MATH_XI1_PREF_SLOPE. |
| `EXPECTS_TARGET` | [[QUERY_G2_XI1]] | Query validation expected target node. |
| `IMPORTS` | [[ANOMALY_G2_COMMON_QUOTIENT]] | g-2 residual can be expressed as transfer-shape / outgoing-prefactor slope. |
| `ORGANIZES` | [[MATH_MONOMIAL_QUOTIENT]] | Xi1 is a direct defect coordinate/prefactor slope within the quotient packet. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_G2_COMMON_CONDITIONAL]] | Firewall preserves this correct status boundary. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_pde_program_compact.md`
- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
