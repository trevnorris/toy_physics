---
id: CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
title: Fluxoid alone gives no universal 3D radial force sign
type: claim
layer: claim_theorem
status: exact_negative_within_3d_closure_analysis
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The circulation package shows that n fixes holonomy but not the current/plumbing coefficient, orientation sign, mutual-inductance closure, or ensemble, so no universal attractio...
future_paper_needed: false
source_links:
- '[[FILE_CIRCULATION_PACKAGE]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
- PHYS_MIXED_EM_CORE
math_ids:
- MATH_FLUXOID
claim_ids:
- CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
- CLAIM_MIXED_RECIRCULATION_OPEN
source_ids:
- FILE_CIRCULATION_PACKAGE
outgoing_edges:
- target: MATH_FLUXOID
  relation: LIMITS_SCOPE_OF
  status: exact_negative_within_3d_closure_analysis
  note: Fluxoid holonomy alone is not a 3D radial-force law.
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: PRESERVES_OPEN_GATE
  status: exact_negative_within_3d_closure_analysis
  note: Current/plumbing coefficients and mixed-sector transport remain external closure data.
incoming_edges:
- source: FILE_CIRCULATION_PACKAGE
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_negative_within_3d_closure_analysis
  note: Circulation package Step 02 anchors the no-universal-force scope result.
- source: CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
  relation: RESPECTS_SCOPE_LIMIT
  status: exact_within_fixed_current_orientation_closure
  note: Facing-mouth attraction statement remains closure-conditional and does not override the no-universal-force result.
- source: DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP
  relation: SUPPLIES_CLOSURE_EXAMPLE
  status: exact_within_fixed_current_far_field_closure
  note: Fixed-current loop law is a named closure, not a consequence of fluxoid alone.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_negative_within_3d_closure_analysis
- topic/maxwell
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Fluxoid alone gives no universal 3D radial force sign

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID`
> **Status:** `exact_negative_within_3d_closure_analysis`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The circulation package shows that n fixes holonomy but not the current/plumbing coefficient, orientation sign, mutual-inductance closure, or ensemble, so no universal attraction/repulsion sign follows from n1*n2 alone in 3D.

## Claim

The circulation package shows that n fixes holonomy but not the current/plumbing coefficient, orientation sign, mutual-inductance closure, or ensemble, so no universal attraction/repulsion sign follows from n1*n2 alone in 3D.

## Physical Meaning

The circulation package shows that n fixes holonomy but not the current/plumbing coefficient, orientation sign, mutual-inductance closure, or ensemble, so no universal attraction/repulsion sign follows from n1*n2 alone in 3D.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_negative_within_3d_closure_analysis`
- Outputs: `CLAIM_MIXED_RECIRCULATION_OPEN`, `OPEN_MIXED_RECIRCULATION`

## Atlas Links

### Related physical nodes
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- [[MATH_FLUXOID]]

### Related equations
- none

### Related claims
- [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_CIRCULATION_PACKAGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `LIMITS_SCOPE_OF` | [[MATH_FLUXOID]] | Fluxoid holonomy alone is not a 3D radial-force law. |
| `PRESERVES_OPEN_GATE` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Current/plumbing coefficients and mixed-sector transport remain external closure data. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_CIRCULATION_PACKAGE]] | Circulation package Step 02 anchors the no-universal-force scope result. |
| `RESPECTS_SCOPE_LIMIT` | [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]] | Facing-mouth attraction statement remains closure-conditional and does not override the no-universal-force ... |
| `SUPPLIES_CLOSURE_EXAMPLE` | [[DERIVATION_COAXIAL_NEUMANN_CURRENT_LOOP]] | Fixed-current loop law is a named closure, not a consequence of fluxoid alone. |

## Source Anchors

### Source anchor notes
- [[FILE_CIRCULATION_PACKAGE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
