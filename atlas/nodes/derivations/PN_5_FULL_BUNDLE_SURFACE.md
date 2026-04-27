---
id: PN_5_FULL_BUNDLE_SURFACE
title: 5PN / full-bundle target surface
type: PN_bridge
layer: derivation
status: reduced_target_surface_open_realization
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Later reduced stack sharpens isotropic full-bundle target and weak-axisymmetric prefactor slope/orbit-lock conditions.
future_paper_needed: false
source_files:
- notes/5pn/5pn_notes_full.md
- notes/pde_audit_full.md
- 5pn_notes_full.md
- pde_audit_full.md
legacy_sources:
- 5pn_notes_full.md
- pde_audit_full.md:V2-19
source_links:
- '[[FILE_5PN_FULL]]'
math_ids:
- MATH_MONOMIAL_QUOTIENT
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
open_gate_ids:
- OPEN_N2_N4_MOMENT_SHAPE
source_ids:
- FILE_5PN_FULL
outgoing_edges:
- target: OPEN_N2_N4_MOMENT_SHAPE
  relation: REQUIRES
  status: open
  note: Higher-order full-bundle surface needs outgoing moment-shape control.
- target: MATH_MONOMIAL_QUOTIENT
  relation: USES
  status: reduced
  note: 5PN weak-axisymmetric transport uses monomial quotient/orbit-lock machinery.
incoming_edges:
- source: PN_4_LOCAL_TAIL
  relation: CONTINUES_TO
  status: reduced
  note: Later full-bundle target surface and weak-axisymmetric packet.
- source: FILE_5PN_FULL
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PN_4_LOCAL_TAIL
  relation: FEEDS
  status: derivation
  note: 4PN tail interface motivates 5PN/full-bundle continuation.
- source: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle_open_branch
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_STAGE6_FULL_GROUPED_BUNDLE
  relation: SUPPORTS
  status: reduced
  note: Full bundle target surface matches later 5PN/orbit-lock compression.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/reduced_target_surface_open_realization
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/pn_bridge
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 5PN / full-bundle target surface

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_5_FULL_BUNDLE_SURFACE`  
> **Status:** `reduced_target_surface_open_realization`  
> **Layer:** `derivation`  
> **Type:** `PN_bridge`

## Summary

Later reduced stack sharpens isotropic full-bundle target and weak-axisymmetric prefactor slope/orbit-lock conditions.

## Physical Meaning

Later reduced stack sharpens isotropic full-bundle target and weak-axisymmetric prefactor slope/orbit-lock conditions.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_bridge`
- Status: `reduced_target_surface_open_realization`

## Equation

$$
D0*C/(3A²)=1
$$

$$
P2=P4=0
$$

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
- none

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]

### Open gates
- [[OPEN_N2_N4_MOMENT_SHAPE]]

### Status firewalls
- none

### Source anchors
- [[FILE_5PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REQUIRES` | [[OPEN_N2_N4_MOMENT_SHAPE]] | Higher-order full-bundle surface needs outgoing moment-shape control. |
| `USES` | [[MATH_MONOMIAL_QUOTIENT]] | 5PN weak-axisymmetric transport uses monomial quotient/orbit-lock machinery. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTINUES_TO` | [[PN_4_LOCAL_TAIL]] | Later full-bundle target surface and weak-axisymmetric packet. |
| `DOCUMENTS` | [[FILE_5PN_FULL]] | File anchor documents this node. |
| `FEEDS` | [[PN_4_LOCAL_TAIL]] | 4PN tail interface motivates 5PN/full-bundle continuation. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Claim feeds this downstream object, output, or open gate. |
| `SUPPORTS` | [[MT_STAGE6_FULL_GROUPED_BUNDLE]] | Full bundle target surface matches later 5PN/orbit-lock compression. |

## Source Anchors

### Source anchor notes
- [[FILE_5PN_FULL]]

### Source files
- `notes/5pn/5pn_notes_full.md`
- `notes/pde_audit_full.md`
- `5pn_notes_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
