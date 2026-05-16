---
id: MT_V2_13_GROUPED_NORMALIZATION
title: V2-13 grouped normalization ratio audit
type: audit_gate
layer: status_audit
status: exact_within_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Defines normalized response, outgoing prefactor, constant-prefactor branch, and isotropic full-bundle surface.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
math_ids:
- MATH_FULL_BUNDLE_TARGET_SURFACE
outgoing_edges:
- target: MATH_FULL_BUNDLE_TARGET_SURFACE
  relation: DEFINES
  status: exact_within_bundle
  note: Defines normalized response, P0, P2, P4, target surface.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/exact_within_bundle
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-13 grouped normalization ratio audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_13_GROUPED_NORMALIZATION`
> **Status:** `exact_within_bundle`
> **Layer:** `status_audit`
> **Type:** `audit_gate`

## Summary

Defines normalized response, outgoing prefactor, constant-prefactor branch, and isotropic full-bundle surface.

## Physical Meaning

Defines normalized response, outgoing prefactor, constant-prefactor branch, and isotropic full-bundle surface.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `exact_within_bundle`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_FULL_BUNDLE_TARGET_SURFACE]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DEFINES` | [[MATH_FULL_BUNDLE_TARGET_SURFACE]] | Defines normalized response, P0, P2, P4, target surface. |

## Incoming Edges

- none

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
