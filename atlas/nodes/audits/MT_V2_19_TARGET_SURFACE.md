---
id: MT_V2_19_TARGET_SURFACE
title: V2-19 isotropic full-bundle target surface
type: target_gate
layer: status_audit
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: 'Freezes full isotropic target packet: one-pole response, quadrupole normalization, constant-prefactor branch, tail transport.'
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
  relation: FREEZES
  status: exact_within_bundle
  note: Freezes the target packet for actual-branch export.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/exact_within_reduced_bundle
- topic/moving_throat
- topic/quadrupole
- type/target_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-19 isotropic full-bundle target surface

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_19_TARGET_SURFACE`  
> **Status:** `exact_within_reduced_bundle`  
> **Layer:** `status_audit`  
> **Type:** `target_gate`

## Summary

Freezes full isotropic target packet: one-pole response, quadrupole normalization, constant-prefactor branch, tail transport.

## Physical Meaning

Freezes full isotropic target packet: one-pole response, quadrupole normalization, constant-prefactor branch, tail transport.

## Mathematical Role

- Layer: `status_audit`
- Type: `target_gate`
- Status: `exact_within_reduced_bundle`

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
| `FREEZES` | [[MATH_FULL_BUNDLE_TARGET_SURFACE]] | Freezes the target packet for actual-branch export. |

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
