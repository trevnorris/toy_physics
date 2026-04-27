---
id: MT_V2_08_BDG_SCHUR
title: V2-08 BdG-wall Schur complement audit
type: audit_gate
layer: status_audit
status: passed_with_stability_gates
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Audits stable-mode Schur complement, positivity/softening bounds, pole stability, and grouped P2 isotropy.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
math_ids:
- MATH_BDG_SCHUR_COMPLEMENT
outgoing_edges:
- target: MATH_BDG_SCHUR_COMPLEMENT
  relation: VALIDATES
  status: reduced
  note: Stable BdG modes produce Schur-complement support moments.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/passed_with_stability_gates
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-08 BdG-wall Schur complement audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_08_BDG_SCHUR`  
> **Status:** `passed_with_stability_gates`  
> **Layer:** `status_audit`  
> **Type:** `audit_gate`

## Summary

Audits stable-mode Schur complement, positivity/softening bounds, pole stability, and grouped P2 isotropy.

## Physical Meaning

Audits stable-mode Schur complement, positivity/softening bounds, pole stability, and grouped P2 isotropy.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `passed_with_stability_gates`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_BDG_SCHUR_COMPLEMENT]]

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
| `VALIDATES` | [[MATH_BDG_SCHUR_COMPLEMENT]] | Stable BdG modes produce Schur-complement support moments. |

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
