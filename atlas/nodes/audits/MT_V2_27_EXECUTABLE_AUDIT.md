---
id: MT_V2_27_EXECUTABLE_AUDIT
title: V2-27 executable audit implementation
type: reproducibility_gate
layer: status_audit
status: referee_harness_pass
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: 'Records executable harness: audit scripts, fixture checks, simulations, Mathematica mirrors, and reproducibility status.'
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
outgoing_edges:
- target: MT_V2_26_STATUS
  relation: SUPPORTS
  status: reproducibility
  note: Executable harness supports audit claims but does not prove branch existence.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/referee_harness_pass
- topic/moving_throat
- type/reproducibility_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-27 executable audit implementation

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_27_EXECUTABLE_AUDIT`
> **Status:** `referee_harness_pass`
> **Layer:** `status_audit`
> **Type:** `reproducibility_gate`

## Summary

Records executable harness: audit scripts, fixture checks, simulations, Mathematica mirrors, and reproducibility status.

## Physical Meaning

Records executable harness: audit scripts, fixture checks, simulations, Mathematica mirrors, and reproducibility status.

## Mathematical Role

- Layer: `status_audit`
- Type: `reproducibility_gate`
- Status: `referee_harness_pass`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

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
| `SUPPORTS` | [[MT_V2_26_STATUS]] | Executable harness supports audit claims but does not prove branch existence. |

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
