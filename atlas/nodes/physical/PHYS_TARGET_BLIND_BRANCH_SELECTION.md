---
id: PHYS_TARGET_BLIND_BRANCH_SELECTION
title: Target-blind branch selection
type: protocol_principle
layer: physical_ontology
status: mandatory
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The actual branch must be frozen before target evaluation; post-hoc fitting to GR/PN residuals is prohibited by the audit.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
source_links:
- '[[FILE_PDE_AUDIT]]'
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
open_gate_ids:
- OPEN_ACTUAL_BRANCH_EXPORTER
source_ids:
- FILE_PDE_AUDIT
outgoing_edges:
- target: OPEN_ACTUAL_BRANCH_EXPORTER
  relation: GOVERNS
  status: mandatory
  note: Exporter must be frozen before comparing to target packets.
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_actual_branch
  note: Physical ontology object grounded by this claim.
- target: MT_V2_16_NO_REFIT
  relation: IMPLEMENTS
  status: mandatory
  note: No-refit protocol enforces target-blind branch freezing.
incoming_edges:
- source: FILE_PDE_AUDIT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/mandatory
- topic/moving_throat
- type/protocol_principle
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Target-blind branch selection

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_TARGET_BLIND_BRANCH_SELECTION`
> **Status:** `mandatory`
> **Layer:** `physical_ontology`
> **Type:** `protocol_principle`

## Summary

The actual branch must be frozen before target evaluation; post-hoc fitting to GR/PN residuals is prohibited by the audit.

## Physical Meaning

The actual branch must be frozen before target evaluation; post-hoc fitting to GR/PN residuals is prohibited by the audit.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `protocol_principle`
- Status: `mandatory`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]

### Open gates
- [[OPEN_ACTUAL_BRANCH_EXPORTER]]

### Status firewalls
- none

### Source anchors
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GOVERNS` | [[OPEN_ACTUAL_BRANCH_EXPORTER]] | Exporter must be frozen before comparing to target packets. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Physical ontology object grounded by this claim. |
| `IMPLEMENTS` | [[MT_V2_16_NO_REFIT]] | No-refit protocol enforces target-blind branch freezing. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_PDE_AUDIT]] | File anchor documents this node. |

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
