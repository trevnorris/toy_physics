---
id: SEC_PDE_PARENT_WALL_EXEC
title: Executive result
type: section_anchor
layer: source_section_anchor
status: v05_first_section_anchor
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Parent wall status split.
source_kind: future_paper_note
future_paper_needed: true
source_files:
- notes/pde_audit_full.md
source_links:
- '[[FILE_PDE_AUDIT]]'
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
open_gate_ids:
- OPEN_PARENT_PROMOTION_S_SIGMA
status_firewall_ids:
- FIREWALL_PARENT_WALL_NOT_STRICT
- FIREWALL_WALL_COEFFS_BRANCH_DATA
source_ids:
- FILE_PDE_AUDIT
outgoing_edges:
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent wall status split.
- target: OPEN_PARENT_PROMOTION_S_SIGMA
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent wall status split.
incoming_edges:
- source: FIREWALL_PARENT_WALL_NOT_STRICT
  relation: ANCHORED_IN
  status: v07
  note: Firewall is anchored in this source/file/section.
- source: FIREWALL_WALL_COEFFS_BRANCH_DATA
  relation: ANCHORED_IN
  status: v07
  note: Firewall is anchored in this source/file/section.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_SECTION_ANCHOR
  status: v06
  note: Paper backlink block references source-section anchor SEC_PDE_PARENT_WALL_EXEC.
- source: FILE_PDE_AUDIT
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: pde_audit_full.md:4 — Executive result
tags:
- atlas/node
- atlas/sources
- layer/source_section_anchor
- status/v05_first_section_anchor
- topic/moving_throat
- type/section_anchor
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Executive result

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `SEC_PDE_PARENT_WALL_EXEC`
> **Status:** `v05_first_section_anchor`
> **Layer:** `source_section_anchor`
> **Type:** `section_anchor`

## Summary

Parent wall status split.

> [!note] Future paper needed
> This node is intentionally anchored to a notes/derivation source until a maintained paper exists.

## Physical Meaning

Parent wall status split.

## Mathematical Role

- Layer: `source_section_anchor`
- Type: `section_anchor`
- Status: `v05_first_section_anchor`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- [[OPEN_PARENT_PROMOTION_S_SIGMA]]

### Status firewalls
- [[FIREWALL_PARENT_WALL_NOT_STRICT]]
- [[FIREWALL_WALL_COEFFS_BRANCH_DATA]]

### Source anchors
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Parent wall status split. |
| `ANCHORS_CLAIM_SECTION` | [[OPEN_PARENT_PROMOTION_S_SIGMA]] | Parent wall status split. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORED_IN` | [[FIREWALL_PARENT_WALL_NOT_STRICT]] | Firewall is anchored in this source/file/section. |
| `ANCHORED_IN` | [[FIREWALL_WALL_COEFFS_BRANCH_DATA]] | Firewall is anchored in this source/file/section. |
| `BACKLINKS_SECTION_ANCHOR` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references source-section anchor SEC_PDE_PARENT_WALL_EXEC. |
| `HAS_SECTION_ANCHOR` | [[FILE_PDE_AUDIT]] | pde_audit_full.md:4 — Executive result |

## Source Anchors

### Source anchor notes
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
