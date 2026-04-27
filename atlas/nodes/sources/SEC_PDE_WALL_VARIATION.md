---
id: SEC_PDE_WALL_VARIATION
title: 3. Variation after adding the quadratic wall action
type: section_anchor
layer: source_section_anchor
status: v05_first_section_anchor
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: S_eta gives wall PDE.
source_kind: future_paper_note
future_paper_needed: true
source_files:
- notes/pde_audit_full.md
source_links:
- '[[FILE_PDE_AUDIT]]'
claim_ids:
- CLAIM_PARENT_WALL_STATUS_SPLIT
- CLAIM_STAGE2_AL_RECOVERY
status_firewall_ids:
- FIREWALL_PARENT_WALL_NOT_STRICT
source_ids:
- FILE_PDE_AUDIT
outgoing_edges:
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: S_eta gives wall PDE.
- target: CLAIM_STAGE2_AL_RECOVERY
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: S_eta gives wall PDE.
incoming_edges:
- source: FIREWALL_PARENT_WALL_NOT_STRICT
  relation: ANCHORED_IN
  status: v07
  note: Firewall is anchored in this source/file/section.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_SECTION_ANCHOR
  status: v06
  note: Paper backlink block references source-section anchor SEC_PDE_WALL_VARIATION.
- source: FILE_PDE_AUDIT
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: pde_audit_full.md:188 — 3. Variation after adding the quadratic wall action
tags:
- atlas/node
- atlas/sources
- layer/source_section_anchor
- status/v05_first_section_anchor
- topic/moving_throat
- type/section_anchor
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 3. Variation after adding the quadratic wall action

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `SEC_PDE_WALL_VARIATION`  
> **Status:** `v05_first_section_anchor`  
> **Layer:** `source_section_anchor`  
> **Type:** `section_anchor`

## Summary

S_eta gives wall PDE.

> [!note] Future paper needed
> This node is intentionally anchored to a notes/derivation source until a maintained paper exists.

## Physical Meaning

S_eta gives wall PDE.

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
- [[CLAIM_STAGE2_AL_RECOVERY]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_PARENT_WALL_NOT_STRICT]]

### Source anchors
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | S_eta gives wall PDE. |
| `ANCHORS_CLAIM_SECTION` | [[CLAIM_STAGE2_AL_RECOVERY]] | S_eta gives wall PDE. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORED_IN` | [[FIREWALL_PARENT_WALL_NOT_STRICT]] | Firewall is anchored in this source/file/section. |
| `BACKLINKS_SECTION_ANCHOR` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references source-section anchor SEC_PDE_WALL_VARIATION. |
| `HAS_SECTION_ANCHOR` | [[FILE_PDE_AUDIT]] | pde_audit_full.md:188 — 3. Variation after adding the quadratic wall action |

## Source Anchors

### Source anchor notes
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
