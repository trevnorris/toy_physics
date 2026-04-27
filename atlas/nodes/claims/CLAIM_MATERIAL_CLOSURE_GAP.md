---
id: CLAIM_MATERIAL_CLOSURE_GAP
title: Material closure must be solved on same branch
type: claim
layer: claim_theorem
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Density, sound speed, effective light cone, flux feedback, and material screening must be exported on the same branch as the response packet.
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_MATERIAL_GAP]]'
physical_ids:
- PHYS_MATERIAL_CLOSURE
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
open_gate_ids:
- MT_V2_29_MATERIAL_GAP
- OPEN_MATERIAL_CLOSURE
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_PDE_MATERIAL_GAP
outgoing_edges:
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: BLOCKS_FULL_BRANCH_REALIZATION_IF_UNSOLVED
  status: active
  note: Claim-level dependency added in v0.4.
- target: MT_V2_29_MATERIAL_GAP
  relation: FEEDS_OR_STATUS_OF
  status: open
  note: Claim feeds this downstream object, output, or open gate.
- target: OPEN_MATERIAL_CLOSURE
  relation: FEEDS_OR_STATUS_OF
  status: open
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_PDE_MATERIAL_GAP
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Material closure gap.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_MATERIAL_CLOSURE_GAP.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: open
  note: 'Claim class: open_gate'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_MATERIAL_CLOSURE
  relation: GROUNDS_PHYSICAL_MEANING
  status: open
  note: Physical ontology object grounded by this claim.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: open
  note: Source artifact anchors this claim.
- source: QUERY_MATERIAL_CLOSURE
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/open
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Material closure must be solved on same branch

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_MATERIAL_CLOSURE_GAP`  
> **Status:** `open`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

Density, sound speed, effective light cone, flux feedback, and material screening must be exported on the same branch as the response packet.

## Claim

Density, sound speed, effective light cone, flux feedback, and material screening must be exported on the same branch as the response packet.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

Density, sound speed, effective light cone, flux feedback, and material screening must be exported on the same branch as the response packet.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `open`
- Outputs: `OPEN_MATERIAL_CLOSURE`, `MT_V2_29_MATERIAL_GAP`

## Atlas Links

### Related physical nodes
- [[PHYS_MATERIAL_CLOSURE]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]

### Open gates
- [[MT_V2_29_MATERIAL_GAP]]
- [[OPEN_MATERIAL_CLOSURE]]

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_MATERIAL_GAP]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `BLOCKS_FULL_BRANCH_REALIZATION_IF_UNSOLVED` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Claim-level dependency added in v0.4. |
| `FEEDS_OR_STATUS_OF` | [[MT_V2_29_MATERIAL_GAP]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_MATERIAL_CLOSURE]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_MATERIAL_GAP]] | Material closure gap. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_MATERIAL_CLOSURE_GAP. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: open_gate |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MATERIAL_CLOSURE]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `STARTS_AT` | [[QUERY_MATERIAL_CLOSURE]] | Query validation start node. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_MATERIAL_GAP]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
