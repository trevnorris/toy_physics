---
id: CLAIM_MIXED_RECIRCULATION_OPEN
title: Mixed recirculation/plumbing law remains open
type: claim
layer: claim_theorem
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The EM ontology needs a closed mixed recirculation law relating circulation, throat intake, mixed transport, brane magnetic fields, and source current.
future_paper_needed: false
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_EM_STATUS]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
- PHYS_MIXED_EM_CORE
- PHYS_OPEN_CONDUIT
equation_ids:
- EQ_MAXWELL_MIXED_TRANSFER
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
open_gate_ids:
- OPEN_MIXED_RECIRCULATION
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_PDE_EM_STATUS
outgoing_edges:
- target: OPEN_MIXED_RECIRCULATION
  relation: FEEDS_OR_STATUS_OF
  status: open
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: REFINES_OPEN_EM_REQUIREMENT
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_PDE_EM_STATUS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: EM ontology and recirculation/plumbing status.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: open
  note: 'Claim class: open_gate'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_MAGNETIC_VORTICAL_CIRCULATION
  relation: GROUNDS_PHYSICAL_MEANING
  status: open
  note: Physical ontology object grounded by this claim.
- source: PHYS_MIXED_EM_CORE
  relation: GROUNDS_PHYSICAL_MEANING
  status: open
  note: Physical ontology object grounded by this claim.
- source: PHYS_OPEN_CONDUIT
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
- source: EQ_MAXWELL_MIXED_TRANSFER
  relation: SUPPORTS_CLAIM
  status: open
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/open
- topic/maxwell
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Mixed recirculation/plumbing law remains open

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_MIXED_RECIRCULATION_OPEN`  
> **Status:** `open`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

The EM ontology needs a closed mixed recirculation law relating circulation, throat intake, mixed transport, brane magnetic fields, and source current.

## Claim

The EM ontology needs a closed mixed recirculation law relating circulation, throat intake, mixed transport, brane magnetic fields, and source current.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

The EM ontology needs a closed mixed recirculation law relating circulation, throat intake, mixed transport, brane magnetic fields, and source current.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `open`
- Outputs: `OPEN_MIXED_RECIRCULATION`

## Atlas Links

### Related physical nodes
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]
- [[PHYS_MIXED_EM_CORE]]
- [[PHYS_OPEN_CONDUIT]]

### Related math nodes
- none

### Related equations
- [[EQ_MAXWELL_MIXED_TRANSFER]]

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]

### Open gates
- [[OPEN_MIXED_RECIRCULATION]]

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_EM_STATUS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[OPEN_MIXED_RECIRCULATION]] | Claim feeds this downstream object, output, or open gate. |
| `REFINES_OPEN_EM_REQUIREMENT` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_EM_STATUS]] | EM ontology and recirculation/plumbing status. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: open_gate |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MIXED_EM_CORE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_OPEN_CONDUIT]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_MAXWELL_MIXED_TRANSFER]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_EM_STATUS]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
