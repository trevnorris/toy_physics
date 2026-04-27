---
id: PHYS_OPEN_CONDUIT
title: Open finite-radius conduit
type: boundary_condition
layer: physical_ontology
status: effective_branch_condition
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Physical branch is an open finite-radius throat/outlet; hard cap is obsolete toy idealization unless declared as negative control.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-28
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
claim_ids:
- CLAIM_MIXED_RECIRCULATION_OPEN
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
outgoing_edges:
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: GROUNDS_PHYSICAL_MEANING
  status: open
  note: Physical ontology object grounded by this claim.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_gauge_invariant_with_reduced_uses
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_projection_plus_controlled_hook
  note: Physical ontology object grounded by this claim.
incoming_edges:
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_OPEN_CONDUIT.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: HAS_BOUNDARY_CLASS
  status: effective branch
  note: Physical branch is open finite-radius conduit.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/effective_branch_condition
- topic/moving_throat
- type/boundary_condition
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Open finite-radius conduit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_OPEN_CONDUIT`  
> **Status:** `effective_branch_condition`  
> **Layer:** `physical_ontology`  
> **Type:** `boundary_condition`

## Summary

Physical branch is an open finite-radius throat/outlet; hard cap is obsolete toy idealization unless declared as negative control.

## Physical Meaning

Physical branch is an open finite-radius throat/outlet; hard cap is obsolete toy idealization unless declared as negative control.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `boundary_condition`
- Status: `effective_branch_condition`

## Equation

$$
R(L)>0
$$

$$
open/radiative branch
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Physical ontology object grounded by this claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references PHYS_OPEN_CONDUIT. |
| `HAS_BOUNDARY_CLASS` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Physical branch is open finite-radius conduit. |

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
