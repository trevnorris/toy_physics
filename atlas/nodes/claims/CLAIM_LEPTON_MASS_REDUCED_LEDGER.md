---
id: CLAIM_LEPTON_MASS_REDUCED_LEDGER
title: Lepton mass partition as reduced defect-energy theorem
type: claim
layer: claim_theorem
status: reduced_closure_theorem_open_scale
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The isolated lepton mass ledger fixes the 11:2:5 partition and total/rest-energy relation within the reduced closure, while the absolute branch scale remains open.
future_paper_needed: false
source_links:
- '[[FILE_LEPTON_MASS]]'
- '[[FILE_LEPTON_WORK]]'
- '[[SEC_LEPTON_MASS_THEOREM]]'
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
- PHYS_INTERIOR_SUPPORT
equation_ids:
- EQ_LEPTON_MASS_PARTITION
claim_ids:
- CLAIM_STAGE1_GEOMETRY_LIFT
source_ids:
- FILE_LEPTON_MASS
- FILE_LEPTON_WORK
- SEC_LEPTON_MASS_THEOREM
outgoing_edges:
- target: LEPTON_MASS_REDUCED_LEDGER
  relation: FEEDS_OR_STATUS_OF
  status: reduced_closure_theorem_open_scale
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: PHYSICAL_THROAT_CONTEXT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_LEPTON_MASS_THEOREM
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mass partition theorem.
- source: BACKLINK_LEPTON_MASS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_LEPTON_MASS_REDUCED_LEDGER.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: reduced_closure_theorem_open_scale
  note: 'Claim class: reduced_sector_consequence'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: GROUNDS_PHYSICAL_MEANING
  status: reduced_closure_theorem_open_scale
  note: Physical ontology object grounded by this claim.
- source: PHYS_INTERIOR_SUPPORT
  relation: GROUNDS_PHYSICAL_MEANING
  status: reduced_closure_theorem_open_scale
  note: Physical ontology object grounded by this claim.
- source: FILE_LEPTON_MASS
  relation: OWNS_OR_ANCHORS_CLAIM
  status: reduced_closure_theorem_open_scale
  note: Source artifact anchors this claim.
- source: FILE_LEPTON_WORK
  relation: OWNS_OR_ANCHORS_CLAIM
  status: reduced_closure_theorem_open_scale
  note: Source artifact anchors this claim.
- source: EQ_LEPTON_MASS_PARTITION
  relation: SUPPORTS_CLAIM
  status: reduced_closure_theorem_open_scale
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/reduced_closure_theorem_open_scale
- topic/lepton
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Lepton mass partition as reduced defect-energy theorem

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_LEPTON_MASS_REDUCED_LEDGER`
> **Status:** `reduced_closure_theorem_open_scale`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The isolated lepton mass ledger fixes the 11:2:5 partition and total/rest-energy relation within the reduced closure, while the absolute branch scale remains open.

## Claim

The isolated lepton mass ledger fixes the 11:2:5 partition and total/rest-energy relation within the reduced closure, while the absolute branch scale remains open.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

The isolated lepton mass ledger fixes the 11:2:5 partition and total/rest-energy relation within the reduced closure, while the absolute branch scale remains open.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `reduced_closure_theorem_open_scale`
- Outputs: `LEPTON_MASS_REDUCED_LEDGER`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]
- [[PHYS_INTERIOR_SUPPORT]]

### Related math nodes
- none

### Related equations
- [[EQ_LEPTON_MASS_PARTITION]]

### Related claims
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_LEPTON_MASS]]
- [[FILE_LEPTON_WORK]]
- [[SEC_LEPTON_MASS_THEOREM]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[LEPTON_MASS_REDUCED_LEDGER]] | Claim feeds this downstream object, output, or open gate. |
| `PHYSICAL_THROAT_CONTEXT_FOR` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_LEPTON_MASS_THEOREM]] | Mass partition theorem. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_MASS]] | Paper backlink block references CLAIM_LEPTON_MASS_REDUCED_LEDGER. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: reduced_sector_consequence |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_INTERIOR_SUPPORT]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_LEPTON_MASS]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_LEPTON_WORK]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_LEPTON_MASS_PARTITION]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_LEPTON_MASS]]
- [[FILE_LEPTON_WORK]]
- [[SEC_LEPTON_MASS_THEOREM]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
