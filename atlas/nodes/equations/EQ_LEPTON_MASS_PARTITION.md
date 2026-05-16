---
id: EQ_LEPTON_MASS_PARTITION
title: Reduced lepton mass partition
type: equation
layer: equation_anchor
status: conditional_reduced_theorem
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Reduced isolated-defect rest-energy partition under declared closure.
future_paper_needed: false
source_files:
- notes/lepton_mass_notes.md
- notes/lepton_work.md
- lepton_mass_notes.md
- lepton_work.md
legacy_sources:
- lepton_mass_notes.md
- lepton_work.md
source_links:
- '[[FILE_LEPTON_MASS]]'
- '[[FILE_LEPTON_WORK]]'
claim_ids:
- CLAIM_LEPTON_MASS_REDUCED_LEDGER
source_ids:
- FILE_LEPTON_MASS
- FILE_LEPTON_WORK
outgoing_edges:
- target: LEPTON_MASS_REDUCED_LEDGER
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_LEPTON_MASS_REDUCED_LEDGER
  relation: SUPPORTS_CLAIM
  status: reduced_closure_theorem_open_scale
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_LEPTON_MASS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_LEPTON_MASS_PARTITION.
- source: FILE_LEPTON_MASS
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_LEPTON_MASS_PARTITION.
- source: FILE_LEPTON_WORK
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_LEPTON_MASS_PARTITION.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/conditional_reduced_theorem
- topic/lepton
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Reduced lepton mass partition

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_LEPTON_MASS_PARTITION`
> **Status:** `conditional_reduced_theorem`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Reduced isolated-defect rest-energy partition under declared closure.

## Physical Meaning

Reduced isolated-defect rest-energy partition under declared closure.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `conditional_reduced_theorem`
- Parent node: `LEPTON_MASS_REDUCED_LEDGER`

## Equation

$$
E_w:E_f:E_PV = 11:2:5,    F_*=(18/11)E_w
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_LEPTON_MASS]]
- [[FILE_LEPTON_WORK]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[LEPTON_MASS_REDUCED_LEDGER]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_MASS]] | Paper backlink block references EQ_LEPTON_MASS_PARTITION. |
| `CONTAINS_EQUATION` | [[FILE_LEPTON_MASS]] | Source artifact contains or supports EQ_LEPTON_MASS_PARTITION. |
| `CONTAINS_EQUATION` | [[FILE_LEPTON_WORK]] | Source artifact contains or supports EQ_LEPTON_MASS_PARTITION. |

## Source Anchors

### Source anchor notes
- [[FILE_LEPTON_MASS]]
- [[FILE_LEPTON_WORK]]

### Source files
- `notes/lepton_mass_notes.md`
- `notes/lepton_work.md`
- `lepton_mass_notes.md`
- `lepton_work.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
