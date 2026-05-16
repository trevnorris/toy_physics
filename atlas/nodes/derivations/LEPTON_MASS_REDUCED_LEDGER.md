---
id: LEPTON_MASS_REDUCED_LEDGER
title: Reduced isolated lepton mass ledger
type: particle_mass_reduction
layer: derivation
status: conditional_reduced_theorem
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Reduced one-parameter mass functional gives 11:2:5 internal partition and D/N support-clock form under declared closure.
future_paper_needed: false
source_files:
- notes/atom_work.md
- notes/lepton_work.md
- notes/lepton_mass_notes.md
- notes/g2/g2_full_output.md
- atom_work.md
- lepton_work.md
- lepton_mass_notes.md
- g2_full_output.md
legacy_sources:
- atom_work.md
- lepton_work.md
- lepton_mass_notes.md
- g2_full_output.md
source_links:
- '[[FILE_LEPTON_MASS]]'
physical_ids:
- PHYS_INTERIOR_SUPPORT
equation_ids:
- EQ_LEPTON_MASS_PARTITION
claim_ids:
- CLAIM_LEPTON_MASS_REDUCED_LEDGER
source_ids:
- FILE_LEPTON_MASS
outgoing_edges:
- target: PHYS_INTERIOR_SUPPORT
  relation: INTERPRETS
  status: reduced
  note: Mass ledger is internal support/rest-energy branch data, not primitive mass.
incoming_edges:
- source: EQ_LEPTON_MASS_PARTITION
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_LEPTON_MASS
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_LEPTON_MASS_REDUCED_LEDGER
  relation: FEEDS_OR_STATUS_OF
  status: reduced_closure_theorem_open_scale
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/conditional_reduced_theorem
- topic/atom
- topic/lepton
- topic/projection
- type/particle_mass_reduction
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Reduced isolated lepton mass ledger

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `LEPTON_MASS_REDUCED_LEDGER`
> **Status:** `conditional_reduced_theorem`
> **Layer:** `derivation`
> **Type:** `particle_mass_reduction`

## Summary

Reduced one-parameter mass functional gives 11:2:5 internal partition and D/N support-clock form under declared closure.

## Physical Meaning

Reduced one-parameter mass functional gives 11:2:5 internal partition and D/N support-clock form under declared closure.

## Mathematical Role

- Layer: `derivation`
- Type: `particle_mass_reduction`
- Status: `conditional_reduced_theorem`

## Atlas Links

### Related physical nodes
- [[PHYS_INTERIOR_SUPPORT]]

### Related math nodes
- none

### Related equations
- [[EQ_LEPTON_MASS_PARTITION]]

### Related claims
- [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_LEPTON_MASS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `INTERPRETS` | [[PHYS_INTERIOR_SUPPORT]] | Mass ledger is internal support/rest-energy branch data, not primitive mass. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_LEPTON_MASS_PARTITION]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_LEPTON_MASS]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_LEPTON_MASS]]

### Source files
- `notes/atom_work.md`
- `notes/lepton_work.md`
- `notes/lepton_mass_notes.md`
- `notes/g2/g2_full_output.md`
- `atom_work.md`
- `lepton_work.md`
- `lepton_mass_notes.md`
- `g2_full_output.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
