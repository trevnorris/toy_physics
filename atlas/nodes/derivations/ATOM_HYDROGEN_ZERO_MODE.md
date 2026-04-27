---
id: ATOM_HYDROGEN_ZERO_MODE
title: Hydrogenic zero-mode sector
type: atomic_reduction
layer: derivation
status: reduced_sector_consequence
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Uses zero-mode Maxwell and lowest transverse matter mode to recover hydrogenic Coulomb/Bohr-scale sector from the 4D action.
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
- '[[FILE_ATOM_WORK]]'
math_ids:
- MATH_QSTAR_QEFF
- MATH_ZERO_MODE_MAXWELL
equation_ids:
- EQ_QEFF_NORMALIZATION
- EQ_ZERO_MODE_MAXWELL
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
status_firewall_ids:
- FIREWALL_ATOM_REDUCED_SECTOR
source_ids:
- FILE_ATOM_WORK
outgoing_edges:
- target: MATH_QSTAR_QEFF
  relation: IMPORTS
  status: reduced
  note: Hydrogenic charge scale uses corrected thickness-controlled q_eff.
- target: MATH_ZERO_MODE_MAXWELL
  relation: IMPORTS
  status: reduced
  note: Hydrogenic reduction uses controlled zero-mode Maxwell/Coulomb sector.
incoming_edges:
- source: FILE_ATOM_WORK
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: EQ_QEFF_NORMALIZATION
  relation: FEEDS
  status: reduced
  note: Sets observed charge scale in the reduced atomic sector.
- source: EQ_ZERO_MODE_MAXWELL
  relation: FEEDS
  status: reduced
  note: Provides the Coulomb sector used in first-pass atomic reduction.
- source: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: FEEDS_OR_STATUS_OF
  status: reduced_sector_consequence
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: FEEDS_OR_STATUS_OF
  status: exact_bookkeeping_firewall
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: FEEDS_OR_STATUS_OF
  status: controlled_reduction
  note: Claim feeds this downstream object, output, or open gate.
- source: FIREWALL_ATOM_REDUCED_SECTOR
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/reduced_sector_consequence
- topic/atom
- topic/lepton
- topic/maxwell
- topic/projection
- type/atomic_reduction
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Hydrogenic zero-mode sector

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `ATOM_HYDROGEN_ZERO_MODE`  
> **Status:** `reduced_sector_consequence`  
> **Layer:** `derivation`  
> **Type:** `atomic_reduction`

## Summary

Uses zero-mode Maxwell and lowest transverse matter mode to recover hydrogenic Coulomb/Bohr-scale sector from the 4D action.

## Physical Meaning

Uses zero-mode Maxwell and lowest transverse matter mode to recover hydrogenic Coulomb/Bohr-scale sector from the 4D action.

## Mathematical Role

- Layer: `derivation`
- Type: `atomic_reduction`
- Status: `reduced_sector_consequence`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_QSTAR_QEFF]]
- [[MATH_ZERO_MODE_MAXWELL]]

### Related equations
- [[EQ_QEFF_NORMALIZATION]]
- [[EQ_ZERO_MODE_MAXWELL]]

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_ATOM_REDUCED_SECTOR]]

### Source anchors
- [[FILE_ATOM_WORK]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `IMPORTS` | [[MATH_QSTAR_QEFF]] | Hydrogenic charge scale uses corrected thickness-controlled q_eff. |
| `IMPORTS` | [[MATH_ZERO_MODE_MAXWELL]] | Hydrogenic reduction uses controlled zero-mode Maxwell/Coulomb sector. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_ATOM_WORK]] | File anchor documents this node. |
| `FEEDS` | [[EQ_QEFF_NORMALIZATION]] | Sets observed charge scale in the reduced atomic sector. |
| `FEEDS` | [[EQ_ZERO_MODE_MAXWELL]] | Provides the Coulomb sector used in first-pass atomic reduction. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Claim feeds this downstream object, output, or open gate. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ATOM_REDUCED_SECTOR]] | Firewall preserves this correct status boundary. |

## Source Anchors

### Source anchor notes
- [[FILE_ATOM_WORK]]

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
