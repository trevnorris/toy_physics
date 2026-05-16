---
id: LEPTON_HALF_INTEGER_CONDITIONAL
title: Conditional half-integer discretizer
type: particle_identity_gate
layer: derivation
status: conditional_open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Half-integer quantizer follows only if selective tau subbundle, central-sign holonomy, and autonomous standing-wave soliton closure hold.
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
- '[[FILE_LEPTON_WORK]]'
equation_ids:
- EQ_LEPTON_HALF_INTEGER
claim_ids:
- CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
open_gate_ids:
- OPEN_LEPTON_SPIN_DISCRETIZER
source_ids:
- FILE_LEPTON_WORK
outgoing_edges:
- target: OPEN_LEPTON_SPIN_DISCRETIZER
  relation: REFINES
  status: open
  note: Half-integer spin route remains conditional on autonomous subbundle closure.
incoming_edges:
- source: EQ_LEPTON_HALF_INTEGER
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_LEPTON_WORK
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
  relation: FEEDS_OR_STATUS_OF
  status: conditional_open
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/conditional_open
- topic/atom
- topic/lepton
- type/particle_identity_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Conditional half-integer discretizer

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `LEPTON_HALF_INTEGER_CONDITIONAL`
> **Status:** `conditional_open`
> **Layer:** `derivation`
> **Type:** `particle_identity_gate`

## Summary

Half-integer quantizer follows only if selective tau subbundle, central-sign holonomy, and autonomous standing-wave soliton closure hold.

## Physical Meaning

Half-integer quantizer follows only if selective tau subbundle, central-sign holonomy, and autonomous standing-wave soliton closure hold.

## Mathematical Role

- Layer: `derivation`
- Type: `particle_identity_gate`
- Status: `conditional_open`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_LEPTON_HALF_INTEGER]]

### Related claims
- [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]]

### Open gates
- [[OPEN_LEPTON_SPIN_DISCRETIZER]]

### Status firewalls
- none

### Source anchors
- [[FILE_LEPTON_WORK]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REFINES` | [[OPEN_LEPTON_SPIN_DISCRETIZER]] | Half-integer spin route remains conditional on autonomous subbundle closure. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_LEPTON_HALF_INTEGER]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_LEPTON_WORK]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_LEPTON_WORK]]

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
