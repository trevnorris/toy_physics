---
id: LEPTON_MIXED_BERRY_ROTOR
title: Mixed-sector Berry rotor
type: particle_identity_corridor
layer: derivation
status: conditional_live_corridor
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Same-charge internal corridor using mixed A_w/F_mu_w sector, transverse odd coordinates, and Berry/Magnus first-order term.
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
physical_ids:
- PHYS_MIXED_EM_CORE
claim_ids:
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
status_firewall_ids:
- FIREWALL_LEPTON_CONDITIONAL
source_ids:
- FILE_LEPTON_WORK
outgoing_edges:
- target: PHYS_MIXED_EM_CORE
  relation: REQUIRES
  status: conditional
  note: Same-charge Berry rotor corridor requires mixed core channels to remain live.
incoming_edges:
- source: FILE_LEPTON_WORK
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: FEEDS_OR_STATUS_OF
  status: exact_bookkeeping_firewall
  note: Claim feeds this downstream object, output, or open gate.
- source: FIREWALL_LEPTON_CONDITIONAL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/conditional_live_corridor
- topic/atom
- topic/charge
- topic/lepton
- topic/maxwell
- type/particle_identity_corridor
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Mixed-sector Berry rotor

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `LEPTON_MIXED_BERRY_ROTOR`
> **Status:** `conditional_live_corridor`
> **Layer:** `derivation`
> **Type:** `particle_identity_corridor`

## Summary

Same-charge internal corridor using mixed A_w/F_mu_w sector, transverse odd coordinates, and Berry/Magnus first-order term.

## Physical Meaning

Same-charge internal corridor using mixed A_w/F_mu_w sector, transverse odd coordinates, and Berry/Magnus first-order term.

## Mathematical Role

- Layer: `derivation`
- Type: `particle_identity_corridor`
- Status: `conditional_live_corridor`

## Atlas Links

### Related physical nodes
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_LEPTON_CONDITIONAL]]

### Source anchors
- [[FILE_LEPTON_WORK]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REQUIRES` | [[PHYS_MIXED_EM_CORE]] | Same-charge Berry rotor corridor requires mixed core channels to remain live. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_LEPTON_WORK]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Claim feeds this downstream object, output, or open gate. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_LEPTON_CONDITIONAL]] | Firewall preserves this correct status boundary. |

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
