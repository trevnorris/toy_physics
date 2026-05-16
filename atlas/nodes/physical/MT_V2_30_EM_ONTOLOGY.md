---
id: MT_V2_30_EM_ONTOLOGY
title: V2-30 EM ontology/status
type: ontology_status
layer: physical_ontology
status: paper_facing
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Records localized Maxwell/mixed sector and charge/circulation firewall.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-30
physical_ids:
- PHYS_CHARGE_BRANCH
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
- PHYS_MIXED_EM_CORE
open_gate_ids:
- OPEN_MIXED_RECIRCULATION
outgoing_edges:
- target: PHYS_CHARGE_BRANCH
  relation: ANCHORS
  status: paper-facing
  note: EM ontology anchors charge dictionary.
- target: PHYS_MAGNETIC_VORTICAL_CIRCULATION
  relation: ANCHORS
  status: paper-facing
  note: EM ontology anchors circulation/magnetism firewall.
- target: OPEN_MIXED_RECIRCULATION
  relation: OPENS_GATE
  status: open
  note: Full recirculation/plumbing law remains future work.
- target: PHYS_MIXED_EM_CORE
  relation: PATCHES_LANGUAGE
  status: paper_facing
  note: EM ontology/status keeps mixed channels live while respecting zero-mode reduction.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/paper_facing
- topic/charge
- topic/maxwell
- topic/moving_throat
- type/ontology_status
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-30 EM ontology/status

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_30_EM_ONTOLOGY`
> **Status:** `paper_facing`
> **Layer:** `physical_ontology`
> **Type:** `ontology_status`

## Summary

Records localized Maxwell/mixed sector and charge/circulation firewall.

## Physical Meaning

Records localized Maxwell/mixed sector and charge/circulation firewall.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `ontology_status`
- Status: `paper_facing`

## Equation

$$
eta_Q,q_*,q_eff distinct from circulation
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_CHARGE_BRANCH]]
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- none

### Open gates
- [[OPEN_MIXED_RECIRCULATION]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[PHYS_CHARGE_BRANCH]] | EM ontology anchors charge dictionary. |
| `ANCHORS` | [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]] | EM ontology anchors circulation/magnetism firewall. |
| `OPENS_GATE` | [[OPEN_MIXED_RECIRCULATION]] | Full recirculation/plumbing law remains future work. |
| `PATCHES_LANGUAGE` | [[PHYS_MIXED_EM_CORE]] | EM ontology/status keeps mixed channels live while respecting zero-mode reduction. |

## Incoming Edges

- none

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
