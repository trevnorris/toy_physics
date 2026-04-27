---
id: PHYS_MAGNETIC_VORTICAL_CIRCULATION
title: Magnetic/vortical circulation sector
type: topological_sector
layer: physical_ontology
status: exact_identity_and_open_plumbing
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Fluxoid/circulation belongs to tangential magnetic/vortical holonomy, not the electric-charge dictionary.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/pde_audit_full.md
- 4d_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- pde_audit_full.md:V2-30
physical_ids:
- MT_V2_30_EM_ONTOLOGY
math_ids:
- MATH_FLUXOID
- MATH_QSTAR_QEFF
claim_ids:
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
- CLAIM_MIXED_RECIRCULATION_OPEN
open_gate_ids:
- OPEN_LEPTON_SPIN_DISCRETIZER
outgoing_edges:
- target: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_bookkeeping_firewall
  note: Physical ontology object grounded by this claim.
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: GROUNDS_PHYSICAL_MEANING
  status: open
  note: Physical ontology object grounded by this claim.
- target: MATH_FLUXOID
  relation: REPRESENTED_BY
  status: exact/open plumbing
  note: Fluxoid/circulation represented by quantized holonomy law.
incoming_edges:
- source: MT_V2_30_EM_ONTOLOGY
  relation: ANCHORS
  status: paper-facing
  note: EM ontology anchors circulation/magnetism firewall.
- source: MATH_QSTAR_QEFF
  relation: DISTINGUISHES_FROM
  status: firewall
  note: Charge dictionary must not be conflated with circulation.
- source: OPEN_LEPTON_SPIN_DISCRETIZER
  relation: DISTINGUISHES_FROM
  status: open
  note: Spin-like state must not reduce to ordinary orbital/vortical circulation.
- source: NEG_QUERY_CHARGE_FROM_CIRCULATION
  relation: STARTS_AT
  status: v07
  note: Negative query starts from PHYS_MAGNETIC_VORTICAL_CIRCULATION.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_identity_and_open_plumbing
- topic/charge
- topic/moving_throat
- type/topological_sector
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Magnetic/vortical circulation sector

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_MAGNETIC_VORTICAL_CIRCULATION`  
> **Status:** `exact_identity_and_open_plumbing`  
> **Layer:** `physical_ontology`  
> **Type:** `topological_sector`

## Summary

Fluxoid/circulation belongs to tangential magnetic/vortical holonomy, not the electric-charge dictionary.

## Physical Meaning

Fluxoid/circulation belongs to tangential magnetic/vortical holonomy, not the electric-charge dictionary.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `topological_sector`
- Status: `exact_identity_and_open_plumbing`

## Equation

$$
∮(∂_i θ - q_* A_i/hbar) dl^i = 2πn
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[MT_V2_30_EM_ONTOLOGY]]

### Related math nodes
- [[MATH_FLUXOID]]
- [[MATH_QSTAR_QEFF]]

### Related equations
- none

### Related claims
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]

### Open gates
- [[OPEN_LEPTON_SPIN_DISCRETIZER]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Physical ontology object grounded by this claim. |
| `REPRESENTED_BY` | [[MATH_FLUXOID]] | Fluxoid/circulation represented by quantized holonomy law. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MT_V2_30_EM_ONTOLOGY]] | EM ontology anchors circulation/magnetism firewall. |
| `DISTINGUISHES_FROM` | [[MATH_QSTAR_QEFF]] | Charge dictionary must not be conflated with circulation. |
| `DISTINGUISHES_FROM` | [[OPEN_LEPTON_SPIN_DISCRETIZER]] | Spin-like state must not reduce to ordinary orbital/vortical circulation. |
| `STARTS_AT` | [[NEG_QUERY_CHARGE_FROM_CIRCULATION]] | Negative query starts from PHYS_MAGNETIC_VORTICAL_CIRCULATION. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
