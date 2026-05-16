---
id: PHYS_OUTGOING_QUADRUPOLE_PORT
title: Passive/outgoing quadrupole port
type: radiative_port
layer: physical_ontology
status: open_actual_branch_data
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Surviving universal dissipative/tail bridge after scalar and dipole demotion; requires actual branch normalization.
future_paper_needed: false
source_files:
- research/4d_2_5pn/paper/4d_2_5pn.tex
- research/4d_4pn/paper/4d_4pn.tex
- notes/pde_audit_full.md
- 4d_2_5pn_summary.md
- 4d_4pn_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_2_5pn_summary.md
- 4d_4pn_summary.md
- pde_audit_full.md
physical_ids:
- PHYS_REG_OUTGOING_PORT
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
- CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
- CLAIM_STAGE022_P0_TARGET
- CLAIM_STAGE025_031_SELECTED_BRANCH
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
outgoing_edges:
- target: CLAIM_25PN_QUAD_NARROWING
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_theorem_open_normalization
  note: Physical ontology object grounded by this claim.
- target: CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
  relation: GROUNDS_PHYSICAL_MEANING
  status: local_closed_tail_conditional
  note: Physical ontology object grounded by this claim.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_gauge_invariant_with_reduced_uses
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_mixed_kernel
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE022_P0_TARGET
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_grouped_bridge
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE025_031_SELECTED_BRANCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_selected_reduced_branch
  note: Physical ontology object grounded by this claim.
- target: OPEN_QUAD_NORMALIZATION
  relation: REQUIRES
  status: open
  note: Outgoing port normalization must match universal quadrupole coefficient.
incoming_edges:
- source: BACKLINK_25PN
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_OUTGOING_QUADRUPOLE_PORT.
- source: BACKLINK_4PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_OUTGOING_QUADRUPOLE_PORT.
- source: PHYS_REG_OUTGOING_PORT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/open_actual_branch_data
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/radiative_port
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Passive/outgoing quadrupole port

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_OUTGOING_QUADRUPOLE_PORT`
> **Status:** `open_actual_branch_data`
> **Layer:** `physical_ontology`
> **Type:** `radiative_port`

## Summary

Surviving universal dissipative/tail bridge after scalar and dipole demotion; requires actual branch normalization.

## Physical Meaning

Surviving universal dissipative/tail bridge after scalar and dipole demotion; requires actual branch normalization.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `radiative_port`
- Status: `open_actual_branch_data`

## Equation

$$
i omega^5
$$

$$
Gamma_2^port=a^5/(27 c_s^5)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_REG_OUTGOING_PORT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_25PN_QUAD_NARROWING]]
- [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]]
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]
- [[CLAIM_STAGE022_P0_TARGET]]
- [[CLAIM_STAGE025_031_SELECTED_BRANCH]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_25PN_QUAD_NARROWING]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE022_P0_TARGET]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE025_031_SELECTED_BRANCH]] | Physical ontology object grounded by this claim. |
| `REQUIRES` | [[OPEN_QUAD_NORMALIZATION]] | Outgoing port normalization must match universal quadrupole coefficient. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_25PN]] | Paper backlink block references PHYS_OUTGOING_QUADRUPOLE_PORT. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4PN_FULL]] | Paper backlink block references PHYS_OUTGOING_QUADRUPOLE_PORT. |
| `LINKS_TO` | [[PHYS_REG_OUTGOING_PORT]] | Physical register entry links to graph object. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d_2_5pn/paper/4d_2_5pn.tex`
- `research/4d_4pn/paper/4d_4pn.tex`
- `notes/pde_audit_full.md`
- `4d_2_5pn_summary.md`
- `4d_4pn_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
