---
id: MATH_STF_SOURCE_MAP
title: STF source-map theorem
type: source_map
layer: math_object
status: exact_angular
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Canonical real STF angular basis maps orbital/worldtube quadrupole sources to grouped P2 ports with mhat_ang=1.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
physical_ids:
- PHYS_REG_P2_QUAD
claim_ids:
- CLAIM_STAGE024_O3_ISOTROPY
open_gate_ids:
- OPEN_SOURCE_PORT_NORMALIZATION
status_firewall_ids:
- FIREWALL_ANGULAR_NOT_RADIAL
outgoing_edges:
- target: OPEN_SOURCE_PORT_NORMALIZATION
  relation: NARROWS
  status: open
  note: Leaves radial/axial source/port normalization only.
incoming_edges:
- source: BACKLINK_25PN
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_STF_SOURCE_MAP.
- source: CLAIM_STAGE024_O3_ISOTROPY
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_REG_P2_QUAD
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FIREWALL_ANGULAR_NOT_RADIAL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: MT_V2_12_STF_SOURCE_MAP
  relation: PROVES
  status: exact_angular
  note: Closes angular source-map identity.
- source: NEG_QUERY_STF_ANGULAR_MAP_CLOSES_RADIAL_NORMALIZATION
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MATH_STF_SOURCE_MAP.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_angular
- topic/moving_throat
- topic/quadrupole
- type/source_map
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# STF source-map theorem

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_STF_SOURCE_MAP`
> **Status:** `exact_angular`
> **Layer:** `math_object`
> **Type:** `source_map`

## Summary

Canonical real STF angular basis maps orbital/worldtube quadrupole sources to grouped P2 ports with mhat_ang=1.

## Physical Meaning

Canonical real STF angular basis maps orbital/worldtube quadrupole sources to grouped P2 ports with mhat_ang=1.

## Mathematical Role

- Layer: `math_object`
- Type: `source_map`
- Status: `exact_angular`

## Equation

$$
mhat_ang=1
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_REG_P2_QUAD]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_STAGE024_O3_ISOTROPY]]

### Open gates
- [[OPEN_SOURCE_PORT_NORMALIZATION]]

### Status firewalls
- [[FIREWALL_ANGULAR_NOT_RADIAL]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `NARROWS` | [[OPEN_SOURCE_PORT_NORMALIZATION]] | Leaves radial/axial source/port normalization only. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_25PN]] | Paper backlink block references MATH_STF_SOURCE_MAP. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE024_O3_ISOTROPY]] | Claim feeds this downstream object, output, or open gate. |
| `LINKS_TO` | [[PHYS_REG_P2_QUAD]] | Physical register entry links to graph object. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ANGULAR_NOT_RADIAL]] | Firewall preserves this correct status boundary. |
| `PROVES` | [[MT_V2_12_STF_SOURCE_MAP]] | Closes angular source-map identity. |
| `STARTS_AT` | [[NEG_QUERY_STF_ANGULAR_MAP_CLOSES_RADIAL_NORMALIZATION]] | Negative query starts from MATH_STF_SOURCE_MAP. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_pde_program_compact.md`
- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
