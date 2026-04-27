---
id: MATH_MAXWELL_MIXED_KERNEL
title: Maxwell/mixed kernel
type: reduced_kernel
layer: math_object
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Conservative localized Maxwell/mixed self-energy and outgoing transfer factor for port-active mixed sector.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
math_ids:
- MATH_MIXED_FIELDS_EW_CA
equation_ids:
- EQ_MAXWELL_MIXED_TRANSFER
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
outgoing_edges:
- target: MT_STAGE4_MAXWELL_MIXED
  relation: FORMALIZES
  status: reduced
  note: Matches Stage-4 port-active mixed-sector bridge.
incoming_edges:
- source: EQ_MAXWELL_MIXED_TRANSFER
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_MAXWELL_MIXED_KERNEL.
- source: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_mixed_kernel
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_10_HAMILTONIAN_STABILITY
  relation: GATES
  status: mandatory
  note: Adds internal-block, wall-softening, ghost/Krein, and dark-port failure modes.
- source: MATH_MIXED_FIELDS_EW_CA
  relation: REDUCES_TO
  status: reduced
  note: Mixed fields feed reduced port-active Maxwell/mixed kernel.
- source: MT_V2_09_MAXWELL_MIXED_KERNEL
  relation: VALIDATES
  status: reduced
  note: Audits the Maxwell/mixed kernel and outgoing transfer factor.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_within_reduced_bundle
- topic/maxwell
- topic/moving_throat
- type/reduced_kernel
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Maxwell/mixed kernel

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_MAXWELL_MIXED_KERNEL`  
> **Status:** `exact_within_reduced_bundle`  
> **Layer:** `math_object`  
> **Type:** `reduced_kernel`

## Summary

Conservative localized Maxwell/mixed self-energy and outgoing transfer factor for port-active mixed sector.

## Physical Meaning

Conservative localized Maxwell/mixed self-energy and outgoing transfer factor for port-active mixed sector.

## Mathematical Role

- Layer: `math_object`
- Type: `reduced_kernel`
- Status: `exact_within_reduced_bundle`

## Equation

$$
Z_n
$$

$$
N_n
$$

$$
Delta=Omega_U^2 Omega_W^2-R^2
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MIXED_FIELDS_EW_CA]]

### Related equations
- [[EQ_MAXWELL_MIXED_TRANSFER]]

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FORMALIZES` | [[MT_STAGE4_MAXWELL_MIXED]] | Matches Stage-4 port-active mixed-sector bridge. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_MAXWELL_MIXED_TRANSFER]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references MATH_MAXWELL_MIXED_KERNEL. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]] | Claim feeds this downstream object, output, or open gate. |
| `GATES` | [[MT_V2_10_HAMILTONIAN_STABILITY]] | Adds internal-block, wall-softening, ghost/Krein, and dark-port failure modes. |
| `REDUCES_TO` | [[MATH_MIXED_FIELDS_EW_CA]] | Mixed fields feed reduced port-active Maxwell/mixed kernel. |
| `VALIDATES` | [[MT_V2_09_MAXWELL_MIXED_KERNEL]] | Audits the Maxwell/mixed kernel and outgoing transfer factor. |

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
