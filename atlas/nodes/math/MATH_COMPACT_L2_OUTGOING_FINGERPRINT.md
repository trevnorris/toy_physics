---
id: MATH_COMPACT_L2_OUTGOING_FINGERPRINT
title: Compact outgoing l=2 fingerprint
type: outgoing_fingerprint
layer: math_object
status: exact_reduced
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Outgoing l=2 normalized response has even coefficients plus the leading i omega^5 fingerprint.
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
- PHYS_REG_OUTGOING_PORT
equation_ids:
- EQ_COMPACT_L2_FINGERPRINT
claim_ids:
- CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
outgoing_edges:
- target: OPEN_QUAD_NORMALIZATION
  relation: FEEDS
  status: conditional
  note: The i omega^5 coefficient feeds the quadrupole normalization target.
incoming_edges:
- source: EQ_COMPACT_L2_FINGERPRINT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: MT_V2_14_L2_OUTGOING
  relation: DERIVES
  status: exact_reduced
  note: Derives compact outgoing l=2 fingerprint.
- source: CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_mixed_kernel
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_REG_OUTGOING_PORT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_reduced
- topic/moving_throat
- type/outgoing_fingerprint
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Compact outgoing l=2 fingerprint

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_COMPACT_L2_OUTGOING_FINGERPRINT`  
> **Status:** `exact_reduced`  
> **Layer:** `math_object`  
> **Type:** `outgoing_fingerprint`

## Summary

Outgoing l=2 normalized response has even coefficients plus the leading i omega^5 fingerprint.

## Physical Meaning

Outgoing l=2 normalized response has even coefficients plus the leading i omega^5 fingerprint.

## Mathematical Role

- Layer: `math_object`
- Type: `outgoing_fingerprint`
- Status: `exact_reduced`

## Equation

$$
1+a^2 omega^2/(9c_s^2)+4a^4 omega^4/(81c_s^4)+i a^5 omega^5/(27 c_s^5)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_REG_OUTGOING_PORT]]

### Related math nodes
- none

### Related equations
- [[EQ_COMPACT_L2_FINGERPRINT]]

### Related claims
- [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[OPEN_QUAD_NORMALIZATION]] | The i omega^5 coefficient feeds the quadrupole normalization target. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_COMPACT_L2_FINGERPRINT]] | Equation anchor belongs to or formalizes this graph node. |
| `DERIVES` | [[MT_V2_14_L2_OUTGOING]] | Derives compact outgoing l=2 fingerprint. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE4_MIXED_OUTGOING_TRANSFER]] | Claim feeds this downstream object, output, or open gate. |
| `LINKS_TO` | [[PHYS_REG_OUTGOING_PORT]] | Physical register entry links to graph object. |

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
