---
id: MATH_WEAK_AXISYM_SPLITTING
title: Weak-axisymmetric splitting law
type: anisotropy_fingerprint
layer: math_object
status: exact_angular_first_order
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Pure weak-axisymmetric l=2 perturbation forces grouped signature (20,21,22)~(1,1/2,-1), equivalently b=3a.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_pde_program_compact.md
- pde_audit_full.md
- moving_throat_pde_program_compact.md
legacy_sources:
- pde_audit_full.md
- moving_throat_pde_program_compact.md
equation_ids:
- EQ_WEAK_AXISYM_SIGNATURE
claim_ids:
- CLAIM_STAGE024_O3_ISOTROPY
open_gate_ids:
- OPEN_WEAK_AXISYM_ORBIT_LOCK
outgoing_edges:
- target: OPEN_WEAK_AXISYM_ORBIT_LOCK
  relation: FEEDS
  status: open
  note: Weak-axisymmetric tangent/orbit-lock packet depends on this signature.
incoming_edges:
- source: EQ_WEAK_AXISYM_SIGNATURE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: CLAIM_STAGE024_O3_ISOTROPY
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_17_WEAK_AXISYM
  relation: FREEZES
  status: exact_first_order
  note: Weak-axisymmetric first-order grouped signature.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_angular_first_order
- topic/moving_throat
- type/anisotropy_fingerprint
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weak-axisymmetric splitting law

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_WEAK_AXISYM_SPLITTING`
> **Status:** `exact_angular_first_order`
> **Layer:** `math_object`
> **Type:** `anisotropy_fingerprint`

## Summary

Pure weak-axisymmetric l=2 perturbation forces grouped signature (20,21,22)~(1,1/2,-1), equivalently b=3a.

## Physical Meaning

Pure weak-axisymmetric l=2 perturbation forces grouped signature (20,21,22)~(1,1/2,-1), equivalently b=3a.

## Mathematical Role

- Layer: `math_object`
- Type: `anisotropy_fingerprint`
- Status: `exact_angular_first_order`

## Equation

$$
lambda20=1
$$

$$
lambda21=1/2
$$

$$
lambda22=-1
$$

$$
b=3a
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_WEAK_AXISYM_SIGNATURE]]

### Related claims
- [[CLAIM_STAGE024_O3_ISOTROPY]]

### Open gates
- [[OPEN_WEAK_AXISYM_ORBIT_LOCK]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[OPEN_WEAK_AXISYM_ORBIT_LOCK]] | Weak-axisymmetric tangent/orbit-lock packet depends on this signature. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_WEAK_AXISYM_SIGNATURE]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE024_O3_ISOTROPY]] | Claim feeds this downstream object, output, or open gate. |
| `FREEZES` | [[MT_V2_17_WEAK_AXISYM]] | Weak-axisymmetric first-order grouped signature. |

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
