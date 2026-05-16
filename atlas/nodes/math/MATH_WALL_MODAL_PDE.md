---
id: MATH_WALL_MODAL_PDE
title: Wall modal PDE
type: modal_equation
layer: math_object
status: exact_if_S_eta_included
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The l-mode equation supplied by the quadratic wall action, with l=0 scalar and l=2 grouped P2 specialization.
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
- PHYS_GROUPED_P2_SHAPE_RESPONSE
math_ids:
- MATH_WALL_BOUNDARY_TRACTION
- MATH_WALL_ENERGY_POSITIVITY
equation_ids:
- EQ_WALL_MODAL_PDE
outgoing_edges:
- target: MATH_WALL_BOUNDARY_TRACTION
  relation: HAS_BOUNDARY_TERM
  status: exact_if_closure
  note: Boundary terms define driven/free/Dirichlet port conditions.
- target: MATH_WALL_ENERGY_POSITIVITY
  relation: IMPLIES_GATE
  status: exact_if_closure
  note: Wall PDE has local positivity gates.
- target: MT_STAGE2_BREATHING_REDUCTION
  relation: REDUCES_TO
  status: exact_if_closure
  note: Axisymmetric two-profile reduction recovers old a,L closure.
- target: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: SUPPORTS
  status: exact_if_closure
  note: l=2 modes are literal wall/support modes.
incoming_edges:
- source: EQ_WALL_MODAL_PDE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_LEPTON_MASS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_WALL_MODAL_PDE.
- source: MT_V2_01_PARENT_WALL_AUDIT
  relation: VALIDATES_IF_INCLUDED
  status: effective_closure_pass
  note: Quadratic wall action supplies modal wall PDE.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_if_s_eta_included
- topic/moving_throat
- type/modal_equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Wall modal PDE

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_WALL_MODAL_PDE`
> **Status:** `exact_if_S_eta_included`
> **Layer:** `math_object`
> **Type:** `modal_equation`

## Summary

The l-mode equation supplied by the quadratic wall action, with l=0 scalar and l=2 grouped P2 specialization.

## Physical Meaning

The l-mode equation supplied by the quadratic wall action, with l=0 scalar and l=2 grouped P2 specialization.

## Mathematical Role

- Layer: `math_object`
- Type: `modal_equation`
- Status: `exact_if_S_eta_included`

## Equation

$$
mu_eta q_tt - d_w(T_w q_w) + (K_eta+l(l+1)T_Omega)q = S_lm
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]

### Related math nodes
- [[MATH_WALL_BOUNDARY_TRACTION]]
- [[MATH_WALL_ENERGY_POSITIVITY]]

### Related equations
- [[EQ_WALL_MODAL_PDE]]

### Related claims
- none

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `HAS_BOUNDARY_TERM` | [[MATH_WALL_BOUNDARY_TRACTION]] | Boundary terms define driven/free/Dirichlet port conditions. |
| `IMPLIES_GATE` | [[MATH_WALL_ENERGY_POSITIVITY]] | Wall PDE has local positivity gates. |
| `REDUCES_TO` | [[MT_STAGE2_BREATHING_REDUCTION]] | Axisymmetric two-profile reduction recovers old a,L closure. |
| `SUPPORTS` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | l=2 modes are literal wall/support modes. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_WALL_MODAL_PDE]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_MASS]] | Paper backlink block references MATH_WALL_MODAL_PDE. |
| `VALIDATES_IF_INCLUDED` | [[MT_V2_01_PARENT_WALL_AUDIT]] | Quadratic wall action supplies modal wall PDE. |

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
