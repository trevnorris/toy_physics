---
id: MATH_PARENT_ACTION_CURRENT
title: 'Current parent action: GNLS + localized Maxwell'
type: action
layer: math_object
status: exact_declared_parent_with_geometry_argument
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Exact parent action currently includes matter and Maxwell sectors; geometry enters through V_conf(Sigma).
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/pde_audit_full.md
- 4d_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- pde_audit_full.md:V2-01
source_links:
- '[[FILE_4D_PARENT]]'
physical_ids:
- PHYS_BULK_ARENA
math_ids:
- MATH_GNLS_PSI
- MATH_LOCALIZED_MAXWELL_AM
- MATH_SIGMA_R_FIELD
equation_ids:
- EQ_PARENT_ACTION_CURRENT
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
source_ids:
- FILE_4D_PARENT
outgoing_edges:
- target: MT_V2_01_PARENT_WALL_AUDIT
  relation: AUDITED_BY
  status: audit
  note: Audit checks whether Sigma/R is autonomous dynamical field.
- target: MATH_GNLS_PSI
  relation: DERIVES
  status: exact
  note: Variation gives GNLS matter equation and continuity.
- target: MATH_LOCALIZED_MAXWELL_AM
  relation: DERIVES
  status: exact
  note: Variation gives localized Maxwell equation.
- target: MATH_SIGMA_R_FIELD
  relation: USES_AS_COUPLING_ARGUMENT
  status: exact coupling argument
  note: Geometry currently enters through V_conf(Sigma).
incoming_edges:
- source: EQ_PARENT_ACTION_CURRENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_4D_PARENT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: FEEDS_OR_STATUS_OF
  status: exact_parent_current
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_BULK_ARENA
  relation: HOSTS
  status: exact
  note: Parent action lives on the 4+1 bulk arena.
- source: NEG_QUERY_CURRENT_PARENT_HAS_AUTONOMOUS_WALL_PDE
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MATH_PARENT_ACTION_CURRENT.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_declared_parent_with_geometry_argument
- topic/maxwell
- topic/moving_throat
- type/action
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Current parent action: GNLS + localized Maxwell

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_PARENT_ACTION_CURRENT`
> **Status:** `exact_declared_parent_with_geometry_argument`
> **Layer:** `math_object`
> **Type:** `action`

## Summary

Exact parent action currently includes matter and Maxwell sectors; geometry enters through V_conf(Sigma).

## Physical Meaning

Exact parent action currently includes matter and Maxwell sectors; geometry enters through V_conf(Sigma).

## Mathematical Role

- Layer: `math_object`
- Type: `action`
- Status: `exact_declared_parent_with_geometry_argument`

## Equation

$$
S_parent=∫(L_psi+L_EM)
$$

$$
V_conf(X;Sigma)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BULK_ARENA]]

### Related math nodes
- [[MATH_GNLS_PSI]]
- [[MATH_LOCALIZED_MAXWELL_AM]]
- [[MATH_SIGMA_R_FIELD]]

### Related equations
- [[EQ_PARENT_ACTION_CURRENT]]

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `AUDITED_BY` | [[MT_V2_01_PARENT_WALL_AUDIT]] | Audit checks whether Sigma/R is autonomous dynamical field. |
| `DERIVES` | [[MATH_GNLS_PSI]] | Variation gives GNLS matter equation and continuity. |
| `DERIVES` | [[MATH_LOCALIZED_MAXWELL_AM]] | Variation gives localized Maxwell equation. |
| `USES_AS_COUPLING_ARGUMENT` | [[MATH_SIGMA_R_FIELD]] | Geometry currently enters through V_conf(Sigma). |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_PARENT_ACTION_CURRENT]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_4D_PARENT]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Claim feeds this downstream object, output, or open gate. |
| `HOSTS` | [[PHYS_BULK_ARENA]] | Parent action lives on the 4+1 bulk arena. |
| `STARTS_AT` | [[NEG_QUERY_CURRENT_PARENT_HAS_AUTONOMOUS_WALL_PDE]] | Negative query starts from MATH_PARENT_ACTION_CURRENT. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]

### Source files
- `research/4d/paper/4d.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
