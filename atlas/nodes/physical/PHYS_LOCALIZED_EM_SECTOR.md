---
id: PHYS_LOCALIZED_EM_SECTOR
title: Localized electromagnetic sector
type: field_sector
layer: physical_ontology
status: exact_declared_parent
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: A genuine localized 4+1 Maxwell field with transverse localization weight.
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
- notes/pde_audit_full.md
- 4d_em_fields_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_em_fields_summary.md
- pde_audit_full.md:V2-30
math_ids:
- MATH_LOCALIZED_MAXWELL_AM
- MATH_Z_PROFILE
claim_ids:
- CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
- CLAIM_PARENT_ACTION_CURRENT_EXACT
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
outgoing_edges:
- target: CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: safe_interpretation_or_structural_patch
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_parent_current
  note: Physical ontology object grounded by this claim.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: GROUNDS_PHYSICAL_MEANING
  status: controlled_reduction
  note: Physical ontology object grounded by this claim.
- target: MATH_Z_PROFILE
  relation: LOCALIZED_BY
  status: exact
  note: Z controls EM localization and effective coupling.
- target: MATH_LOCALIZED_MAXWELL_AM
  relation: REPRESENTED_BY
  status: exact
  note: Localized EM sector represented by A_M,F_MN with Z(w).
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_declared_parent
- topic/maxwell
- topic/moving_throat
- type/field_sector
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Localized electromagnetic sector

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_LOCALIZED_EM_SECTOR`
> **Status:** `exact_declared_parent`
> **Layer:** `physical_ontology`
> **Type:** `field_sector`

## Summary

A genuine localized 4+1 Maxwell field with transverse localization weight.

## Physical Meaning

A genuine localized 4+1 Maxwell field with transverse localization weight.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `field_sector`
- Status: `exact_declared_parent`

## Equation

$$
A_M
$$

$$
F_MN
$$

$$
Z(w)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_LOCALIZED_MAXWELL_AM]]
- [[MATH_Z_PROFILE]]

### Related equations
- none

### Related claims
- [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]]
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Physical ontology object grounded by this claim. |
| `LOCALIZED_BY` | [[MATH_Z_PROFILE]] | Z controls EM localization and effective coupling. |
| `REPRESENTED_BY` | [[MATH_LOCALIZED_MAXWELL_AM]] | Localized EM sector represented by A_M,F_MN with Z(w). |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`
- `notes/pde_audit_full.md`
- `4d_em_fields_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
