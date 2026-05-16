---
id: MATH_LOCALIZED_MAXWELL_AM
title: Localized 4+1 Maxwell field
type: field_equation
layer: math_object
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Gauge field sector with localization weight and source consistency.
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
- notes/pde_audit_full.md
- 4d_em_fields_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_em_fields_summary.md
- pde_audit_full.md:V2-02
source_links:
- '[[FILE_EM_FIELDS]]'
physical_ids:
- PHYS_LOCALIZED_EM_SECTOR
math_ids:
- MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
- MATH_MIXED_FIELDS_EW_CA
- MATH_PARENT_ACTION_CURRENT
- MATH_ZERO_MODE_MAXWELL
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
source_ids:
- FILE_EM_FIELDS
outgoing_edges:
- target: MT_V2_02_MAXWELL_GAUGE_AUDIT
  relation: AUDITED_BY
  status: audit
  note: Gauge localization and gauge-fixing conventions audited.
- target: MATH_MIXED_FIELDS_EW_CA
  relation: CONTAINS
  status: exact
  note: Mixed gauge-invariant fields remain in parent ontology.
- target: MATH_ZERO_MODE_MAXWELL
  relation: REDUCES_TO
  status: controlled
  note: Axial gauge/zero-mode/far-field brane assumptions give 3+1 Maxwell.
incoming_edges:
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_LOCALIZED_MAXWELL_AM.
- source: MATH_PARENT_ACTION_CURRENT
  relation: DERIVES
  status: exact
  note: Variation gives localized Maxwell equation.
- source: FILE_EM_FIELDS
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: FEEDS_OR_STATUS_OF
  status: exact_parent_current
  note: Claim feeds this downstream object, output, or open gate.
- source: MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
  relation: PATCHES_OR_INTERPRETS
  status: patch_or_safe_interpretation
  note: Keeps localized Maxwell sector safe under chosen gauge-fixing convention.
- source: PHYS_LOCALIZED_EM_SECTOR
  relation: REPRESENTED_BY
  status: exact
  note: Localized EM sector represented by A_M,F_MN with Z(w).
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact
- topic/maxwell
- topic/moving_throat
- type/field_equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Localized 4+1 Maxwell field

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_LOCALIZED_MAXWELL_AM`
> **Status:** `exact`
> **Layer:** `math_object`
> **Type:** `field_equation`

## Summary

Gauge field sector with localization weight and source consistency.

## Physical Meaning

Gauge field sector with localization weight and source consistency.

## Mathematical Role

- Layer: `math_object`
- Type: `field_equation`
- Status: `exact`

## Equation

$$
∂_M(ZF^{MN})+(1/ξ)∂^N(∂·A)=μ0J^N
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_LOCALIZED_EM_SECTOR]]

### Related math nodes
- [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]]
- [[MATH_MIXED_FIELDS_EW_CA]]
- [[MATH_PARENT_ACTION_CURRENT]]
- [[MATH_ZERO_MODE_MAXWELL]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_EM_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `AUDITED_BY` | [[MT_V2_02_MAXWELL_GAUGE_AUDIT]] | Gauge localization and gauge-fixing conventions audited. |
| `CONTAINS` | [[MATH_MIXED_FIELDS_EW_CA]] | Mixed gauge-invariant fields remain in parent ontology. |
| `REDUCES_TO` | [[MATH_ZERO_MODE_MAXWELL]] | Axial gauge/zero-mode/far-field brane assumptions give 3+1 Maxwell. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references MATH_LOCALIZED_MAXWELL_AM. |
| `DERIVES` | [[MATH_PARENT_ACTION_CURRENT]] | Variation gives localized Maxwell equation. |
| `DOCUMENTS` | [[FILE_EM_FIELDS]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Claim feeds this downstream object, output, or open gate. |
| `PATCHES_OR_INTERPRETS` | [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]] | Keeps localized Maxwell sector safe under chosen gauge-fixing convention. |
| `REPRESENTED_BY` | [[PHYS_LOCALIZED_EM_SECTOR]] | Localized EM sector represented by A_M,F_MN with Z(w). |

## Source Anchors

### Source anchor notes
- [[FILE_EM_FIELDS]]

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`
- `notes/pde_audit_full.md`
- `4d_em_fields_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
