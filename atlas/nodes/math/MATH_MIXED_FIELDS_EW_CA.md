---
id: MATH_MIXED_FIELDS_EW_CA
title: Mixed gauge invariants
type: gauge_invariant_fields
layer: math_object
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Gauge-invariant mixed slots used for brane-bulk EM exchange and outgoing bridge.
future_paper_needed: false
source_files:
- research/4d_plasma/paper/4d_plasma.tex
- notes/pde_audit_full.md
- 4d_plasma_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_plasma_summary.md
- pde_audit_full.md:V2-30
source_links:
- '[[FILE_EM_FIELDS]]'
physical_ids:
- PHYS_MIXED_EM_CORE
math_ids:
- MATH_LOCALIZED_MAXWELL_AM
- MATH_MAXWELL_MIXED_KERNEL
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
source_ids:
- FILE_EM_FIELDS
outgoing_edges:
- target: MT_STAGE4_MAXWELL_MIXED
  relation: ENABLES
  status: exact/reduced
  note: Mixed fields are the microscopic place for outgoing bridge.
- target: MATH_MAXWELL_MIXED_KERNEL
  relation: REDUCES_TO
  status: reduced
  note: Mixed fields feed reduced port-active Maxwell/mixed kernel.
incoming_edges:
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_MIXED_FIELDS_EW_CA.
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_MIXED_FIELDS_EW_CA.
- source: MATH_LOCALIZED_MAXWELL_AM
  relation: CONTAINS
  status: exact
  note: Mixed gauge-invariant fields remain in parent ontology.
- source: FILE_EM_FIELDS
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_MIXED_EM_CORE
  relation: REPRESENTED_BY
  status: exact
  note: Mixed physical channels represented by gauge-invariant E_w,C_a and A_w/F_mu_w/J^w.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact
- topic/maxwell
- topic/moving_throat
- type/gauge_invariant_fields
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Mixed gauge invariants

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_MIXED_FIELDS_EW_CA`  
> **Status:** `exact`  
> **Layer:** `math_object`  
> **Type:** `gauge_invariant_fields`

## Summary

Gauge-invariant mixed slots used for brane-bulk EM exchange and outgoing bridge.

## Physical Meaning

Gauge-invariant mixed slots used for brane-bulk EM exchange and outgoing bridge.

## Mathematical Role

- Layer: `math_object`
- Type: `gauge_invariant_fields`
- Status: `exact`

## Equation

$$
E_w=-∂_t A_w-∂_w A_0
$$

$$
C_a=∂_a A_w-∂_w A_a
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- [[MATH_LOCALIZED_MAXWELL_AM]]
- [[MATH_MAXWELL_MIXED_KERNEL]]

### Related equations
- none

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_EM_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ENABLES` | [[MT_STAGE4_MAXWELL_MIXED]] | Mixed fields are the microscopic place for outgoing bridge. |
| `REDUCES_TO` | [[MATH_MAXWELL_MIXED_KERNEL]] | Mixed fields feed reduced port-active Maxwell/mixed kernel. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references MATH_MIXED_FIELDS_EW_CA. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references MATH_MIXED_FIELDS_EW_CA. |
| `CONTAINS` | [[MATH_LOCALIZED_MAXWELL_AM]] | Mixed gauge-invariant fields remain in parent ontology. |
| `DOCUMENTS` | [[FILE_EM_FIELDS]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim feeds this downstream object, output, or open gate. |
| `REPRESENTED_BY` | [[PHYS_MIXED_EM_CORE]] | Mixed physical channels represented by gauge-invariant E_w,C_a and A_w/F_mu_w/J^w. |

## Source Anchors

### Source anchor notes
- [[FILE_EM_FIELDS]]

### Source files
- `research/4d_plasma/paper/4d_plasma.tex`
- `notes/pde_audit_full.md`
- `4d_plasma_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
