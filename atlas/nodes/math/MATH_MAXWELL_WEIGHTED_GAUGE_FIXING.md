---
id: MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
title: Weighted Maxwell gauge fixing
type: gauge_fixing_audit
layer: math_object
status: safe_interpretation_or_patch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: General H(w)-weighted gauge-fixing audit for localized Maxwell theory; H=1 is a bulk gauge device, H=Z is localized gauge fixing.
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
- MATH_LOCALIZED_MAXWELL_AM
claim_ids:
- CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
status_firewall_ids:
- FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED
outgoing_edges:
- target: MATH_LOCALIZED_MAXWELL_AM
  relation: PATCHES_OR_INTERPRETS
  status: patch_or_safe_interpretation
  note: Keeps localized Maxwell sector safe under chosen gauge-fixing convention.
incoming_edges:
- source: MT_V2_02_MAXWELL_GAUGE_AUDIT
  relation: AUDITS
  status: status
  note: Weighted gauge-fixing clarifies H=1 versus H=Z interpretations.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_MAXWELL_WEIGHTED_GAUGE_FIXING.
- source: CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
  relation: FEEDS_OR_STATUS_OF
  status: safe_interpretation_or_structural_patch
  note: Claim feeds this downstream object, output, or open gate.
- source: FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_GAUGE_FIXING_LOCALIZATION_IGNORED
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MATH_MAXWELL_WEIGHTED_GAUGE_FIXING.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/safe_interpretation_or_patch
- topic/maxwell
- topic/moving_throat
- type/gauge_fixing_audit
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weighted Maxwell gauge fixing

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_MAXWELL_WEIGHTED_GAUGE_FIXING`  
> **Status:** `safe_interpretation_or_patch`  
> **Layer:** `math_object`  
> **Type:** `gauge_fixing_audit`

## Summary

General H(w)-weighted gauge-fixing audit for localized Maxwell theory; H=1 is a bulk gauge device, H=Z is localized gauge fixing.

## Physical Meaning

General H(w)-weighted gauge-fixing audit for localized Maxwell theory; H=1 is a bulk gauge device, H=Z is localized gauge fixing.

## Mathematical Role

- Layer: `math_object`
- Type: `gauge_fixing_audit`
- Status: `safe_interpretation_or_patch`

## Equation

$$
d_M(ZF^MN)+(1/xi)d^N(HB)=mu0 J^N
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_LOCALIZED_MAXWELL_AM]]

### Related equations
- none

### Related claims
- [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `PATCHES_OR_INTERPRETS` | [[MATH_LOCALIZED_MAXWELL_AM]] | Keeps localized Maxwell sector safe under chosen gauge-fixing convention. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `AUDITS` | [[MT_V2_02_MAXWELL_GAUGE_AUDIT]] | Weighted gauge-fixing clarifies H=1 versus H=Z interpretations. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references MATH_MAXWELL_WEIGHTED_GAUGE_FIXING. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]] | Claim feeds this downstream object, output, or open gate. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_GAUGE_FIXING_LOCALIZATION_IGNORED]] | Negative query starts from MATH_MAXWELL_WEIGHTED_GAUGE_FIXING. |

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
