---
id: MT_V2_02_MAXWELL_GAUGE_AUDIT
title: V2-02 Maxwell gauge-localization audit
type: audit_gate
layer: status_audit
status: patch_required_or_safe_interpretation
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Audits localized Maxwell kinetic term versus gauge-fixing weight and downstream mixed-sector safety.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-02
math_ids:
- MATH_LOCALIZED_MAXWELL_AM
- MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
claim_ids:
- CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
outgoing_edges:
- target: MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
  relation: AUDITS
  status: status
  note: Weighted gauge-fixing clarifies H=1 versus H=Z interpretations.
incoming_edges:
- source: MATH_LOCALIZED_MAXWELL_AM
  relation: AUDITED_BY
  status: audit
  note: Gauge localization and gauge-fixing conventions audited.
- source: CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
  relation: FEEDS_OR_STATUS_OF
  status: safe_interpretation_or_structural_patch
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/patch_required_or_safe_interpretation
- topic/maxwell
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-02 Maxwell gauge-localization audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_02_MAXWELL_GAUGE_AUDIT`  
> **Status:** `patch_required_or_safe_interpretation`  
> **Layer:** `status_audit`  
> **Type:** `audit_gate`

## Summary

Audits localized Maxwell kinetic term versus gauge-fixing weight and downstream mixed-sector safety.

## Physical Meaning

Audits localized Maxwell kinetic term versus gauge-fixing weight and downstream mixed-sector safety.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `patch_required_or_safe_interpretation`

## Equation

$$
∂_M(ZF^{MN})+(1/ξ)∂^N(HB)=μ0J^N
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_LOCALIZED_MAXWELL_AM]]
- [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]]

### Related equations
- none

### Related claims
- [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `AUDITS` | [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]] | Weighted gauge-fixing clarifies H=1 versus H=Z interpretations. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `AUDITED_BY` | [[MATH_LOCALIZED_MAXWELL_AM]] | Gauge localization and gauge-fixing conventions audited. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
