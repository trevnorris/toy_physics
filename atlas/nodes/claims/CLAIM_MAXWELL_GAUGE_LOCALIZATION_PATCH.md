---
id: CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
title: Maxwell gauge-fixing localization patch
type: claim
layer: claim_theorem
status: safe_interpretation_or_structural_patch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The unweighted Lorenz term is consistent as a bulk gauge device but singular as a noncompact zero-mode gauge-fixed action; the atlas records either impose Lorenz before reductio...
future_paper_needed: false
source_links:
- '[[FILE_EM_FIELDS]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_PDE_MAXWELL_GAUGE]]'
physical_ids:
- PHYS_LOCALIZED_EM_SECTOR
- PHYS_MIXED_EM_CORE
math_ids:
- MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
equation_ids:
- EQ_ZERO_MODE_MAXWELL
status_firewall_ids:
- FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED
source_ids:
- FILE_EM_FIELDS
- FILE_PDE_AUDIT
- SEC_PDE_MAXWELL_GAUGE
outgoing_edges:
- target: MATH_MAXWELL_WEIGHTED_GAUGE_FIXING
  relation: FEEDS_OR_STATUS_OF
  status: safe_interpretation_or_structural_patch
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_V2_02_MAXWELL_GAUGE_AUDIT
  relation: FEEDS_OR_STATUS_OF
  status: safe_interpretation_or_structural_patch
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_PDE_MAXWELL_GAUGE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Gauge localization audit.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: safe_interpretation_or_structural_patch
  note: 'Claim class: patched_required'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_LOCALIZED_EM_SECTOR
  relation: GROUNDS_PHYSICAL_MEANING
  status: safe_interpretation_or_structural_patch
  note: Physical ontology object grounded by this claim.
- source: PHYS_MIXED_EM_CORE
  relation: GROUNDS_PHYSICAL_MEANING
  status: safe_interpretation_or_structural_patch
  note: Physical ontology object grounded by this claim.
- source: FILE_EM_FIELDS
  relation: OWNS_OR_ANCHORS_CLAIM
  status: safe_interpretation_or_structural_patch
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: safe_interpretation_or_structural_patch
  note: Source artifact anchors this claim.
- source: FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: EQ_ZERO_MODE_MAXWELL
  relation: SUPPORTS_CLAIM
  status: safe_interpretation_or_structural_patch
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/safe_interpretation_or_structural_patch
- topic/maxwell
- topic/moving_throat
- topic/projection
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Maxwell gauge-fixing localization patch

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH`  
> **Status:** `safe_interpretation_or_structural_patch`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

The unweighted Lorenz term is consistent as a bulk gauge device but singular as a noncompact zero-mode gauge-fixed action; the atlas records either impose Lorenz before reduction or use localized H(w), preferably H=Z.

## Claim

The unweighted Lorenz term is consistent as a bulk gauge device but singular as a noncompact zero-mode gauge-fixed action; the atlas records either impose Lorenz before reduction or use localized H(w), preferably H=Z.

## Physical Meaning

The unweighted Lorenz term is consistent as a bulk gauge device but singular as a noncompact zero-mode gauge-fixed action; the atlas records either impose Lorenz before reduction or use localized H(w), preferably H=Z.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `safe_interpretation_or_structural_patch`
- Outputs: `MT_V2_02_MAXWELL_GAUGE_AUDIT`, `MATH_MAXWELL_WEIGHTED_GAUGE_FIXING`

## Atlas Links

### Related physical nodes
- [[PHYS_LOCALIZED_EM_SECTOR]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]]

### Related equations
- [[EQ_ZERO_MODE_MAXWELL]]

### Related claims
- none

### Open gates
- none

### Status firewalls
- [[FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED]]

### Source anchors
- [[FILE_EM_FIELDS]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_MAXWELL_GAUGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_MAXWELL_WEIGHTED_GAUGE_FIXING]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_V2_02_MAXWELL_GAUGE_AUDIT]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_MAXWELL_GAUGE]] | Gauge localization audit. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: patched_required |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_LOCALIZED_EM_SECTOR]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MIXED_EM_CORE]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_EM_FIELDS]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_MAXWELL_GAUGE_PATCH_REQUIRED]] | Firewall preserves this correct status boundary. |
| `SUPPORTS_CLAIM` | [[EQ_ZERO_MODE_MAXWELL]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_EM_FIELDS]]
- [[FILE_PDE_AUDIT]]
- [[SEC_PDE_MAXWELL_GAUGE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
