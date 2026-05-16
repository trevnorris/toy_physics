---
id: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
title: Zero-mode localized Maxwell gives brane Maxwell only as a controlled limit
type: claim
layer: claim_theorem
status: controlled_reduction
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The localized 4+1 Maxwell action reduces to ordinary 3+1 Maxwell with mu0_eff=mu0/Z_int only under axial/zero-mode/far-field assumptions; mixed core channels are suppressed, not...
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_EM_FIELDS]]'
- '[[FILE_PLASMA]]'
- '[[SEC_4D_ZERO_MODE_MAXWELL]]'
- '[[SEC_EM_LOCALIZED_ACTION]]'
- '[[SEC_EM_QEFF]]'
- '[[SEC_EM_ZERO_MODE]]'
- '[[SEC_MT_PROJECTION_ZERO_MODE]]'
- '[[SEC_PLASMA_Z_VS_W]]'
tex_anchor:
  file: research/4d_em_fields/paper/4d_em_fields.tex
  line: 245
  heading_level: section
  heading: Conventions and localized 5D Maxwell action
  nearest_label:
    name: sec:conventions_action
    line: 246
  nearby_labels:
  - name: sec:conventions_action
    line: 246
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 0.534
  confidence: medium
  source_anchor_node: SEC_EM_LOCALIZED_ACTION
physical_ids:
- PHYS_BRANE_OBSERVER
- PHYS_LOCALIZED_EM_SECTOR
- PHYS_MIXED_EM_CORE
math_ids:
- MATH_ZERO_MODE_MAXWELL
- MATH_Z_PROFILE
equation_ids:
- EQ_QEFF_NORMALIZATION
- EQ_ZERO_MODE_MAXWELL
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
- CLAIM_PARENT_ACTION_CURRENT_EXACT
status_firewall_ids:
- FIREWALL_PROJECTION_NOT_REDUCTION
- FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE
source_ids:
- FILE_4D_PARENT
- FILE_EM_FIELDS
- FILE_PLASMA
- SEC_4D_ZERO_MODE_MAXWELL
- SEC_EM_LOCALIZED_ACTION
- SEC_EM_QEFF
- SEC_EM_ZERO_MODE
- SEC_MT_PROJECTION_ZERO_MODE
- SEC_PLASMA_Z_VS_W
outgoing_edges:
- target: ATOM_HYDROGEN_ZERO_MODE
  relation: FEEDS_OR_STATUS_OF
  status: controlled_reduction
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_ZERO_MODE_MAXWELL
  relation: FEEDS_OR_STATUS_OF
  status: controlled_reduction
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_Z_PROFILE
  relation: FEEDS_OR_STATUS_OF
  status: controlled_reduction
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: SUPPORTS_REDUCED_SECTOR
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_4D_ZERO_MODE_MAXWELL
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Controlled EM zero-mode reduction.
- source: SEC_EM_LOCALIZED_ACTION
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Maxwell action anchor.
- source: SEC_EM_QEFF
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Thickness-controlled charge.
- source: SEC_EM_ZERO_MODE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Zero-mode Maxwell limit.
- source: SEC_MT_PROJECTION_ZERO_MODE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Projection and zero-mode hooks.
- source: SEC_PLASMA_Z_VS_W
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Z controls action; W controls observation.
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_ZERO_MODE_MAXWELL_REDUCTION.
- source: BACKLINK_ATOM_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_ZERO_MODE_MAXWELL_REDUCTION.
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_ZERO_MODE_MAXWELL_REDUCTION.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: controlled_reduction
  note: 'Claim class: controlled_reduction'
- source: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: CONSTRAINS_NOTATION_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_BRANE_OBSERVER
  relation: GROUNDS_PHYSICAL_MEANING
  status: controlled_reduction
  note: Physical ontology object grounded by this claim.
- source: PHYS_LOCALIZED_EM_SECTOR
  relation: GROUNDS_PHYSICAL_MEANING
  status: controlled_reduction
  note: Physical ontology object grounded by this claim.
- source: PHYS_MIXED_EM_CORE
  relation: GROUNDS_PHYSICAL_MEANING
  status: controlled_reduction
  note: Physical ontology object grounded by this claim.
- source: FILE_4D_PARENT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: controlled_reduction
  note: Source artifact anchors this claim.
- source: FILE_EM_FIELDS
  relation: OWNS_OR_ANCHORS_CLAIM
  status: controlled_reduction
  note: Source artifact anchors this claim.
- source: FILE_PLASMA
  relation: OWNS_OR_ANCHORS_CLAIM
  status: controlled_reduction
  note: Source artifact anchors this claim.
- source: FIREWALL_PROJECTION_NOT_REDUCTION
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_ZERO_MODE_ERASES_MIXED_CORE
  relation: STARTS_AT
  status: v07
  note: Negative query starts from CLAIM_ZERO_MODE_MAXWELL_REDUCTION.
- source: QUERY_ZERO_MODE_DOWNSTREAM
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
- source: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPLIES_PARENT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_QEFF_NORMALIZATION
  relation: SUPPORTS_CLAIM
  status: controlled_reduction
  note: Equation anchor supports this named claim.
- source: EQ_ZERO_MODE_MAXWELL
  relation: SUPPORTS_CLAIM
  status: controlled_reduction
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/controlled_reduction
- topic/atom
- topic/maxwell
- topic/projection
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Zero-mode localized Maxwell gives brane Maxwell only as a controlled limit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_ZERO_MODE_MAXWELL_REDUCTION`
> **Status:** `controlled_reduction`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The localized 4+1 Maxwell action reduces to ordinary 3+1 Maxwell with mu0_eff=mu0/Z_int only under axial/zero-mode/far-field assumptions; mixed core channels are suppressed, not erased.

## Claim

The localized 4+1 Maxwell action reduces to ordinary 3+1 Maxwell with mu0_eff=mu0/Z_int only under axial/zero-mode/far-field assumptions; mixed core channels are suppressed, not erased.

## Physical Meaning

The localized 4+1 Maxwell action reduces to ordinary 3+1 Maxwell with mu0_eff=mu0/Z_int only under axial/zero-mode/far-field assumptions; mixed core channels are suppressed, not erased.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `controlled_reduction`
- Outputs: `MATH_ZERO_MODE_MAXWELL`, `MATH_Z_PROFILE`, `ATOM_HYDROGEN_ZERO_MODE`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_OBSERVER]]
- [[PHYS_LOCALIZED_EM_SECTOR]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- [[MATH_ZERO_MODE_MAXWELL]]
- [[MATH_Z_PROFILE]]

### Related equations
- [[EQ_QEFF_NORMALIZATION]]
- [[EQ_ZERO_MODE_MAXWELL]]

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_PROJECTION_NOT_REDUCTION]]
- [[FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE]]

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]
- [[FILE_PLASMA]]
- [[SEC_4D_ZERO_MODE_MAXWELL]]
- [[SEC_EM_LOCALIZED_ACTION]]
- [[SEC_EM_QEFF]]
- [[SEC_EM_ZERO_MODE]]
- [[SEC_MT_PROJECTION_ZERO_MODE]]
- [[SEC_PLASMA_Z_VS_W]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[ATOM_HYDROGEN_ZERO_MODE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_ZERO_MODE_MAXWELL]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_Z_PROFILE]] | Claim feeds this downstream object, output, or open gate. |
| `SUPPORTS_REDUCED_SECTOR` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_ZERO_MODE_MAXWELL]] | Controlled EM zero-mode reduction. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_LOCALIZED_ACTION]] | Maxwell action anchor. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_QEFF]] | Thickness-controlled charge. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_ZERO_MODE]] | Zero-mode Maxwell limit. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_PROJECTION_ZERO_MODE]] | Projection and zero-mode hooks. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PLASMA_Z_VS_W]] | Z controls action; W controls observation. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references CLAIM_ZERO_MODE_MAXWELL_REDUCTION. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_ATOM_WORK]] | Paper backlink block references CLAIM_ZERO_MODE_MAXWELL_REDUCTION. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references CLAIM_ZERO_MODE_MAXWELL_REDUCTION. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: controlled_reduction |
| `CONSTRAINS_NOTATION_FOR` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Claim-level dependency added in v0.4. |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_BRANE_OBSERVER]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_LOCALIZED_EM_SECTOR]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MIXED_EM_CORE]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_4D_PARENT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_EM_FIELDS]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PLASMA]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_PROJECTION_NOT_REDUCTION]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_ZERO_MODE_ERASES_MIXED_CORE]] | Negative query starts from CLAIM_ZERO_MODE_MAXWELL_REDUCTION. |
| `STARTS_AT` | [[QUERY_ZERO_MODE_DOWNSTREAM]] | Query validation start node. |
| `SUPPLIES_PARENT_FOR` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_QEFF_NORMALIZATION]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_ZERO_MODE_MAXWELL]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]
- [[FILE_PLASMA]]
- [[SEC_4D_ZERO_MODE_MAXWELL]]
- [[SEC_EM_LOCALIZED_ACTION]]
- [[SEC_EM_QEFF]]
- [[SEC_EM_ZERO_MODE]]
- [[SEC_MT_PROJECTION_ZERO_MODE]]
- [[SEC_PLASMA_Z_VS_W]]

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`

### TeX anchor
- File: `research/4d_em_fields/paper/4d_em_fields.tex`
- Line: `245`
- Heading: Conventions and localized 5D Maxwell action
- Nearest label: `sec:conventions_action` at line `246`
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
