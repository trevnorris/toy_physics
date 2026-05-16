---
id: CLAIM_PARENT_ACTION_CURRENT_EXACT
title: Current parent action and exact bulk equations
type: claim
layer: claim_theorem
status: exact_parent_current
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The declared 4+1 parent stack supplies gauged GNLS matter, localized Maxwell, exact bulk current/continuity, and mixed-sector observables, with Sigma/R currently entering throug...
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_EM_FIELDS]]'
- '[[FILE_PLASMA]]'
- '[[SEC_4D_EXACT_EQUATIONS]]'
- '[[SEC_4D_PARENT_ACTION]]'
- '[[SEC_EM_BULK_EQUATIONS]]'
- '[[SEC_EM_LOCALIZED_ACTION]]'
- '[[SEC_MT_PARENT_FIELDS]]'
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
- PHYS_BULK_ARENA
- PHYS_LOCALIZED_EM_SECTOR
- PHYS_SUPERFLUID_MEDIUM
math_ids:
- MATH_GNLS_PSI
- MATH_LOCALIZED_MAXWELL_AM
- MATH_PARENT_ACTION_CURRENT
equation_ids:
- EQ_BULK_CONTINUITY
- EQ_GNLS_BULK
- EQ_PARENT_ACTION_CURRENT
- EQ_ZERO_MODE_MAXWELL
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
- CLAIM_STAGE1_GEOMETRY_LIFT
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
source_ids:
- FILE_4D_PARENT
- FILE_EM_FIELDS
- FILE_PLASMA
- SEC_4D_EXACT_EQUATIONS
- SEC_4D_PARENT_ACTION
- SEC_EM_BULK_EQUATIONS
- SEC_EM_LOCALIZED_ACTION
- SEC_MT_PARENT_FIELDS
outgoing_edges:
- target: MATH_GNLS_PSI
  relation: FEEDS_OR_STATUS_OF
  status: exact_parent_current
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_LOCALIZED_MAXWELL_AM
  relation: FEEDS_OR_STATUS_OF
  status: exact_parent_current
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_PARENT_ACTION_CURRENT
  relation: FEEDS_OR_STATUS_OF
  status: exact_parent_current
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: SUPPLIES_PARENT_CONTEXT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: SUPPLIES_PARENT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: SUPPLIES_PARENT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_4D_EXACT_EQUATIONS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Exact E-L and hydrodynamic identities.
- source: SEC_4D_PARENT_ACTION
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent action and sectors.
- source: SEC_EM_BULK_EQUATIONS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Localized Maxwell E-L.
- source: SEC_EM_LOCALIZED_ACTION
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Maxwell action anchor.
- source: SEC_MT_PARENT_FIELDS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Parent/moving-throat fields.
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_PARENT_ACTION_CURRENT_EXACT.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_parent_current
  note: 'Claim class: exact'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_BULK_ARENA
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_parent_current
  note: Physical ontology object grounded by this claim.
- source: PHYS_LOCALIZED_EM_SECTOR
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_parent_current
  note: Physical ontology object grounded by this claim.
- source: PHYS_SUPERFLUID_MEDIUM
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_parent_current
  note: Physical ontology object grounded by this claim.
- source: FILE_4D_PARENT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_parent_current
  note: Source artifact anchors this claim.
- source: FILE_EM_FIELDS
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_parent_current
  note: Source artifact anchors this claim.
- source: FILE_PLASMA
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_parent_current
  note: Source artifact anchors this claim.
- source: EQ_BULK_CONTINUITY
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
- source: EQ_GNLS_BULK
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
- source: EQ_PARENT_ACTION_CURRENT
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
- source: EQ_ZERO_MODE_MAXWELL
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_parent_current
- topic/maxwell
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Current parent action and exact bulk equations

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_PARENT_ACTION_CURRENT_EXACT`
> **Status:** `exact_parent_current`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The declared 4+1 parent stack supplies gauged GNLS matter, localized Maxwell, exact bulk current/continuity, and mixed-sector observables, with Sigma/R currently entering through confinement.

## Claim

The declared 4+1 parent stack supplies gauged GNLS matter, localized Maxwell, exact bulk current/continuity, and mixed-sector observables, with Sigma/R currently entering through confinement.

## Physical Meaning

The declared 4+1 parent stack supplies gauged GNLS matter, localized Maxwell, exact bulk current/continuity, and mixed-sector observables, with Sigma/R currently entering through confinement.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_parent_current`
- Outputs: `MATH_PARENT_ACTION_CURRENT`, `MATH_GNLS_PSI`, `MATH_LOCALIZED_MAXWELL_AM`

## Atlas Links

### Related physical nodes
- [[PHYS_BULK_ARENA]]
- [[PHYS_LOCALIZED_EM_SECTOR]]
- [[PHYS_SUPERFLUID_MEDIUM]]

### Related math nodes
- [[MATH_GNLS_PSI]]
- [[MATH_LOCALIZED_MAXWELL_AM]]
- [[MATH_PARENT_ACTION_CURRENT]]

### Related equations
- [[EQ_BULK_CONTINUITY]]
- [[EQ_GNLS_BULK]]
- [[EQ_PARENT_ACTION_CURRENT]]
- [[EQ_ZERO_MODE_MAXWELL]]

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]
- [[FILE_PLASMA]]
- [[SEC_4D_EXACT_EQUATIONS]]
- [[SEC_4D_PARENT_ACTION]]
- [[SEC_EM_BULK_EQUATIONS]]
- [[SEC_EM_LOCALIZED_ACTION]]
- [[SEC_MT_PARENT_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_GNLS_PSI]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_LOCALIZED_MAXWELL_AM]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_PARENT_ACTION_CURRENT]] | Claim feeds this downstream object, output, or open gate. |
| `SUPPLIES_PARENT_CONTEXT_FOR` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Claim-level dependency added in v0.4. |
| `SUPPLIES_PARENT_FOR` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Claim-level dependency added in v0.4. |
| `SUPPLIES_PARENT_FOR` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_EXACT_EQUATIONS]] | Exact E-L and hydrodynamic identities. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_PARENT_ACTION]] | Parent action and sectors. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_BULK_EQUATIONS]] | Localized Maxwell E-L. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_LOCALIZED_ACTION]] | Maxwell action anchor. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_PARENT_FIELDS]] | Parent/moving-throat fields. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references CLAIM_PARENT_ACTION_CURRENT_EXACT. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_BULK_ARENA]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_LOCALIZED_EM_SECTOR]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_SUPERFLUID_MEDIUM]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_4D_PARENT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_EM_FIELDS]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PLASMA]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_BULK_CONTINUITY]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_GNLS_BULK]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_PARENT_ACTION_CURRENT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_ZERO_MODE_MAXWELL]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]
- [[FILE_PLASMA]]
- [[SEC_4D_EXACT_EQUATIONS]]
- [[SEC_4D_PARENT_ACTION]]
- [[SEC_EM_BULK_EQUATIONS]]
- [[SEC_EM_LOCALIZED_ACTION]]
- [[SEC_MT_PARENT_FIELDS]]

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
