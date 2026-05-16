---
id: CLAIM_2PN_ADM_WITHIN_HIERARCHY
title: 2PN ADM equality within declared hierarchy
type: claim
layer: claim_theorem
status: exact_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The conservative 2PN two-body ledger Legendre-transforms to the ADM target within the declared hierarchy; it also exposes dipole, P0/P2 support, and geometry-closure structure.
future_paper_needed: false
source_files:
- research/4d_2pn/paper/4d_2pn.tex
source_links:
- '[[FILE_2PN_FULL]]'
- '[[SEC_2PN_CLAIMS]]'
- '[[SEC_2PN_MAIN_THEOREM]]'
- '[[SEC_2PN_OPEN_DYNAMIC]]'
- '[[SEC_2PN_P0P2_APPENDIX]]'
tex_anchor:
  file: research/4d_2pn/paper/4d_2pn.tex
  line: 179
  heading_level: subsection
  heading: What this paper claims
  nearest_label:
    name: sec:intro-claims
    line: 180
  nearby_labels:
  - name: sec:intro-claims
    line: 180
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 0.855
  confidence: medium
  source_anchor_node: SEC_2PN_CLAIMS
physical_ids:
- PHYS_GROUPED_P2_SHAPE_RESPONSE
- PHYS_INTERIOR_SUPPORT
equation_ids:
- EQ_2PN_ADM_EQUALITY
claim_ids:
- CLAIM_1PN_EIH_WITHIN_HIERARCHY
- CLAIM_25PN_QUAD_NARROWING
source_ids:
- FILE_2PN_FULL
- SEC_2PN_CLAIMS
- SEC_2PN_MAIN_THEOREM
- SEC_2PN_OPEN_DYNAMIC
- SEC_2PN_P0P2_APPENDIX
outgoing_edges:
- target: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
- target: PN_2_ADM_FULL
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_25PN_QUAD_NARROWING
  relation: STRUCTURAL_INPUT_TO
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_2PN_CLAIMS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: 2PN claims and non-claims.
- source: SEC_2PN_MAIN_THEOREM
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: 2PN exact ADM match.
- source: SEC_2PN_OPEN_DYNAMIC
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Open moving-throat dynamics.
- source: SEC_2PN_P0P2_APPENDIX
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: P0/P2 response structure.
- source: BACKLINK_2PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_2PN_ADM_WITHIN_HIERARCHY.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_assembly_within_closure
  note: 'Claim class: exact_within_closure'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_assembly_within_closure
  note: Physical ontology object grounded by this claim.
- source: PHYS_INTERIOR_SUPPORT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_assembly_within_closure
  note: Physical ontology object grounded by this claim.
- source: CLAIM_1PN_EIH_WITHIN_HIERARCHY
  relation: LOWER_ORDER_ANCHOR_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: FILE_2PN_FULL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_assembly_within_closure
  note: Source artifact anchors this claim.
- source: EQ_2PN_ADM_EQUALITY
  relation: SUPPORTS_CLAIM
  status: exact_assembly_within_closure
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_assembly_within_closure
- topic/pn_chain
- topic/quadrupole
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 2PN ADM equality within declared hierarchy

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_2PN_ADM_WITHIN_HIERARCHY`
> **Status:** `exact_assembly_within_closure`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The conservative 2PN two-body ledger Legendre-transforms to the ADM target within the declared hierarchy; it also exposes dipole, P0/P2 support, and geometry-closure structure.

## Claim

The conservative 2PN two-body ledger Legendre-transforms to the ADM target within the declared hierarchy; it also exposes dipole, P0/P2 support, and geometry-closure structure.

## Physical Meaning

The conservative 2PN two-body ledger Legendre-transforms to the ADM target within the declared hierarchy; it also exposes dipole, P0/P2 support, and geometry-closure structure.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_assembly_within_closure`
- Outputs: `PN_2_ADM_FULL`, `PHYS_GROUPED_P2_SHAPE_RESPONSE`

## Atlas Links

### Related physical nodes
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]
- [[PHYS_INTERIOR_SUPPORT]]

### Related math nodes
- none

### Related equations
- [[EQ_2PN_ADM_EQUALITY]]

### Related claims
- [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]]
- [[CLAIM_25PN_QUAD_NARROWING]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_2PN_FULL]]
- [[SEC_2PN_CLAIMS]]
- [[SEC_2PN_MAIN_THEOREM]]
- [[SEC_2PN_OPEN_DYNAMIC]]
- [[SEC_2PN_P0P2_APPENDIX]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[PN_2_ADM_FULL]] | Claim feeds this downstream object, output, or open gate. |
| `STRUCTURAL_INPUT_TO` | [[CLAIM_25PN_QUAD_NARROWING]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_2PN_CLAIMS]] | 2PN claims and non-claims. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_2PN_MAIN_THEOREM]] | 2PN exact ADM match. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_2PN_OPEN_DYNAMIC]] | Open moving-throat dynamics. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_2PN_P0P2_APPENDIX]] | P0/P2 response structure. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_2PN_FULL]] | Paper backlink block references CLAIM_2PN_ADM_WITHIN_HIERARCHY. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_INTERIOR_SUPPORT]] | Physical ontology object grounded by this claim. |
| `LOWER_ORDER_ANCHOR_FOR` | [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]] | Claim-level dependency added in v0.4. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_2PN_FULL]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_2PN_ADM_EQUALITY]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_2PN_FULL]]
- [[SEC_2PN_CLAIMS]]
- [[SEC_2PN_MAIN_THEOREM]]
- [[SEC_2PN_OPEN_DYNAMIC]]
- [[SEC_2PN_P0P2_APPENDIX]]

### Source files
- `research/4d_2pn/paper/4d_2pn.tex`

### TeX anchor
- File: `research/4d_2pn/paper/4d_2pn.tex`
- Line: `179`
- Heading: What this paper claims
- Nearest label: `sec:intro-claims` at line `180`
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
