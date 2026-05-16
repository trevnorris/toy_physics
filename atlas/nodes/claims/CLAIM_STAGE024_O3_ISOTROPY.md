---
id: CLAIM_STAGE024_O3_ISOTROPY
title: Stage 024 O(3) isotropy and weak-axisymmetric splitting theorem
type: claim
layer: claim_theorem
status: exact_angular_reduced
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The angular source map is identity on the canonical real STF basis; O(3)-invariant kernels force grouped-lane isotropy, while weak axisymmetric splitting has signature (1,1/2,-1...
future_paper_needed: false
source_files:
- research/4d_2_5pn/paper/4d_2_5pn.tex
source_links:
- '[[FILE_5PN_FULL]]'
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[SEC_25PN_STF_MAP]]'
- '[[SEC_PDE_P2_PROJECTORS]]'
- '[[SEC_PDE_STF_SOURCE]]'
tex_anchor:
  file: research/4d_2_5pn/paper/4d_2_5pn.tex
  line: 2248
  heading_level: subsection
  heading: The universal source is the orbital/worldtube STF quadrupole
  nearest_label:
    name: sec:quadrupole-surviving-source-map
    line: 2249
  nearby_labels:
  - name: sec:quadrupole-surviving-source-map
    line: 2249
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 0.583
  confidence: medium
  source_anchor_node: SEC_25PN_STF_MAP
physical_ids:
- PHYS_FINITE_MOUTH_SHAPE
- PHYS_GROUPED_P2_SHAPE_RESPONSE
- PHYS_REG_P2_QUAD
math_ids:
- MATH_STF_SOURCE_MAP
- MATH_WEAK_AXISYM_SPLITTING
equation_ids:
- EQ_WEAK_AXISYM_SIGNATURE
claim_ids:
- CLAIM_STAGE023_FULL_BUNDLE_RATIO
- CLAIM_STAGE025_031_SELECTED_BRANCH
open_gate_ids:
- OPEN_WEAK_AXISYM_ORBIT_LOCK
source_ids:
- FILE_5PN_FULL
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- SEC_25PN_STF_MAP
- SEC_PDE_P2_PROJECTORS
- SEC_PDE_STF_SOURCE
outgoing_edges:
- target: MATH_STF_SOURCE_MAP
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_WEAK_AXISYM_SPLITTING
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_STAGE024_OVERLAP_ISOTROPY
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- target: OPEN_WEAK_AXISYM_ORBIT_LOCK
  relation: FEEDS_OR_STATUS_OF
  status: exact_angular_reduced
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE025_031_SELECTED_BRANCH
  relation: RADIAL_AXIAL_CONTINUATION_FOR
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_25PN_STF_MAP
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: STF/source map.
- source: SEC_PDE_P2_PROJECTORS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Grouped P2 projectors.
- source: SEC_PDE_STF_SOURCE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: STF angular source-map theorem.
- source: CLAIM_STAGE023_FULL_BUNDLE_RATIO
  relation: ANGULARLY_REFINED_BY
  status: active
  note: Claim-level dependency added in v0.4.
- source: BACKLINK_25PN
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE024_O3_ISOTROPY.
- source: BACKLINK_PDE_AUDIT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_STAGE024_O3_ISOTROPY.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_angular_reduced
  note: 'Claim class: exact_within_closure'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: QUERY_P2_BECOMES_PHYSICAL
  relation: EXPECTS_TARGET
  status: v06
  note: Query validation expected target node.
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_angular_reduced
  note: Physical ontology object grounded by this claim.
- source: PHYS_GROUPED_P2_SHAPE_RESPONSE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_angular_reduced
  note: Physical ontology object grounded by this claim.
- source: PHYS_REG_P2_QUAD
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FILE_5PN_FULL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_angular_reduced
  note: Source artifact anchors this claim.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_angular_reduced
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_angular_reduced
  note: Source artifact anchors this claim.
- source: EQ_WEAK_AXISYM_SIGNATURE
  relation: SUPPORTS_CLAIM
  status: exact_angular_reduced
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_angular_reduced
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 024 O(3) isotropy and weak-axisymmetric splitting theorem

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_STAGE024_O3_ISOTROPY`
> **Status:** `exact_angular_reduced`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The angular source map is identity on the canonical real STF basis; O(3)-invariant kernels force grouped-lane isotropy, while weak axisymmetric splitting has signature (1,1/2,-1) or b=3a.

## Claim

The angular source map is identity on the canonical real STF basis; O(3)-invariant kernels force grouped-lane isotropy, while weak axisymmetric splitting has signature (1,1/2,-1) or b=3a.

## Physical Meaning

The angular source map is identity on the canonical real STF basis; O(3)-invariant kernels force grouped-lane isotropy, while weak axisymmetric splitting has signature (1,1/2,-1) or b=3a.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_angular_reduced`
- Outputs: `MT_STAGE024_OVERLAP_ISOTROPY`, `MATH_STF_SOURCE_MAP`, `MATH_WEAK_AXISYM_SPLITTING`, `OPEN_WEAK_AXISYM_ORBIT_LOCK`

## Atlas Links

### Related physical nodes
- [[PHYS_FINITE_MOUTH_SHAPE]]
- [[PHYS_GROUPED_P2_SHAPE_RESPONSE]]
- [[PHYS_REG_P2_QUAD]]

### Related math nodes
- [[MATH_STF_SOURCE_MAP]]
- [[MATH_WEAK_AXISYM_SPLITTING]]

### Related equations
- [[EQ_WEAK_AXISYM_SIGNATURE]]

### Related claims
- [[CLAIM_STAGE023_FULL_BUNDLE_RATIO]]
- [[CLAIM_STAGE025_031_SELECTED_BRANCH]]

### Open gates
- [[OPEN_WEAK_AXISYM_ORBIT_LOCK]]

### Status firewalls
- none

### Source anchors
- [[FILE_5PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_25PN_STF_MAP]]
- [[SEC_PDE_P2_PROJECTORS]]
- [[SEC_PDE_STF_SOURCE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[MATH_STF_SOURCE_MAP]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_WEAK_AXISYM_SPLITTING]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_STAGE024_OVERLAP_ISOTROPY]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_WEAK_AXISYM_ORBIT_LOCK]] | Claim feeds this downstream object, output, or open gate. |
| `RADIAL_AXIAL_CONTINUATION_FOR` | [[CLAIM_STAGE025_031_SELECTED_BRANCH]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_25PN_STF_MAP]] | STF/source map. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_P2_PROJECTORS]] | Grouped P2 projectors. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_STF_SOURCE]] | STF angular source-map theorem. |
| `ANGULARLY_REFINED_BY` | [[CLAIM_STAGE023_FULL_BUNDLE_RATIO]] | Claim-level dependency added in v0.4. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_25PN]] | Paper backlink block references CLAIM_STAGE024_O3_ISOTROPY. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block references CLAIM_STAGE024_O3_ISOTROPY. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `EXPECTS_TARGET` | [[QUERY_P2_BECOMES_PHYSICAL]] | Query validation expected target node. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_FINITE_MOUTH_SHAPE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_GROUPED_P2_SHAPE_RESPONSE]] | Physical ontology object grounded by this claim. |
| `LINKS_TO` | [[PHYS_REG_P2_QUAD]] | Physical register entry links to graph object. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_5PN_FULL]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `SUPPORTS_CLAIM` | [[EQ_WEAK_AXISYM_SIGNATURE]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_5PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[SEC_25PN_STF_MAP]]
- [[SEC_PDE_P2_PROJECTORS]]
- [[SEC_PDE_STF_SOURCE]]

### Source files
- `research/4d_2_5pn/paper/4d_2_5pn.tex`

### TeX anchor
- File: `research/4d_2_5pn/paper/4d_2_5pn.tex`
- Line: `2248`
- Heading: The universal source is the orbital/worldtube STF quadrupole
- Nearest label: `sec:quadrupole-surviving-source-map` at line `2249`
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
