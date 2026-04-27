---
id: PN_2_5_QUAD_NARROWING
title: 2.5PN quadrupole-normalization narrowing
type: PN_bridge
layer: derivation
status: conditional_theorem_open_normalization
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Demotes scalar/dipole dangers and isolates orbital/worldtube STF quadrupole branch plus normalization gap.
future_paper_needed: false
source_files:
- research/4d_2_5pn/paper/4d_2_5pn.tex
- 4d_2_5pn_summary.md
legacy_sources:
- 4d_2_5pn_summary.md
source_links:
- '[[FILE_2_5PN]]'
tex_anchor:
  file: research/4d_2_5pn/paper/4d_2_5pn.tex
  line: 4995
  heading_level: section
  heading: Quadrupole normalization package
  nearest_label:
    name: app:quadrupole-normalization-package
    line: 4996
  nearby_labels:
  - name: app:quadrupole-normalization-package
    line: 4996
  match_basis: semantic_heading_match
  match_score: 0.603
  confidence: medium
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
status_firewall_ids:
- FIREWALL_25PN_CONDITIONAL
source_ids:
- FILE_2_5PN
outgoing_edges:
- target: PN_3_GROUPED_P2_FULL
  relation: FEEDS
  status: derivation
  note: 2.5PN identifies grouped real P2 as higher-order payload.
- target: PN_3_GROUPED_P2_FULL
  relation: FEEDS
  status: closure
  note: Grouped real P2 bundle becomes 3PN conservative payload.
- target: PN_4_LOCAL_TAIL
  relation: FEEDS
  status: conditional
  note: Same quadrupole normalization controls 4PN tail bridge.
- target: OPEN_QUAD_NORMALIZATION
  relation: OPENS_GATE
  status: open
  note: Surviving universal branch needs passive/outgoing normalization.
- target: OPEN_QUAD_NORMALIZATION
  relation: REQUIRES
  status: open
  note: Final 2.5PN theorem remains conditional on quadrupole normalization.
incoming_edges:
- source: FILE_2_5PN
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PN_2_ADM_FULL
  relation: FEEDS
  status: derivation
  note: 2PN P0/P2 mouth/support structure feeds 2.5PN quadrupole narrowing.
- source: PN_2_ADM_FULL
  relation: FEEDS
  status: conditional
  note: 2PN P0⊕P2 throat response sets up quadrupole narrowing.
- source: CLAIM_25PN_QUAD_NARROWING
  relation: FEEDS_OR_STATUS_OF
  status: conditional_theorem_open_normalization
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_15_25PN_4PN_INTERFACE
  relation: INTERFACES
  status: conditional
  note: Relates moving-throat quadrupole branch to 2.5PN channel.
- source: FIREWALL_25PN_CONDITIONAL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: MT_STAGE5_GROUPED_P2_BRIDGE
  relation: SUPPORTS
  status: reduced
  note: Moving-throat bridge supplies exact reduced target mhat0²P0.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/conditional_theorem_open_normalization
- topic/pn_chain
- topic/quadrupole
- type/pn_bridge
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 2.5PN quadrupole-normalization narrowing

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_2_5_QUAD_NARROWING`  
> **Status:** `conditional_theorem_open_normalization`  
> **Layer:** `derivation`  
> **Type:** `PN_bridge`

## Summary

Demotes scalar/dipole dangers and isolates orbital/worldtube STF quadrupole branch plus normalization gap.

## Physical Meaning

Demotes scalar/dipole dangers and isolates orbital/worldtube STF quadrupole branch plus normalization gap.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_bridge`
- Status: `conditional_theorem_open_normalization`

## Equation

$$
surviving branch = STF quadrupole
$$

$$
passive/outgoing normalization open
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_25PN_QUAD_NARROWING]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- [[FIREWALL_25PN_CONDITIONAL]]

### Source anchors
- [[FILE_2_5PN]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[PN_3_GROUPED_P2_FULL]] | 2.5PN identifies grouped real P2 as higher-order payload. |
| `FEEDS` | [[PN_3_GROUPED_P2_FULL]] | Grouped real P2 bundle becomes 3PN conservative payload. |
| `FEEDS` | [[PN_4_LOCAL_TAIL]] | Same quadrupole normalization controls 4PN tail bridge. |
| `OPENS_GATE` | [[OPEN_QUAD_NORMALIZATION]] | Surviving universal branch needs passive/outgoing normalization. |
| `REQUIRES` | [[OPEN_QUAD_NORMALIZATION]] | Final 2.5PN theorem remains conditional on quadrupole normalization. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_2_5PN]] | File anchor documents this node. |
| `FEEDS` | [[PN_2_ADM_FULL]] | 2PN P0/P2 mouth/support structure feeds 2.5PN quadrupole narrowing. |
| `FEEDS` | [[PN_2_ADM_FULL]] | 2PN P0⊕P2 throat response sets up quadrupole narrowing. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_25PN_QUAD_NARROWING]] | Claim feeds this downstream object, output, or open gate. |
| `INTERFACES` | [[MT_V2_15_25PN_4PN_INTERFACE]] | Relates moving-throat quadrupole branch to 2.5PN channel. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_25PN_CONDITIONAL]] | Firewall preserves this correct status boundary. |
| `SUPPORTS` | [[MT_STAGE5_GROUPED_P2_BRIDGE]] | Moving-throat bridge supplies exact reduced target mhat0²P0. |

## Source Anchors

### Source anchor notes
- [[FILE_2_5PN]]

### Source files
- `research/4d_2_5pn/paper/4d_2_5pn.tex`
- `4d_2_5pn_summary.md`

### TeX anchor
- File: `research/4d_2_5pn/paper/4d_2_5pn.tex`
- Line: `4995`
- Heading: Quadrupole normalization package
- Nearest label: `app:quadrupole-normalization-package` at line `4996`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
