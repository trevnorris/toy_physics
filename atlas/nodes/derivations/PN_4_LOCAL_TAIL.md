---
id: PN_4_LOCAL_TAIL
title: 4PN local + tail split
type: PN_bridge
layer: derivation
status: local_closed_tail_conditional
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Local instantaneous 4PN sector assembled; full 4PN theorem inherits same quadrupole normalization gap.
future_paper_needed: false
source_files:
- research/4d_4pn/paper/4d_4pn.tex
- 4d_4pn_summary.md
legacy_sources:
- 4d_4pn_summary.md
source_links:
- '[[FILE_4PN_FULL]]'
tex_anchor:
  file: research/4d_4pn/paper/4d_4pn.tex
  line: 4891
  heading_level: subsection
  heading: Exchange-symmetric local interaction scaffold
  nearest_label:
    name: eq:app-4pn-local-scaffold-split
    line: 4891
  nearby_labels:
  - name: eq:app-4pn-local-scaffold-split
    line: 4891
  match_basis: semantic_label_match
  match_score: 0.73
  confidence: high
equation_ids:
- EQ_4PN_TAIL_BRIDGE
claim_ids:
- CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
source_ids:
- FILE_4PN_FULL
outgoing_edges:
- target: PN_5_FULL_BUNDLE_SURFACE
  relation: CONTINUES_TO
  status: reduced
  note: Later full-bundle target surface and weak-axisymmetric packet.
- target: PN_5_FULL_BUNDLE_SURFACE
  relation: FEEDS
  status: derivation
  note: 4PN tail interface motivates 5PN/full-bundle continuation.
- target: OPEN_QUAD_NORMALIZATION
  relation: INHERITS_GATE
  status: open
  note: Full 4PN theorem inherits same quadrupole normalization condition.
- target: OPEN_QUAD_NORMALIZATION
  relation: REQUIRES
  status: open
  note: 4PN tail theorem inherits same quadrupole normalization gap.
incoming_edges:
- source: EQ_4PN_TAIL_BRIDGE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_4PN_FULL
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PN_2_5_QUAD_NARROWING
  relation: FEEDS
  status: conditional
  note: Same quadrupole normalization controls 4PN tail bridge.
- source: PN_3_GROUPED_P2_FULL
  relation: FEEDS
  status: derivation
  note: 3PN grouped P2 closure feeds 4PN local/tail interface.
- source: PN_3_GROUPED_P2_FULL
  relation: FEEDS
  status: closure
  note: Lower conservative ledger carried into local 4PN assembly.
- source: CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
  relation: FEEDS_OR_STATUS_OF
  status: local_closed_tail_conditional
  note: Claim feeds this downstream object, output, or open gate.
- source: MT_V2_15_25PN_4PN_INTERFACE
  relation: INTERFACES
  status: conditional
  note: Relates same quadrupole branch to 4PN tail coefficient.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/local_closed_tail_conditional
- topic/pn_chain
- topic/quadrupole
- type/pn_bridge
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 4PN local + tail split

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_4_LOCAL_TAIL`
> **Status:** `local_closed_tail_conditional`
> **Layer:** `derivation`
> **Type:** `PN_bridge`

## Summary

Local instantaneous 4PN sector assembled; full 4PN theorem inherits same quadrupole normalization gap.

## Physical Meaning

Local instantaneous 4PN sector assembled; full 4PN theorem inherits same quadrupole normalization gap.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_bridge`
- Status: `local_closed_tail_conditional`

## Equation

$$
L_4PN=L_local+L_tail
$$

$$
C_tail=(GM/2c^3)gamma_quad_eff
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_4PN_TAIL_BRIDGE]]

### Related claims
- [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- none

### Source anchors
- [[FILE_4PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `CONTINUES_TO` | [[PN_5_FULL_BUNDLE_SURFACE]] | Later full-bundle target surface and weak-axisymmetric packet. |
| `FEEDS` | [[PN_5_FULL_BUNDLE_SURFACE]] | 4PN tail interface motivates 5PN/full-bundle continuation. |
| `INHERITS_GATE` | [[OPEN_QUAD_NORMALIZATION]] | Full 4PN theorem inherits same quadrupole normalization condition. |
| `REQUIRES` | [[OPEN_QUAD_NORMALIZATION]] | 4PN tail theorem inherits same quadrupole normalization gap. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_4PN_TAIL_BRIDGE]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_4PN_FULL]] | File anchor documents this node. |
| `FEEDS` | [[PN_2_5_QUAD_NARROWING]] | Same quadrupole normalization controls 4PN tail bridge. |
| `FEEDS` | [[PN_3_GROUPED_P2_FULL]] | 3PN grouped P2 closure feeds 4PN local/tail interface. |
| `FEEDS` | [[PN_3_GROUPED_P2_FULL]] | Lower conservative ledger carried into local 4PN assembly. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]] | Claim feeds this downstream object, output, or open gate. |
| `INTERFACES` | [[MT_V2_15_25PN_4PN_INTERFACE]] | Relates same quadrupole branch to 4PN tail coefficient. |

## Source Anchors

### Source anchor notes
- [[FILE_4PN_FULL]]

### Source files
- `research/4d_4pn/paper/4d_4pn.tex`
- `4d_4pn_summary.md`

### TeX anchor
- File: `research/4d_4pn/paper/4d_4pn.tex`
- Line: `4891`
- Heading: Exchange-symmetric local interaction scaffold
- Nearest label: `eq:app-4pn-local-scaffold-split` at line `4891`
- Match basis: `semantic_label_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
