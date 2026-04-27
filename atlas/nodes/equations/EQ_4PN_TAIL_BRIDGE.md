---
id: EQ_4PN_TAIL_BRIDGE
title: 4PN tail coefficient bridge
type: equation
layer: equation_anchor
status: conditional_on_quadrupole_normalization
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: 4PN hereditary coefficient is not new; it is fixed by the same quadrupole normalization as 2.5PN.
future_paper_needed: false
source_files:
- research/4d_4pn/paper/4d_4pn.tex
- 4d_4pn_summary.md
- 4d_4pn_full_notes.md
legacy_sources:
- 4d_4pn_summary.md
- 4d_4pn_full_notes.md
source_links:
- '[[FILE_4PN_FULL]]'
tex_anchor:
  file: research/4d_4pn/paper/4d_4pn.tex
  line: 7473
  heading_level: subsection
  heading: Exact bridge to the carried 2.5PN quadrupole coefficient
  nearest_label:
    name: eq:app-tail-exact-coefficient-bridge
    line: 7473
  nearby_labels:
  - name: eq:app-tail-exact-coefficient-bridge
    line: 7473
  match_basis: semantic_label_match
  match_score: 0.77
  confidence: high
equation_ids:
- EQ_P0_TARGET
claim_ids:
- CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
status_firewall_ids:
- FIREWALL_4PN_LOCAL_NOT_FULL_TAIL
source_ids:
- FILE_4PN_FULL
outgoing_edges:
- target: PN_4_LOCAL_TAIL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
  relation: SUPPORTS_CLAIM
  status: local_closed_tail_conditional
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_4PN_TAIL_BRIDGE.
- source: FILE_4PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_4PN_TAIL_BRIDGE.
- source: FIREWALL_4PN_LOCAL_NOT_FULL_TAIL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: EQ_P0_TARGET
  relation: REQUIRED_BY
  status: conditional
  note: The 4PN tail bridge uses the same quadrupole normalization.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/conditional_on_quadrupole_normalization
- topic/pn_chain
- topic/quadrupole
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 4PN tail coefficient bridge

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_4PN_TAIL_BRIDGE`  
> **Status:** `conditional_on_quadrupole_normalization`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

4PN hereditary coefficient is not new; it is fixed by the same quadrupole normalization as 2.5PN.

## Physical Meaning

4PN hereditary coefficient is not new; it is fixed by the same quadrupole normalization as 2.5PN.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `conditional_on_quadrupole_normalization`
- Parent node: `PN_4_LOCAL_TAIL`

## Equation

$$
C_tail = (GM/(2c³)) γ_quad^eff
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_P0_TARGET]]

### Related claims
- [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_4PN_LOCAL_NOT_FULL_TAIL]]

### Source anchors
- [[FILE_4PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[PN_4_LOCAL_TAIL]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4PN_FULL]] | Paper backlink block references EQ_4PN_TAIL_BRIDGE. |
| `CONTAINS_EQUATION` | [[FILE_4PN_FULL]] | Source artifact contains or supports EQ_4PN_TAIL_BRIDGE. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_4PN_LOCAL_NOT_FULL_TAIL]] | Firewall preserves this correct status boundary. |
| `REQUIRED_BY` | [[EQ_P0_TARGET]] | The 4PN tail bridge uses the same quadrupole normalization. |

## Source Anchors

### Source anchor notes
- [[FILE_4PN_FULL]]

### Source files
- `research/4d_4pn/paper/4d_4pn.tex`
- `4d_4pn_summary.md`
- `4d_4pn_full_notes.md`

### TeX anchor
- File: `research/4d_4pn/paper/4d_4pn.tex`
- Line: `7473`
- Heading: Exact bridge to the carried 2.5PN quadrupole coefficient
- Nearest label: `eq:app-tail-exact-coefficient-bridge` at line `7473`
- Match basis: `semantic_label_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
