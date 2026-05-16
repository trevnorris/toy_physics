---
id: CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL
title: 4PN local sector closed; tail inherits quadrupole gate
type: claim
layer: claim_theorem
status: local_closed_tail_conditional
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The local instantaneous 4PN sector is assembled in the declared hierarchy; the full 4PN theorem depends on the same passive/outgoing STF quadrupole normalization isolated at 2.5PN.
future_paper_needed: false
source_files:
- research/4d_4pn/paper/4d_4pn.tex
source_links:
- '[[FILE_4PN_FULL]]'
- '[[SEC_4PN_FINAL_THEOREM]]'
- '[[SEC_4PN_NO_NEW_GAP]]'
- '[[SEC_4PN_TAIL_BRIDGE]]'
tex_anchor:
  file: research/4d_4pn/paper/4d_4pn.tex
  line: 2161
  heading_level: section
  heading: Tail / hereditary 4PN bridge
  nearest_label:
    name: sec:tail-bridge
    line: 2162
  nearby_labels:
  - name: sec:tail-bridge
    line: 2162
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 1.0
  confidence: medium
  source_anchor_node: SEC_4PN_TAIL_BRIDGE
physical_ids:
- PHYS_OUTGOING_QUADRUPOLE_PORT
equation_ids:
- EQ_4PN_TAIL_BRIDGE
- EQ_P0_TARGET
claim_ids:
- CLAIM_25PN_QUAD_NARROWING
open_gate_ids:
- OPEN_QUAD_NORMALIZATION
status_firewall_ids:
- FIREWALL_4PN_LOCAL_NOT_FULL_TAIL
source_ids:
- FILE_4PN_FULL
- SEC_4PN_FINAL_THEOREM
- SEC_4PN_NO_NEW_GAP
- SEC_4PN_TAIL_BRIDGE
outgoing_edges:
- target: OPEN_QUAD_NORMALIZATION
  relation: FEEDS_OR_STATUS_OF
  status: local_closed_tail_conditional
  note: Claim feeds this downstream object, output, or open gate.
- target: PN_4_LOCAL_TAIL
  relation: FEEDS_OR_STATUS_OF
  status: local_closed_tail_conditional
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_4PN_FINAL_THEOREM
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Final 4PN theorem split.
- source: SEC_4PN_NO_NEW_GAP
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: 4PN inherits 2.5PN normalization.
- source: SEC_4PN_TAIL_BRIDGE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: 4PN tail interface.
- source: BACKLINK_4PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: local_closed_tail_conditional
  note: 'Claim class: exact_within_closure_plus_open_gate'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_OUTGOING_QUADRUPOLE_PORT
  relation: GROUNDS_PHYSICAL_MEANING
  status: local_closed_tail_conditional
  note: Physical ontology object grounded by this claim.
- source: FILE_4PN_FULL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: local_closed_tail_conditional
  note: Source artifact anchors this claim.
- source: FIREWALL_4PN_LOCAL_NOT_FULL_TAIL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_LOCAL_4PN_CLOSES_FULL_TAIL
  relation: STARTS_AT
  status: v07
  note: Negative query starts from CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL.
- source: QUERY_4PN_TAIL_BLOCKER
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
- source: CLAIM_25PN_QUAD_NARROWING
  relation: SUPPLIES_TAIL_GATE_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_4PN_TAIL_BRIDGE
  relation: SUPPORTS_CLAIM
  status: local_closed_tail_conditional
  note: Equation anchor supports this named claim.
- source: EQ_P0_TARGET
  relation: SUPPORTS_CLAIM
  status: local_closed_tail_conditional
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/local_closed_tail_conditional
- topic/pn_chain
- topic/quadrupole
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 4PN local sector closed; tail inherits quadrupole gate

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL`
> **Status:** `local_closed_tail_conditional`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The local instantaneous 4PN sector is assembled in the declared hierarchy; the full 4PN theorem depends on the same passive/outgoing STF quadrupole normalization isolated at 2.5PN.

## Claim

The local instantaneous 4PN sector is assembled in the declared hierarchy; the full 4PN theorem depends on the same passive/outgoing STF quadrupole normalization isolated at 2.5PN.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

The local instantaneous 4PN sector is assembled in the declared hierarchy; the full 4PN theorem depends on the same passive/outgoing STF quadrupole normalization isolated at 2.5PN.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `local_closed_tail_conditional`
- Outputs: `PN_4_LOCAL_TAIL`, `OPEN_QUAD_NORMALIZATION`

## Atlas Links

### Related physical nodes
- [[PHYS_OUTGOING_QUADRUPOLE_PORT]]

### Related math nodes
- none

### Related equations
- [[EQ_4PN_TAIL_BRIDGE]]
- [[EQ_P0_TARGET]]

### Related claims
- [[CLAIM_25PN_QUAD_NARROWING]]

### Open gates
- [[OPEN_QUAD_NORMALIZATION]]

### Status firewalls
- [[FIREWALL_4PN_LOCAL_NOT_FULL_TAIL]]

### Source anchors
- [[FILE_4PN_FULL]]
- [[SEC_4PN_FINAL_THEOREM]]
- [[SEC_4PN_NO_NEW_GAP]]
- [[SEC_4PN_TAIL_BRIDGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[OPEN_QUAD_NORMALIZATION]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[PN_4_LOCAL_TAIL]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_4PN_FINAL_THEOREM]] | Final 4PN theorem split. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4PN_NO_NEW_GAP]] | 4PN inherits 2.5PN normalization. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4PN_TAIL_BRIDGE]] | 4PN tail interface. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4PN_FULL]] | Paper backlink block references CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure_plus_open_gate |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_OUTGOING_QUADRUPOLE_PORT]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_4PN_FULL]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_4PN_LOCAL_NOT_FULL_TAIL]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_LOCAL_4PN_CLOSES_FULL_TAIL]] | Negative query starts from CLAIM_4PN_LOCAL_CLOSED_TAIL_CONDITIONAL. |
| `STARTS_AT` | [[QUERY_4PN_TAIL_BLOCKER]] | Query validation start node. |
| `SUPPLIES_TAIL_GATE_FOR` | [[CLAIM_25PN_QUAD_NARROWING]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_4PN_TAIL_BRIDGE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_P0_TARGET]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_4PN_FULL]]
- [[SEC_4PN_FINAL_THEOREM]]
- [[SEC_4PN_NO_NEW_GAP]]
- [[SEC_4PN_TAIL_BRIDGE]]

### Source files
- `research/4d_4pn/paper/4d_4pn.tex`

### TeX anchor
- File: `research/4d_4pn/paper/4d_4pn.tex`
- Line: `2161`
- Heading: Tail / hereditary 4PN bridge
- Nearest label: `sec:tail-bridge` at line `2162`
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
