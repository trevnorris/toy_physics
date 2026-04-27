---
id: PN_2_ADM_FULL
title: Full conservative 2PN ADM assembly
type: PN_backbone
layer: derivation
status: full_assembly_within_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: 2PN generic-frame ADM Hamiltonian equality within declared hierarchy; reveals P0⊕P2 mouth/support layer.
future_paper_needed: false
source_files:
- research/4d_2pn/paper/4d_2pn.tex
- 4d_2pn_summary.md
legacy_sources:
- 4d_2pn_summary.md
source_links:
- '[[FILE_2PN_FULL]]'
tex_anchor:
  file: research/4d_2pn/paper/4d_2pn.tex
  line: 3380
  heading_level: subsection
  heading: Full carried conservative 1PN assembly
  nearest_label:
    name: app:carry-forward-full-1pn
    line: 3381
  nearby_labels:
  - name: app:carry-forward-full-1pn
    line: 3381
  - name: eq:app-cf-full-1pn-lagrangian
    line: 3401
  - name: eq:app-cf-exact-eih
    line: 3406
  - name: eq:app-cf-test-mass-1pn
    line: 3424
  - name: eq:app-cf-orbit-form
    line: 3433
  - name: eq:app-cf-perihelion
    line: 3443
  match_basis: semantic_heading_match
  match_score: 0.719
  confidence: medium
equation_ids:
- EQ_2PN_ADM_EQUALITY
claim_ids:
- CLAIM_2PN_ADM_WITHIN_HIERARCHY
source_ids:
- FILE_2PN_FULL
outgoing_edges:
- target: PN_2_5_QUAD_NARROWING
  relation: FEEDS
  status: derivation
  note: 2PN P0/P2 mouth/support structure feeds 2.5PN quadrupole narrowing.
- target: PN_2_5_QUAD_NARROWING
  relation: FEEDS
  status: conditional
  note: 2PN P0⊕P2 throat response sets up quadrupole narrowing.
incoming_edges:
- source: EQ_2PN_ADM_EQUALITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_2PN_FULL
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PN_1_EIH_FULL
  relation: FEEDS
  status: closure
  note: 1PN ledger carried into 2PN ADM assembly.
- source: PN_1_EIH_FULL
  relation: FEEDS
  status: derivation
  note: 2PN assembly imports 1PN closure hierarchy and EIH ledger.
- source: CLAIM_2PN_ADM_WITHIN_HIERARCHY
  relation: FEEDS_OR_STATUS_OF
  status: exact_assembly_within_closure
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/full_assembly_within_closure
- topic/pn_chain
- topic/quadrupole
- type/pn_backbone
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Full conservative 2PN ADM assembly

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PN_2_ADM_FULL`  
> **Status:** `full_assembly_within_closure`  
> **Layer:** `derivation`  
> **Type:** `PN_backbone`

## Summary

2PN generic-frame ADM Hamiltonian equality within declared hierarchy; reveals P0⊕P2 mouth/support layer.

## Physical Meaning

2PN generic-frame ADM Hamiltonian equality within declared hierarchy; reveals P0⊕P2 mouth/support layer.

## Mathematical Role

- Layer: `derivation`
- Type: `PN_backbone`
- Status: `full_assembly_within_closure`

## Equation

$$
H_2=H_2PN^ADM
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_2PN_ADM_EQUALITY]]

### Related claims
- [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_2PN_FULL]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[PN_2_5_QUAD_NARROWING]] | 2PN P0/P2 mouth/support structure feeds 2.5PN quadrupole narrowing. |
| `FEEDS` | [[PN_2_5_QUAD_NARROWING]] | 2PN P0⊕P2 throat response sets up quadrupole narrowing. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_2PN_ADM_EQUALITY]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_2PN_FULL]] | File anchor documents this node. |
| `FEEDS` | [[PN_1_EIH_FULL]] | 1PN ledger carried into 2PN ADM assembly. |
| `FEEDS` | [[PN_1_EIH_FULL]] | 2PN assembly imports 1PN closure hierarchy and EIH ledger. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_2PN_FULL]]

### Source files
- `research/4d_2pn/paper/4d_2pn.tex`
- `4d_2pn_summary.md`

### TeX anchor
- File: `research/4d_2pn/paper/4d_2pn.tex`
- Line: `3380`
- Heading: Full carried conservative 1PN assembly
- Nearest label: `app:carry-forward-full-1pn` at line `3381`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
