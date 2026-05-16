---
id: EQ_G2_COMMON_TANGENT
title: g-2 common tangent target
type: equation
layer: equation_anchor
status: conditional_open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: First tangent projection required by the minimal common quotient correction to match the residual.
future_paper_needed: false
source_files:
- notes/g2/g2_full_output.md
- g2_full_output.md
legacy_sources:
- g2_full_output.md
source_links:
- '[[FILE_G2_OUTPUT]]'
equation_ids:
- EQ_MONOMIAL_QUOTIENT
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
source_ids:
- FILE_G2_OUTPUT
outgoing_edges:
- target: ANOMALY_G2_COMMON_QUOTIENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: SUPPORTS_CLAIM
  status: conditional_reduced_residual
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_G2_OUTPUT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_G2_COMMON_TANGENT.
- source: FILE_G2_OUTPUT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_G2_COMMON_TANGENT.
- source: EQ_MONOMIAL_QUOTIENT
  relation: FEEDS
  status: conditional
  note: g-2 common tangent target is expressed in quotient coordinates.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/conditional_open
- topic/g2
- topic/projection
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# g-2 common tangent target

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_G2_COMMON_TANGENT`
> **Status:** `conditional_open`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

First tangent projection required by the minimal common quotient correction to match the residual.

## Physical Meaning

First tangent projection required by the minimal common quotient correction to match the residual.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `conditional_open`
- Parent node: `ANOMALY_G2_COMMON_QUOTIENT`

## Equation

$$
Λ_1 ≈ 0.279605891931464
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_MONOMIAL_QUOTIENT]]

### Related claims
- [[CLAIM_G2_COMMON_QUOTIENT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_G2_OUTPUT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[ANOMALY_G2_COMMON_QUOTIENT]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_G2_COMMON_QUOTIENT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_G2_OUTPUT]] | Paper backlink block references EQ_G2_COMMON_TANGENT. |
| `CONTAINS_EQUATION` | [[FILE_G2_OUTPUT]] | Source artifact contains or supports EQ_G2_COMMON_TANGENT. |
| `FEEDS` | [[EQ_MONOMIAL_QUOTIENT]] | g-2 common tangent target is expressed in quotient coordinates. |

## Source Anchors

### Source anchor notes
- [[FILE_G2_OUTPUT]]

### Source files
- `notes/g2/g2_full_output.md`
- `g2_full_output.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
