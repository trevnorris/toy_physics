---
id: EQ_XI1_PREF_SLOPE
title: Weak-axisymmetric prefactor slope
type: equation
layer: equation_anchor
status: exact_within_reduced_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The weak-axisymmetric defect also measures outgoing-prefactor / transfer-shape slope.
future_paper_needed: false
source_files:
- notes/5pn/5pn_notes_full.md
- notes/g2/g2_full_output.md
- notes/pde_audit_full.md
- 5pn_notes_full.md
- g2_full_output.md
- pde_audit_full.md
legacy_sources:
- 5pn_notes_full.md
- g2_full_output.md
- pde_audit_full.md
source_links:
- '[[FILE_5PN_FULL]]'
- '[[FILE_G2_OUTPUT]]'
- '[[FILE_PDE_AUDIT]]'
math_ids:
- MATH_XI1_PREF_SLOPE
equation_ids:
- EQ_MONOMIAL_QUOTIENT
- EQ_WEAK_AXISYM_SIGNATURE
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
- CLAIM_PACKET_A_PACKET_B_SPLIT
source_ids:
- FILE_5PN_FULL
- FILE_G2_OUTPUT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_XI1_PREF_SLOPE
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: SUPPORTS_CLAIM
  status: conditional_reduced_residual
  note: Equation anchor supports this named claim.
- target: CLAIM_PACKET_A_PACKET_B_SPLIT
  relation: SUPPORTS_CLAIM
  status: open_branch_packets
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_5PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_XI1_PREF_SLOPE.
- source: BACKLINK_G2_OUTPUT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_XI1_PREF_SLOPE.
- source: FILE_5PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_XI1_PREF_SLOPE.
- source: FILE_G2_OUTPUT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_XI1_PREF_SLOPE.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_XI1_PREF_SLOPE.
- source: EQ_WEAK_AXISYM_SIGNATURE
  relation: ORGANIZES
  status: exact_first_order
  note: Weak-axisymmetric signature carries prefactor slope Xi1.
- source: EQ_MONOMIAL_QUOTIENT
  relation: PARAMETERIZES
  status: exact_within_branch
  note: Quotient coordinates parameterize the Xi1 branch defect.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_reduced_branch
- topic/moving_throat
- topic/pn_chain
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weak-axisymmetric prefactor slope

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_XI1_PREF_SLOPE`
> **Status:** `exact_within_reduced_branch`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

The weak-axisymmetric defect also measures outgoing-prefactor / transfer-shape slope.

## Physical Meaning

The weak-axisymmetric defect also measures outgoing-prefactor / transfer-shape slope.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_reduced_branch`
- Parent node: [[MATH_XI1_PREF_SLOPE]]

## Equation

$$
Ξ₁ = P₁/P₀ = δln(T_eff²)/(ελ_A)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_XI1_PREF_SLOPE]]

### Related equations
- [[EQ_MONOMIAL_QUOTIENT]]
- [[EQ_WEAK_AXISYM_SIGNATURE]]

### Related claims
- [[CLAIM_G2_COMMON_QUOTIENT]]
- [[CLAIM_PACKET_A_PACKET_B_SPLIT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_5PN_FULL]]
- [[FILE_G2_OUTPUT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_XI1_PREF_SLOPE]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_G2_COMMON_QUOTIENT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_PACKET_A_PACKET_B_SPLIT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_5PN_FULL]] | Paper backlink block references EQ_XI1_PREF_SLOPE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_G2_OUTPUT]] | Paper backlink block references EQ_XI1_PREF_SLOPE. |
| `CONTAINS_EQUATION` | [[FILE_5PN_FULL]] | Source artifact contains or supports EQ_XI1_PREF_SLOPE. |
| `CONTAINS_EQUATION` | [[FILE_G2_OUTPUT]] | Source artifact contains or supports EQ_XI1_PREF_SLOPE. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_XI1_PREF_SLOPE. |
| `ORGANIZES` | [[EQ_WEAK_AXISYM_SIGNATURE]] | Weak-axisymmetric signature carries prefactor slope Xi1. |
| `PARAMETERIZES` | [[EQ_MONOMIAL_QUOTIENT]] | Quotient coordinates parameterize the Xi1 branch defect. |

## Source Anchors

### Source anchor notes
- [[FILE_5PN_FULL]]
- [[FILE_G2_OUTPUT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `notes/5pn/5pn_notes_full.md`
- `notes/g2/g2_full_output.md`
- `notes/pde_audit_full.md`
- `5pn_notes_full.md`
- `g2_full_output.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
