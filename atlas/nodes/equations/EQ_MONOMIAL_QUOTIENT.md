---
id: EQ_MONOMIAL_QUOTIENT
title: Monomial quotient drift map
type: equation
layer: equation_anchor
status: exact_within_coherent_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Separates true branch motion from five-dimensional similarity-orbit motion.
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
- MATH_MONOMIAL_QUOTIENT
equation_ids:
- EQ_G2_COMMON_TANGENT
- EQ_XI1_PREF_SLOPE
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_G2_COMMON_QUOTIENT
source_ids:
- FILE_5PN_FULL
- FILE_G2_OUTPUT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_MONOMIAL_QUOTIENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_G2_COMMON_TANGENT
  relation: FEEDS
  status: conditional
  note: g-2 common tangent target is expressed in quotient coordinates.
- target: EQ_XI1_PREF_SLOPE
  relation: PARAMETERIZES
  status: exact_within_branch
  note: Quotient coordinates parameterize the Xi1 branch defect.
- target: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Equation anchor supports this named claim.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: SUPPORTS_CLAIM
  status: conditional_reduced_residual
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_5PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MONOMIAL_QUOTIENT.
- source: FILE_G2_OUTPUT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MONOMIAL_QUOTIENT.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MONOMIAL_QUOTIENT.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_coherent_branch
- topic/moving_throat
- topic/pn_chain
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Monomial quotient drift map

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_MONOMIAL_QUOTIENT`
> **Status:** `exact_within_coherent_branch`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Separates true branch motion from five-dimensional similarity-orbit motion.

## Physical Meaning

Separates true branch motion from five-dimensional similarity-orbit motion.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_coherent_branch`
- Parent node: [[MATH_MONOMIAL_QUOTIENT]]

## Equation

$$
(q_tr,q_nt,q_eta)^T = M_* δx,    rank(M_*)=3
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MONOMIAL_QUOTIENT]]

### Related equations
- [[EQ_G2_COMMON_TANGENT]]
- [[EQ_XI1_PREF_SLOPE]]

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_G2_COMMON_QUOTIENT]]

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
| `ANCHORS` | [[MATH_MONOMIAL_QUOTIENT]] | Equation anchor belongs to or formalizes this graph node. |
| `FEEDS` | [[EQ_G2_COMMON_TANGENT]] | g-2 common tangent target is expressed in quotient coordinates. |
| `PARAMETERIZES` | [[EQ_XI1_PREF_SLOPE]] | Quotient coordinates parameterize the Xi1 branch defect. |
| `SUPPORTS_CLAIM` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_G2_COMMON_QUOTIENT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_5PN_FULL]] | Source artifact contains or supports EQ_MONOMIAL_QUOTIENT. |
| `CONTAINS_EQUATION` | [[FILE_G2_OUTPUT]] | Source artifact contains or supports EQ_MONOMIAL_QUOTIENT. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_MONOMIAL_QUOTIENT. |

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
