---
id: EQ_WEAK_AXISYM_SIGNATURE
title: Weak-axisymmetric grouped signature
type: equation
layer: equation_anchor
status: exact_angular_first_order
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Diagnostic first-order splitting law for a weak axisymmetric quadrupole perturbation.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_007.tex
- notes/pde_audit_full.md
- notes/5pn/5pn_notes_full.md
- moving_throat_pde_stage007_overlap_isotropy.md
- pde_audit_full.md
- 5pn_notes_full.md
legacy_sources:
- moving_throat_pde_stage007_overlap_isotropy.md
- pde_audit_full.md
- 5pn_notes_full.md
source_links:
- '[[FILE_5PN_FULL]]'
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
tex_anchor:
  file: research/pde_ledger/paper/stages/stage_007.tex
  line: 149
  heading_level: paragraph
  heading: 4. Weak-axisymmetric quadrupole splitting.
  nearest_label:
    name: eq:app-stage007-deltaK-axisym
    line: 152
  nearby_labels:
  - name: eq:app-stage007-deltaK-axisym
    line: 152
  - name: eq:app-stage007-MAB20
    line: 157
  - name: eq:app-stage007-axisym-matrix
    line: 162
  - name: eq:app-stage007-lambda-signature
    line: 169
  - name: eq:app-stage007-x-signature
    line: 179
  - name: eq:app-stage007-b-3a
    line: 188
  match_basis: semantic_heading_match
  match_score: 0.612
  confidence: medium
math_ids:
- MATH_WEAK_AXISYM_SPLITTING
equation_ids:
- EQ_XI1_PREF_SLOPE
claim_ids:
- CLAIM_5PN_FULL_BUNDLE_SURFACE
- CLAIM_STAGE7_O3_ISOTROPY
source_ids:
- FILE_5PN_FULL
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_WEAK_AXISYM_SPLITTING
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_XI1_PREF_SLOPE
  relation: ORGANIZES
  status: exact_first_order
  note: Weak-axisymmetric signature carries prefactor slope Xi1.
- target: CLAIM_5PN_FULL_BUNDLE_SURFACE
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Equation anchor supports this named claim.
- target: CLAIM_STAGE7_O3_ISOTROPY
  relation: SUPPORTS_CLAIM
  status: exact_angular_reduced
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_5PN_FULL
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE.
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_angular_first_order
- topic/moving_throat
- topic/pn_chain
- topic/quadrupole
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Weak-axisymmetric grouped signature

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_WEAK_AXISYM_SIGNATURE`  
> **Status:** `exact_angular_first_order`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Diagnostic first-order splitting law for a weak axisymmetric quadrupole perturbation.

## Physical Meaning

Diagnostic first-order splitting law for a weak axisymmetric quadrupole perturbation.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_angular_first_order`
- Parent node: [[MATH_WEAK_AXISYM_SPLITTING]]

## Equation

$$
(20,21,22) ~ (1, 1/2, -1), equivalently b=3a
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_WEAK_AXISYM_SPLITTING]]

### Related equations
- [[EQ_XI1_PREF_SLOPE]]

### Related claims
- [[CLAIM_5PN_FULL_BUNDLE_SURFACE]]
- [[CLAIM_STAGE7_O3_ISOTROPY]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_5PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_WEAK_AXISYM_SPLITTING]] | Equation anchor belongs to or formalizes this graph node. |
| `ORGANIZES` | [[EQ_XI1_PREF_SLOPE]] | Weak-axisymmetric signature carries prefactor slope Xi1. |
| `SUPPORTS_CLAIM` | [[CLAIM_5PN_FULL_BUNDLE_SURFACE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_STAGE7_O3_ISOTROPY]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_5PN_FULL]] | Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE. |
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_WEAK_AXISYM_SIGNATURE. |

## Source Anchors

### Source anchor notes
- [[FILE_5PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/pde_ledger/paper/stages/stage_007.tex`
- `notes/pde_audit_full.md`
- `notes/5pn/5pn_notes_full.md`
- `moving_throat_pde_stage007_overlap_isotropy.md`
- `pde_audit_full.md`
- `5pn_notes_full.md`

### TeX anchor
- File: `research/pde_ledger/paper/stages/stage_007.tex`
- Line: `149`
- Heading: 4. Weak-axisymmetric quadrupole splitting.
- Nearest label: `eq:app-stage007-deltaK-axisym` at line `152`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
