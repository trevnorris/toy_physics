---
id: EQ_PARENT_ACTION_CURRENT
title: Current parent action
type: equation
layer: equation_anchor
status: exact_declared_parent_with_geometry_argument
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Current declared parent has GNLS + localized Maxwell, with geometry through V_conf(X;Sigma).
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/pde_audit_full.md
- 4d_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- pde_audit_full.md
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_PDE_AUDIT]]'
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 1860
  heading_level: paragraph
  heading: Matter current from the action.
  nearest_label:
    name: eq:Jpsi-variational
    line: 1867
  nearby_labels:
  - name: eq:Jpsi-variational
    line: 1867
  - name: eq:Jpsi-spatial
    line: 1875
  - name: eq:Jpsi-time
    line: 1880
  match_basis: semantic_heading_match
  match_score: 0.667
  confidence: medium
math_ids:
- MATH_PARENT_ACTION_CURRENT
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
- CLAIM_PARENT_WALL_STATUS_SPLIT
source_ids:
- FILE_4D_PARENT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_PARENT_ACTION_CURRENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPORTS_CLAIM
  status: exact_parent_current
  note: Equation anchor supports this named claim.
- target: CLAIM_PARENT_WALL_STATUS_SPLIT
  relation: SUPPORTS_CLAIM
  status: strict_parent_fail_effective_wall_pass
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_PARENT_ACTION_CURRENT.
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_PARENT_ACTION_CURRENT.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_PARENT_ACTION_CURRENT.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_declared_parent_with_geometry_argument
- topic/maxwell
- topic/moving_throat
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Current parent action

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_PARENT_ACTION_CURRENT`  
> **Status:** `exact_declared_parent_with_geometry_argument`  
> **Layer:** `equation_anchor`  
> **Type:** `equation`

## Summary

Current declared parent has GNLS + localized Maxwell, with geometry through V_conf(X;Sigma).

## Physical Meaning

Current declared parent has GNLS + localized Maxwell, with geometry through V_conf(X;Sigma).

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_declared_parent_with_geometry_argument`
- Parent node: [[MATH_PARENT_ACTION_CURRENT]]

## Equation

$$
S_current = ∫ dt d^4X [L_psi(psi,A;Sigma)+L_EM(A)]
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_PARENT_ACTION_CURRENT]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]
- [[CLAIM_PARENT_WALL_STATUS_SPLIT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_PARENT_ACTION_CURRENT]] | Equation anchor belongs to or formalizes this graph node. |
| `SUPPORTS_CLAIM` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_PARENT_WALL_STATUS_SPLIT]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references EQ_PARENT_ACTION_CURRENT. |
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_PARENT_ACTION_CURRENT. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_PARENT_ACTION_CURRENT. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/4d/paper/4d.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `1860`
- Heading: Matter current from the action.
- Nearest label: `eq:Jpsi-variational` at line `1867`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
