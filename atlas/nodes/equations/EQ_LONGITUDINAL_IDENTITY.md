---
id: EQ_LONGITUDINAL_IDENTITY
title: Exact longitudinal identity
type: equation
layer: equation_anchor
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Exact projected identity that contains the Poisson-hook Laplacian.
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
  line: 907
  heading_level: subsection
  heading: Exact longitudinal identity from projected continuity
  nearest_label:
    name: sec:poisson-exact
    line: 908
  nearby_labels:
  - name: sec:poisson-exact
    line: 908
  - name: eq:proj-cont-again
    line: 913
  - name: eq:longitudinal-identity
    line: 946
  match_basis: semantic_heading_match
  match_score: 0.962
  confidence: high
math_ids:
- MATH_LONGITUDINAL_IDENTITY
- MATH_POISSON_HOOK
equation_ids:
- EQ_PROJECTED_CONTINUITY
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
source_ids:
- FILE_4D_PARENT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_LONGITUDINAL_IDENTITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: MATH_POISSON_HOOK
  relation: REDUCES_TO
  status: controlled
  note: Quasi-static longitudinal regime gives Poisson hook.
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: SUPPORTS_CLAIM
  status: exact_projection_plus_controlled_hook
  note: Equation anchor supports this named claim.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references EQ_LONGITUDINAL_IDENTITY.
- source: FILE_4D_PARENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_LONGITUDINAL_IDENTITY.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_LONGITUDINAL_IDENTITY.
- source: EQ_PROJECTED_CONTINUITY
  relation: FEEDS
  status: exact
  note: Longitudinal Helmholtz split of brane velocity gives identity.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact
- topic/moving_throat
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Exact longitudinal identity

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_LONGITUDINAL_IDENTITY`
> **Status:** `exact`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Exact projected identity that contains the Poisson-hook Laplacian.

## Physical Meaning

Exact projected identity that contains the Poisson-hook Laplacian.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact`
- Parent node: [[MATH_LONGITUDINAL_IDENTITY]]

## Equation

$$
ρ_brane ∇_3² φ = S_leak - ∂_t ρ_brane - (∇_3ρ_brane)·(∇_3φ+v_T)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_LONGITUDINAL_IDENTITY]]
- [[MATH_POISSON_HOOK]]

### Related equations
- [[EQ_PROJECTED_CONTINUITY]]

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

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
| `ANCHORS` | [[MATH_LONGITUDINAL_IDENTITY]] | Equation anchor belongs to or formalizes this graph node. |
| `REDUCES_TO` | [[MATH_POISSON_HOOK]] | Quasi-static longitudinal regime gives Poisson hook. |
| `SUPPORTS_CLAIM` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references EQ_LONGITUDINAL_IDENTITY. |
| `CONTAINS_EQUATION` | [[FILE_4D_PARENT]] | Source artifact contains or supports EQ_LONGITUDINAL_IDENTITY. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_LONGITUDINAL_IDENTITY. |
| `FEEDS` | [[EQ_PROJECTED_CONTINUITY]] | Longitudinal Helmholtz split of brane velocity gives identity. |

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
- Line: `907`
- Heading: Exact longitudinal identity from projected continuity
- Nearest label: `sec:poisson-exact` at line `908`
- Match basis: `semantic_heading_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
