---
id: MATH_LONGITUDINAL_IDENTITY
title: Exact longitudinal identity
type: identity
layer: math_object
status: exact
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Exact brane scalar-potential identity containing a 3D Laplacian.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- 4d_summary.md
legacy_sources:
- 4d_summary.md
source_links:
- '[[FILE_4D_PARENT]]'
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
- MATH_POISSON_HOOK
- MATH_PROJECTED_CONTINUITY
equation_ids:
- EQ_LONGITUDINAL_IDENTITY
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
source_ids:
- FILE_4D_PARENT
outgoing_edges:
- target: MATH_POISSON_HOOK
  relation: REDUCES_TO
  status: controlled
  note: Quasi-static/longitudinal regime gives Poisson hook.
incoming_edges:
- source: EQ_LONGITUDINAL_IDENTITY
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: MATH_PROJECTED_CONTINUITY
  relation: DERIVES
  status: exact
  note: Helmholtz split yields exact longitudinal identity.
- source: FILE_4D_PARENT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact
- type/identity
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Exact longitudinal identity

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_LONGITUDINAL_IDENTITY`
> **Status:** `exact`
> **Layer:** `math_object`
> **Type:** `identity`

## Summary

Exact brane scalar-potential identity containing a 3D Laplacian.

## Physical Meaning

Exact brane scalar-potential identity containing a 3D Laplacian.

## Mathematical Role

- Layer: `math_object`
- Type: `identity`
- Status: `exact`

## Equation

$$
ρ_brane ∇3²φ = S_leak - ∂tρ_brane - (∇ρ)·(∇φ+v_T)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_POISSON_HOOK]]
- [[MATH_PROJECTED_CONTINUITY]]

### Related equations
- [[EQ_LONGITUDINAL_IDENTITY]]

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_4D_PARENT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REDUCES_TO` | [[MATH_POISSON_HOOK]] | Quasi-static/longitudinal regime gives Poisson hook. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_LONGITUDINAL_IDENTITY]] | Equation anchor belongs to or formalizes this graph node. |
| `DERIVES` | [[MATH_PROJECTED_CONTINUITY]] | Helmholtz split yields exact longitudinal identity. |
| `DOCUMENTS` | [[FILE_4D_PARENT]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]

### Source files
- `research/4d/paper/4d.tex`
- `4d_summary.md`

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
