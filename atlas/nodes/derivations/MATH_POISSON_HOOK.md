---
id: MATH_POISSON_HOOK
title: Controlled Poisson hook
type: controlled_reduction
layer: derivation
status: controlled_reduction
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Quasi-static/longitudinal regime yielding inverse-square scalar sector.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_1pn_full/paper/4d_1pn_full.tex
- notes/pde_audit_full.md
- 4d_summary.md
- 4d_1pn_full_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- 4d_1pn_full_summary.md
- pde_audit_full.md:V2-06/V2-07
tex_anchor:
  file: research/4d_1pn_full/paper/4d_1pn_full.tex
  line: 514
  heading_level: subsection
  heading: The exact longitudinal identity and the Poisson hook
  nearest_label:
    name: sec:poisson-hook
    line: 514
  nearby_labels:
  - name: sec:poisson-hook
    line: 514
  match_basis: semantic_label_match
  match_score: 0.69
  confidence: medium
math_ids:
- MATH_LONGITUDINAL_IDENTITY
equation_ids:
- EQ_LONGITUDINAL_IDENTITY
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
outgoing_edges:
- target: PN_0_NEWTONIAN_HOOK
  relation: FEEDS
  status: controlled
  note: Poisson hook feeds Newtonian point-particle limit.
incoming_edges:
- source: BACKLINK_1PN_BRIDGE
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_POISSON_HOOK.
- source: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
- source: EQ_LONGITUDINAL_IDENTITY
  relation: REDUCES_TO
  status: controlled
  note: Quasi-static longitudinal regime gives Poisson hook.
- source: MATH_LONGITUDINAL_IDENTITY
  relation: REDUCES_TO
  status: controlled
  note: Quasi-static/longitudinal regime gives Poisson hook.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/controlled_reduction
- topic/moving_throat
- topic/pn_chain
- topic/projection
- type/controlled_reduction
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Controlled Poisson hook

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_POISSON_HOOK`
> **Status:** `controlled_reduction`
> **Layer:** `derivation`
> **Type:** `controlled_reduction`

## Summary

Quasi-static/longitudinal regime yielding inverse-square scalar sector.

## Physical Meaning

Quasi-static/longitudinal regime yielding inverse-square scalar sector.

## Mathematical Role

- Layer: `derivation`
- Type: `controlled_reduction`
- Status: `controlled_reduction`

## Equation

$$
∇3²φ ~ source
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_LONGITUDINAL_IDENTITY]]

### Related equations
- [[EQ_LONGITUDINAL_IDENTITY]]

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[PN_0_NEWTONIAN_HOOK]] | Poisson hook feeds Newtonian point-particle limit. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_1PN_BRIDGE]] | Paper backlink block references MATH_POISSON_HOOK. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Claim feeds this downstream object, output, or open gate. |
| `REDUCES_TO` | [[EQ_LONGITUDINAL_IDENTITY]] | Quasi-static longitudinal regime gives Poisson hook. |
| `REDUCES_TO` | [[MATH_LONGITUDINAL_IDENTITY]] | Quasi-static/longitudinal regime gives Poisson hook. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_1pn_full/paper/4d_1pn_full.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `4d_1pn_full_summary.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/4d_1pn_full/paper/4d_1pn_full.tex`
- Line: `514`
- Heading: The exact longitudinal identity and the Poisson hook
- Nearest label: `sec:poisson-hook` at line `514`
- Match basis: `semantic_label_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
