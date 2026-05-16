---
id: MATH_Z_PROFILE
title: Localization profile Z(w)
type: kernel
layer: math_object
status: exact_input_profile
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Controls EM localization and effective brane coupling.
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
- research/4d_plasma/paper/4d_plasma.tex
- 4d_em_fields_summary.md
- 4d_plasma_summary.md
legacy_sources:
- 4d_em_fields_summary.md
- 4d_plasma_summary.md
tex_anchor:
  file: research/4d_plasma/paper/4d_plasma.tex
  line: 263
  heading_level: paragraph
  heading: Localization profile Z(w).
  nearest_label:
    name: eq:Zint_def
    line: 268
  nearby_labels:
  - name: eq:Zint_def
    line: 268
  - name: eq:Z_gaussian_def
    line: 275
  match_basis: semantic_heading_match
  match_score: 1.0
  confidence: high
physical_ids:
- PHYS_LOCALIZED_EM_SECTOR
claim_ids:
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
incoming_edges:
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_Z_PROFILE.
- source: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: FEEDS_OR_STATUS_OF
  status: controlled_reduction
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_LOCALIZED_EM_SECTOR
  relation: LOCALIZED_BY
  status: exact
  note: Z controls EM localization and effective coupling.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_input_profile
- type/kernel
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Localization profile Z(w)

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_Z_PROFILE`
> **Status:** `exact_input_profile`
> **Layer:** `math_object`
> **Type:** `kernel`

## Summary

Controls EM localization and effective brane coupling.

## Physical Meaning

Controls EM localization and effective brane coupling.

## Mathematical Role

- Layer: `math_object`
- Type: `kernel`
- Status: `exact_input_profile`

## Equation

$$
Z_int=∫Z(w)dw
$$

$$
μ0_eff=μ0/Z_int
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_LOCALIZED_EM_SECTOR]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references MATH_Z_PROFILE. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Claim feeds this downstream object, output, or open gate. |
| `LOCALIZED_BY` | [[PHYS_LOCALIZED_EM_SECTOR]] | Z controls EM localization and effective coupling. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`
- `research/4d_plasma/paper/4d_plasma.tex`
- `4d_em_fields_summary.md`
- `4d_plasma_summary.md`

### TeX anchor
- File: `research/4d_plasma/paper/4d_plasma.tex`
- Line: `263`
- Heading: Localization profile Z(w).
- Nearest label: `eq:Zint_def` at line `268`
- Match basis: `semantic_heading_match`
- Confidence: `high`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
