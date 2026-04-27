---
id: ANOMALY_G2_COMMON_QUOTIENT
title: g-2 common quotient correction
type: anomaly_reduction
layer: derivation
status: conditional_open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Common correction organized by q_tr,q_nt,q_eta monomial quotient; Xi1/P1/P0 bridge supplies microscopic target for residual.
future_paper_needed: false
source_files:
- notes/atom_work.md
- notes/lepton_work.md
- notes/lepton_mass_notes.md
- notes/g2/g2_full_output.md
- atom_work.md
- lepton_work.md
- lepton_mass_notes.md
- g2_full_output.md
legacy_sources:
- atom_work.md
- lepton_work.md
- lepton_mass_notes.md
- g2_full_output.md
source_links:
- '[[FILE_G2_OUTPUT]]'
math_ids:
- MATH_MONOMIAL_QUOTIENT
- MATH_XI1_PREF_SLOPE
equation_ids:
- EQ_G2_COMMON_TANGENT
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
source_ids:
- FILE_G2_OUTPUT
outgoing_edges:
- target: MATH_MONOMIAL_QUOTIENT
  relation: IMPORTS
  status: conditional
  note: g-2 common layer uses the same quotient coordinates.
- target: MATH_XI1_PREF_SLOPE
  relation: IMPORTS
  status: conditional
  note: g-2 residual can be expressed as transfer-shape / outgoing-prefactor slope.
incoming_edges:
- source: EQ_G2_COMMON_TANGENT
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: FILE_G2_OUTPUT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_G2_COMMON_QUOTIENT
  relation: FEEDS_OR_STATUS_OF
  status: conditional_reduced_residual
  note: Claim feeds this downstream object, output, or open gate.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/conditional_open
- topic/atom
- topic/g2
- topic/lepton
- topic/projection
- topic/quadrupole
- type/anomaly_reduction
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# g-2 common quotient correction

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `ANOMALY_G2_COMMON_QUOTIENT`  
> **Status:** `conditional_open`  
> **Layer:** `derivation`  
> **Type:** `anomaly_reduction`

## Summary

Common correction organized by q_tr,q_nt,q_eta monomial quotient; Xi1/P1/P0 bridge supplies microscopic target for residual.

## Physical Meaning

Common correction organized by q_tr,q_nt,q_eta monomial quotient; Xi1/P1/P0 bridge supplies microscopic target for residual.

## Mathematical Role

- Layer: `derivation`
- Type: `anomaly_reduction`
- Status: `conditional_open`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MONOMIAL_QUOTIENT]]
- [[MATH_XI1_PREF_SLOPE]]

### Related equations
- [[EQ_G2_COMMON_TANGENT]]

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
| `IMPORTS` | [[MATH_MONOMIAL_QUOTIENT]] | g-2 common layer uses the same quotient coordinates. |
| `IMPORTS` | [[MATH_XI1_PREF_SLOPE]] | g-2 residual can be expressed as transfer-shape / outgoing-prefactor slope. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_G2_COMMON_TANGENT]] | Equation anchor belongs to or formalizes this graph node. |
| `DOCUMENTS` | [[FILE_G2_OUTPUT]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_G2_COMMON_QUOTIENT]] | Claim feeds this downstream object, output, or open gate. |

## Source Anchors

### Source anchor notes
- [[FILE_G2_OUTPUT]]

### Source files
- `notes/atom_work.md`
- `notes/lepton_work.md`
- `notes/lepton_mass_notes.md`
- `notes/g2/g2_full_output.md`
- `atom_work.md`
- `lepton_work.md`
- `lepton_mass_notes.md`
- `g2_full_output.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
