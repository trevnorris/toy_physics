---
id: FILE_G2_OUTPUT
title: g2_full_output.md
type: source_file
layer: file_anchor
status: anomaly_anchor
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Anomaly/g-2 output, staggered baseline, common quotient path, and outgoing-prefactor bridge.
source_kind: future_paper_note
future_paper_needed: true
source_files:
- notes/g2/g2_full_output.md
- g2_full_output.md
legacy_sources:
- g2_full_output.md
source_links:
- '[[SEC_G2_PREF_SLOPE]]'
- '[[SEC_G2_QUOTIENT]]'
equation_ids:
- EQ_G2_COMMON_TANGENT
- EQ_MONOMIAL_QUOTIENT
- EQ_XI1_PREF_SLOPE
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
status_firewall_ids:
- FIREWALL_G2_COMMON_CONDITIONAL
source_ids:
- SEC_G2_PREF_SLOPE
- SEC_G2_QUOTIENT
outgoing_edges:
- target: EQ_G2_COMMON_TANGENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_G2_COMMON_TANGENT.
- target: EQ_MONOMIAL_QUOTIENT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MONOMIAL_QUOTIENT.
- target: EQ_XI1_PREF_SLOPE
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_XI1_PREF_SLOPE.
- target: ANOMALY_G2_COMMON_QUOTIENT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- target: ANOMALY_G2_STAGGERED_BASELINE
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- target: BACKLINK_G2_OUTPUT
  relation: HAS_PAPER_BACKLINK_BLOCK
  status: v06
  note: Paste-ready atlas backlink block prepared in v0.6.
- target: SEC_G2_PREF_SLOPE
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: g2_full_output.md:1144 — Step 4 — Exact transfer-shape / outgoing-prefactor bridge
- target: SEC_G2_QUOTIENT
  relation: HAS_SECTION_ANCHOR
  status: v05
  note: g2_full_output.md:1 — Step 1 — Exact moving-throat quotient bridge and staggered anomaly benchmark
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: conditional_reduced_residual
  note: Source artifact anchors this claim.
incoming_edges:
- source: FIREWALL_G2_COMMON_CONDITIONAL
  relation: ANCHORED_IN
  status: v07
  note: Firewall is anchored in this source/file/section.
tags:
- atlas/node
- atlas/sources
- layer/file_anchor
- status/anomaly_anchor
- topic/g2
- type/source_file
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# g2_full_output.md

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `FILE_G2_OUTPUT`
> **Status:** `anomaly_anchor`
> **Layer:** `file_anchor`
> **Type:** `source_file`

## Summary

Anomaly/g-2 output, staggered baseline, common quotient path, and outgoing-prefactor bridge.

> [!note] Future paper needed
> This node is intentionally anchored to a notes/derivation source until a maintained paper exists.

## Physical Meaning

Anomaly/g-2 output, staggered baseline, common quotient path, and outgoing-prefactor bridge.

## Mathematical Role

- Layer: `file_anchor`
- Type: `source_file`
- Status: `anomaly_anchor`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- [[EQ_G2_COMMON_TANGENT]]
- [[EQ_MONOMIAL_QUOTIENT]]
- [[EQ_XI1_PREF_SLOPE]]

### Related claims
- [[CLAIM_G2_COMMON_QUOTIENT]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_G2_COMMON_CONDITIONAL]]

### Source anchors
- [[SEC_G2_PREF_SLOPE]]
- [[SEC_G2_QUOTIENT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[EQ_G2_COMMON_TANGENT]] | Source artifact contains or supports EQ_G2_COMMON_TANGENT. |
| `CONTAINS_EQUATION` | [[EQ_MONOMIAL_QUOTIENT]] | Source artifact contains or supports EQ_MONOMIAL_QUOTIENT. |
| `CONTAINS_EQUATION` | [[EQ_XI1_PREF_SLOPE]] | Source artifact contains or supports EQ_XI1_PREF_SLOPE. |
| `DOCUMENTS` | [[ANOMALY_G2_COMMON_QUOTIENT]] | File anchor documents this node. |
| `DOCUMENTS` | [[ANOMALY_G2_STAGGERED_BASELINE]] | File anchor documents this node. |
| `HAS_PAPER_BACKLINK_BLOCK` | [[BACKLINK_G2_OUTPUT]] | Paste-ready atlas backlink block prepared in v0.6. |
| `HAS_SECTION_ANCHOR` | [[SEC_G2_PREF_SLOPE]] | g2_full_output.md:1144 — Step 4 — Exact transfer-shape / outgoing-prefactor bridge |
| `HAS_SECTION_ANCHOR` | [[SEC_G2_QUOTIENT]] | g2_full_output.md:1 — Step 1 — Exact moving-throat quotient bridge and staggered anomaly benchmark |
| `OWNS_OR_ANCHORS_CLAIM` | [[CLAIM_G2_COMMON_QUOTIENT]] | Source artifact anchors this claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORED_IN` | [[FIREWALL_G2_COMMON_CONDITIONAL]] | Firewall is anchored in this source/file/section. |

## Source Anchors

### Source anchor notes
- [[SEC_G2_PREF_SLOPE]]
- [[SEC_G2_QUOTIENT]]

### Source files
- `notes/g2/g2_full_output.md`
- `g2_full_output.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
