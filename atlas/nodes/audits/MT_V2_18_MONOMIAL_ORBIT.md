---
id: MT_V2_18_MONOMIAL_ORBIT
title: V2-18 monomial quotient and similarity-orbit audit
type: audit_gate
layer: status_audit
status: exact_within_coherent_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Derives q_tr, q_nt, q_eta quotient coordinates, similarity orbit, and Xi1 bridge.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
math_ids:
- MATH_MONOMIAL_QUOTIENT
status_firewall_ids:
- FIREWALL_SIMILARITY_NOT_FULL_5PN
outgoing_edges:
- target: MATH_MONOMIAL_QUOTIENT
  relation: DERIVES
  status: exact_within_branch
  note: Derives quotient coordinates and similarity-orbit split.
incoming_edges:
- source: FIREWALL_SIMILARITY_NOT_FULL_5PN
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_SIMILARITY_ORBIT_SOLVES_FULL_5PN
  relation: STARTS_AT
  status: v07
  note: Negative query starts from MT_V2_18_MONOMIAL_ORBIT.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/exact_within_coherent_branch
- topic/moving_throat
- type/audit_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# V2-18 monomial quotient and similarity-orbit audit

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_V2_18_MONOMIAL_ORBIT`
> **Status:** `exact_within_coherent_branch`
> **Layer:** `status_audit`
> **Type:** `audit_gate`

## Summary

Derives q_tr, q_nt, q_eta quotient coordinates, similarity orbit, and Xi1 bridge.

## Physical Meaning

Derives q_tr, q_nt, q_eta quotient coordinates, similarity orbit, and Xi1 bridge.

## Mathematical Role

- Layer: `status_audit`
- Type: `audit_gate`
- Status: `exact_within_coherent_branch`

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MONOMIAL_QUOTIENT]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- [[FIREWALL_SIMILARITY_NOT_FULL_5PN]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DERIVES` | [[MATH_MONOMIAL_QUOTIENT]] | Derives quotient coordinates and similarity-orbit split. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `PROTECTS_STATUS_OF` | [[FIREWALL_SIMILARITY_NOT_FULL_5PN]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_SIMILARITY_ORBIT_SOLVES_FULL_5PN]] | Negative query starts from MT_V2_18_MONOMIAL_ORBIT. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
