---
id: PHYS_SUPERFLUID_INTAKE_OUTPUT
title: Superfluid intake/output channels
type: open_system_flux
layer: physical_ontology
status: open_physical_closure
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Open throat has flux-like intake/output channels that must not be identified by fiat; they feed leakage, work, and export-kernel questions.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- research/4d_plasma/paper/4d_plasma.tex
- pde_audit_full.md
- 4d_plasma_summary.md
legacy_sources:
- pde_audit_full.md
- 4d_plasma_summary.md
source_links:
- '[[FILE_PLASMA]]'
physical_ids:
- PHYS_BRANCH_EXPORTER
- PHYS_INTERIOR_SUPPORT
claim_ids:
- CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
source_ids:
- FILE_PLASMA
outgoing_edges:
- target: PHYS_BRANCH_EXPORTER
  relation: FEEDS
  status: open
  note: Open-system flux and support data must be exported to response packets.
- target: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_projection_plus_controlled_hook
  note: Physical ontology object grounded by this claim.
incoming_edges:
- source: FILE_PLASMA
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PHYS_INTERIOR_SUPPORT
  relation: SUPPORTS
  status: physical
  note: Interior/open conduit hosts intake/output and support channels.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/open_physical_closure
- topic/moving_throat
- type/open_system_flux
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Superfluid intake/output channels

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_SUPERFLUID_INTAKE_OUTPUT`  
> **Status:** `open_physical_closure`  
> **Layer:** `physical_ontology`  
> **Type:** `open_system_flux`

## Summary

Open throat has flux-like intake/output channels that must not be identified by fiat; they feed leakage, work, and export-kernel questions.

## Physical Meaning

Open throat has flux-like intake/output channels that must not be identified by fiat; they feed leakage, work, and export-kernel questions.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `open_system_flux`
- Status: `open_physical_closure`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANCH_EXPORTER]]
- [[PHYS_INTERIOR_SUPPORT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_PLASMA]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[PHYS_BRANCH_EXPORTER]] | Open-system flux and support data must be exported to response packets. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PROJECTION_OPEN_BRANE_SYSTEM]] | Physical ontology object grounded by this claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_PLASMA]] | File anchor documents this node. |
| `SUPPORTS` | [[PHYS_INTERIOR_SUPPORT]] | Interior/open conduit hosts intake/output and support channels. |

## Source Anchors

### Source anchor notes
- [[FILE_PLASMA]]

### Source files
- `notes/pde_audit_full.md`
- `research/4d_plasma/paper/4d_plasma.tex`
- `pde_audit_full.md`
- `4d_plasma_summary.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
