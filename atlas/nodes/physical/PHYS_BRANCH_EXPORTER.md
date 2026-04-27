---
id: PHYS_BRANCH_EXPORTER
title: Actual branch exporter
type: exporter_requirement
layer: physical_ontology
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The missing physical map that takes a solved stationary open-throat branch to K,M,B_n,Z_n,N_n and the frozen response packet.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md
source_links:
- '[[FILE_MOVING_THROAT_OUTPUT]]'
physical_ids:
- PHYS_REG_BRANCH_EXPORT
- PHYS_SUPERFLUID_INTAKE_OUTPUT
claim_ids:
- CLAIM_BRANCH_EXPORTER_REQUIRED
source_ids:
- FILE_MOVING_THROAT_OUTPUT
outgoing_edges:
- target: CLAIM_BRANCH_EXPORTER_REQUIRED
  relation: GROUNDS_PHYSICAL_MEANING
  status: open_actual_branch
  note: Physical ontology object grounded by this claim.
incoming_edges:
- source: FILE_MOVING_THROAT_OUTPUT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PHYS_SUPERFLUID_INTAKE_OUTPUT
  relation: FEEDS
  status: open
  note: Open-system flux and support data must be exported to response packets.
- source: PHYS_REG_BRANCH_EXPORT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/open
- topic/moving_throat
- topic/quadrupole
- type/exporter_requirement
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Actual branch exporter

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_BRANCH_EXPORTER`  
> **Status:** `open`  
> **Layer:** `physical_ontology`  
> **Type:** `exporter_requirement`

## Summary

The missing physical map that takes a solved stationary open-throat branch to K,M,B_n,Z_n,N_n and the frozen response packet.

## Physical Meaning

The missing physical map that takes a solved stationary open-throat branch to K,M,B_n,Z_n,N_n and the frozen response packet.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `exporter_requirement`
- Status: `open`

## Equation

$$
branch -> (K,M,B_n,Z_n,N_n)
$$

$$
branch -> P0,N2,N4
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_REG_BRANCH_EXPORT]]
- [[PHYS_SUPERFLUID_INTAKE_OUTPUT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_BRANCH_EXPORTER_REQUIRED]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_OUTPUT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_BRANCH_EXPORTER_REQUIRED]] | Physical ontology object grounded by this claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `DOCUMENTS` | [[FILE_MOVING_THROAT_OUTPUT]] | File anchor documents this node. |
| `FEEDS` | [[PHYS_SUPERFLUID_INTAKE_OUTPUT]] | Open-system flux and support data must be exported to response packets. |
| `LINKS_TO` | [[PHYS_REG_BRANCH_EXPORT]] | Physical register entry links to graph object. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_OUTPUT]]

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
