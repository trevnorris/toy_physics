---
id: PHYS_SUPERFLUID_MEDIUM
title: Bulk superfluid/GNLS medium
type: medium
layer: physical_ontology
status: exact_declared_parent
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The substrate carrying density, current, EOS, sound speed, and support dynamics.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/pde_audit_full.md
- 4d_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- pde_audit_full.md
math_ids:
- MATH_GNLS_PSI
claim_ids:
- CLAIM_PARENT_ACTION_CURRENT_EXACT
outgoing_edges:
- target: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_parent_current
  note: Physical ontology object grounded by this claim.
- target: MATH_GNLS_PSI
  relation: REPRESENTED_BY
  status: exact
  note: Superfluid medium represented by gauged GNLS field and density/current.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_declared_parent
- topic/moving_throat
- type/medium
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Bulk superfluid/GNLS medium

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_SUPERFLUID_MEDIUM`
> **Status:** `exact_declared_parent`
> **Layer:** `physical_ontology`
> **Type:** `medium`

## Summary

The substrate carrying density, current, EOS, sound speed, and support dynamics.

## Physical Meaning

The substrate carrying density, current, EOS, sound speed, and support dynamics.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `medium`
- Status: `exact_declared_parent`

## Equation

$$
psi
$$

$$
rho=|psi|^2
$$

$$
P(rho)=K rho^5
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_GNLS_PSI]]

### Related equations
- none

### Related claims
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Physical ontology object grounded by this claim. |
| `REPRESENTED_BY` | [[MATH_GNLS_PSI]] | Superfluid medium represented by gauged GNLS field and density/current. |

## Incoming Edges

- none

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
