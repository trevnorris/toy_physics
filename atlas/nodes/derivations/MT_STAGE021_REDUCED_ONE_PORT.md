---
id: MT_STAGE021_REDUCED_ONE_PORT
title: Stage 021 reduced Maxwell/mixed one-port adapter
type: moving_throat_stage
layer: derivation
status: retained_reduced_adapter
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Retains the earlier reduced Maxwell/mixed calculation as the one-port normal form used downstream after the projection-first EM packet.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_021.tex
- notes/moving_throat_pde_program_compact.md
- moving_throat_pde_stage021_reduced_one_port_normal_form.md
legacy_sources:
- moving_throat_pde_stage021_reduced_one_port_normal_form.md
claim_ids:
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
outgoing_edges:
- target: MT_STAGE022_GROUPED_P2_BRIDGE
  relation: FEEDS
  status: exact_closure_adapter
  note: One-lane outgoing transfer lifted to grouped P2 normalization language.
incoming_edges:
- source: MT_STAGE004_020_PROJECTED_MAXWELL
  relation: CLOSED_BY
  status: retained_reduced_adapter
  note: Projection-first EM derivation is closed for downstream response work by the retained one-port normal form.
- source: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: FEEDS_OR_STATUS_OF
  status: retained_reduced_adapter
  note: Retained reduced one-port adapter closes the projection-first EM packet for downstream grouped normalization.
tags:
- atlas/derivations
- atlas/node
- layer/derivation
- status/retained_reduced_adapter
- topic/maxwell
- topic/moving_throat
- topic/projection
- type/moving_throat_stage
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 021 reduced Maxwell/mixed one-port adapter

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MT_STAGE021_REDUCED_ONE_PORT`
> **Status:** `retained_reduced_adapter`
> **Layer:** `derivation`
> **Type:** `moving_throat_stage`

## Summary

Retains the earlier reduced Maxwell/mixed calculation as the one-port normal form used downstream after the projection-first EM packet.

## Physical Meaning

Retains the earlier reduced Maxwell/mixed calculation as the one-port normal form used downstream after the projection-first EM packet.

## Mathematical Role

- Layer: `derivation`
- Type: `moving_throat_stage`
- Status: `retained_reduced_adapter`

## Equation

$$
N_l(0)=[Ω_A²g_W+Rg_A]²/[Ω_A²Ω_W²-R²]²
$$

$$
Y_2^out=1+a²ω²/(9c_s²)+4a⁴ω⁴/(81c_s⁴)+i a⁵ω⁵/(27c_s⁵)+O(ω⁶)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS` | [[MT_STAGE022_GROUPED_P2_BRIDGE]] | One-lane outgoing transfer lifted to grouped P2 normalization language. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CLOSED_BY` | [[MT_STAGE004_020_PROJECTED_MAXWELL]] | Projection-first EM derivation is closed for downstream response work by the retained one-port normal form. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Retained reduced one-port adapter closes the projection-first EM packet for downstream grouped normalization. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/pde_ledger/paper/stages/stage_021.tex`
- `notes/moving_throat_pde_program_compact.md`
- `moving_throat_pde_stage021_reduced_one_port_normal_form.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
