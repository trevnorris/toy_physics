---
id: OPEN_MIXED_RECIRCULATION
title: 'Open gate: mixed recirculation/plumbing law'
type: EM_plumbing_gate
layer: open_gate
status: open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Need closed law connecting circulation, throat intake, mixed transport, brane magnetic fields, and current source.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/lepton_work.md
- notes/circulation/step_06_mixed_sector_plumbing_closure.md
- notes/circulation/step_06_mixed_sector_plumbing_closure_sympy.py
- pde_audit_full.md
- lepton_work.md
legacy_sources:
- pde_audit_full.md:V2-30
- lepton_work.md
source_links:
- '[[SEC_PDE_EM_STATUS]]'
physical_ids:
- MT_V2_30_EM_ONTOLOGY
claim_ids:
- CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
- CLAIM_MIXED_RECIRCULATION_OPEN
source_ids:
- SEC_PDE_EM_STATUS
incoming_edges:
- source: SEC_PDE_EM_STATUS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: EM ontology and recirculation/plumbing status.
- source: CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
  relation: FEEDS_OR_STATUS_OF
  status: conditional_open_plumbing
  note: Circulation package refines the open gate to include Lambda_A sign and N0 magnitude status.
- source: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: FEEDS_OR_STATUS_OF
  status: open
  note: Claim feeds this downstream object, output, or open gate.
- source: BACKLINK_EM_FIELDS
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION.
- source: BACKLINK_LEPTON_WORK
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION.
- source: BACKLINK_PDE_AUDIT
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION.
- source: BACKLINK_PLASMA
  relation: FLAGS_OPEN_GATE
  status: v06
  note: Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION.
- source: ATLAS_CURRENT_READINESS_V05
  relation: FLAGS_REMAINING_GATE
  status: v05
  note: Still open after atlas organization; atlas tracks it but does not solve it.
- source: MT_V2_30_EM_ONTOLOGY
  relation: OPENS_GATE
  status: open
  note: Full recirculation/plumbing law remains future work.
tags:
- atlas/node
- atlas/open_gates
- layer/open_gate
- status/open
- topic/lepton
- topic/maxwell
- topic/moving_throat
- type/em_plumbing_gate
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Open gate: mixed recirculation/plumbing law

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `OPEN_MIXED_RECIRCULATION`
> **Status:** `open`
> **Layer:** `open_gate`
> **Type:** `EM_plumbing_gate`

## Summary

Need closed law connecting circulation, throat intake, mixed transport, brane magnetic fields, and current source.

## What Remains Open

Need closed law connecting circulation, throat intake, mixed transport, brane magnetic fields, and current source.

## What Would Close It

A source-backed derivation, branch computation, theorem, or paper update must change the graph source of truth before this note can change status.

## Physical Meaning

Need closed law connecting circulation, throat intake, mixed transport, brane magnetic fields, and current source.

## Mathematical Role

- Layer: `open_gate`
- Type: `EM_plumbing_gate`
- Status: `open`

## Equation

$$
J^M source
$$

$$
A_w/F_muw/J^w transport
$$

$$
fluxoid
$$

$$
Lambda_A sign/plumbing map
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[MT_V2_30_EM_ONTOLOGY]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[SEC_PDE_EM_STATUS]]

## Outgoing Edges

- none

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_EM_STATUS]] | EM ontology and recirculation/plumbing status. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]] | Circulation package refines the open gate to include Lambda_A sign and N0 magnitude status. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Claim feeds this downstream object, output, or open gate. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_LEPTON_WORK]] | Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_PDE_AUDIT]] | Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION. |
| `FLAGS_OPEN_GATE` | [[BACKLINK_PLASMA]] | Paper backlink block flags open gate OPEN_MIXED_RECIRCULATION. |
| `FLAGS_REMAINING_GATE` | [[ATLAS_CURRENT_READINESS_V05]] | Still open after atlas organization; atlas tracks it but does not solve it. |
| `OPENS_GATE` | [[MT_V2_30_EM_ONTOLOGY]] | Full recirculation/plumbing law remains future work. |

## Source Anchors

### Source anchor notes
- [[SEC_PDE_EM_STATUS]]

### Source files
- `notes/pde_audit_full.md`
- `notes/lepton_work.md`
- `notes/circulation/step_06_mixed_sector_plumbing_closure.md`
- `notes/circulation/step_06_mixed_sector_plumbing_closure_sympy.py`
- `pde_audit_full.md`
- `lepton_work.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
