---
id: CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
title: Mixed-sector circulation force sign remains Lambda_A-conditional
type: claim
layer: claim_theorem
status: conditional_open_plumbing
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Step 06 verifies that N0 is a nonnegative square-over-square transfer magnitude and includes it in the force, but Lambda_A remains real/open; facing-mouth opposite local swirl a...
future_paper_needed: false
source_links:
- '[[FILE_CIRCULATION_PACKAGE]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
- PHYS_MIXED_EM_CORE
math_ids:
- MATH_FLUXOID
equation_ids:
- EQ_MAXWELL_MIXED_TRANSFER
claim_ids:
- CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
- CLAIM_MIXED_RECIRCULATION_OPEN
open_gate_ids:
- OPEN_MIXED_RECIRCULATION
source_ids:
- FILE_CIRCULATION_PACKAGE
outgoing_edges:
- target: PHYS_MIXED_EM_CORE
  relation: DEPENDS_ON
  status: conditional_open_plumbing
  note: Plumbing condition lives in the mixed brane-bulk EM channels.
- target: OPEN_MIXED_RECIRCULATION
  relation: FEEDS_OR_STATUS_OF
  status: conditional_open_plumbing
  note: Circulation package refines the open gate to include Lambda_A sign and N0 magnitude status.
- target: MATH_FLUXOID
  relation: IMPORTS
  status: conditional_open_plumbing
  note: Effective current map starts from the fluxoid/circulation integer but does not derive from it alone.
- target: CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL
  relation: IMPORTS_SIGN_TABLE_CONDITIONALLY
  status: conditional_open_plumbing
  note: Mixed-sector statement imports the facing-mouth sign rule only after a current-like plumbing map is assumed.
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: PRESERVES_OPEN_GATE
  status: conditional_open_plumbing
  note: Lambda_A sign remains an open plumbing-law datum.
incoming_edges:
- source: FILE_CIRCULATION_PACKAGE
  relation: OWNS_OR_ANCHORS_CLAIM
  status: conditional_open_plumbing
  note: Circulation package Step 06 anchors the Lambda_A-conditional mixed plumbing status.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/conditional_open_plumbing
- topic/maxwell
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Mixed-sector circulation force sign remains Lambda_A-conditional

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL`
> **Status:** `conditional_open_plumbing`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

Step 06 verifies that N0 is a nonnegative square-over-square transfer magnitude and includes it in the force, but Lambda_A remains real/open; facing-mouth opposite local swirl attraction requires Lambda1*Lambda2*N0>0.

## Claim

Step 06 verifies that N0 is a nonnegative square-over-square transfer magnitude and includes it in the force, but Lambda_A remains real/open; facing-mouth opposite local swirl attraction requires Lambda1*Lambda2*N0>0.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

Step 06 verifies that N0 is a nonnegative square-over-square transfer magnitude and includes it in the force, but Lambda_A remains real/open; facing-mouth opposite local swirl attraction requires Lambda1*Lambda2*N0>0.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `conditional_open_plumbing`
- Outputs: `OPEN_MIXED_RECIRCULATION`

## Atlas Links

### Related physical nodes
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- [[MATH_FLUXOID]]

### Related equations
- [[EQ_MAXWELL_MIXED_TRANSFER]]

### Related claims
- [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]

### Open gates
- [[OPEN_MIXED_RECIRCULATION]]

### Status firewalls
- none

### Source anchors
- [[FILE_CIRCULATION_PACKAGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DEPENDS_ON` | [[PHYS_MIXED_EM_CORE]] | Plumbing condition lives in the mixed brane-bulk EM channels. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_MIXED_RECIRCULATION]] | Circulation package refines the open gate to include Lambda_A sign and N0 magnitude status. |
| `IMPORTS` | [[MATH_FLUXOID]] | Effective current map starts from the fluxoid/circulation integer but does not derive from it alone. |
| `IMPORTS_SIGN_TABLE_CONDITIONALLY` | [[CLAIM_FACING_MOUTH_SWIRL_CONDITIONAL]] | Mixed-sector statement imports the facing-mouth sign rule only after a current-like plumbing map is assumed. |
| `PRESERVES_OPEN_GATE` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Lambda_A sign remains an open plumbing-law datum. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_CIRCULATION_PACKAGE]] | Circulation package Step 06 anchors the Lambda_A-conditional mixed plumbing status. |

## Source Anchors

### Source anchor notes
- [[FILE_CIRCULATION_PACKAGE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
