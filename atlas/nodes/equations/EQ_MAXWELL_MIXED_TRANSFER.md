---
id: EQ_MAXWELL_MIXED_TRANSFER
title: Stage 021 one-port outgoing transfer factor
type: equation
layer: equation_anchor
status: exact_within_reduced_bundle
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Positive static transfer factor from the retained Stage 021 reduced one-port normal form; the current ledger treats it as the adapter after projection-first EM Stages 004--020.
future_paper_needed: false
source_files:
- research/pde_ledger/paper/stages/stage_021.tex
- notes/moving_throat_pde_program_compact.md
- moving_throat_pde_stage021_reduced_one_port_normal_form.md
- pde_audit_full.md
legacy_sources:
- moving_throat_pde_stage021_reduced_one_port_normal_form.md
- pde_audit_full.md
source_links:
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
tex_anchor:
  file: research/pde_ledger/paper/stages/stage_021.tex
  line: 43
  heading_level: paragraph
  heading: Outgoing transfer.
  nearest_label:
    name: eq:app-stage021-transfer-factor
    line: 47
  nearby_labels:
  - name: eq:app-stage021-transfer-factor
    line: 47
  - name: eq:app-stage021-n0-positive
    line: 54
  match_basis: semantic_heading_match
  match_score: 0.606
  confidence: medium
math_ids:
- MATH_MAXWELL_MIXED_KERNEL
equation_ids:
- EQ_COMPACT_L2_FINGERPRINT
claim_ids:
- CLAIM_MIXED_RECIRCULATION_OPEN
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
source_ids:
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
outgoing_edges:
- target: MATH_MAXWELL_MIXED_KERNEL
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- target: EQ_COMPACT_L2_FINGERPRINT
  relation: MULTIPLIES
  status: reduced
  note: Mixed transfer carries the outgoing fingerprint to the wall branch.
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: SUPPORTS_CLAIM
  status: open
  note: Equation anchor supports this named claim.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: SUPPORTS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Equation anchor supports this named claim.
- target: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_mixed_kernel
  note: Equation anchor supports this named claim.
incoming_edges:
- source: FILE_MOVING_THROAT_COMPACT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MAXWELL_MIXED_TRANSFER.
- source: FILE_PDE_AUDIT
  relation: CONTAINS_EQUATION
  status: source_anchor
  note: Source artifact contains or supports EQ_MAXWELL_MIXED_TRANSFER.
tags:
- atlas/equations
- atlas/node
- layer/equation_anchor
- status/exact_within_reduced_bundle
- topic/maxwell
- topic/moving_throat
- topic/projection
- type/equation
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Stage 021 one-port outgoing transfer factor

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `EQ_MAXWELL_MIXED_TRANSFER`
> **Status:** `exact_within_reduced_bundle`
> **Layer:** `equation_anchor`
> **Type:** `equation`

## Summary

Positive static transfer factor from the retained Stage 021 reduced one-port normal form; the current ledger treats it as the adapter after projection-first EM Stages 004--020.

## Physical Meaning

Positive static transfer factor from the retained Stage 021 reduced one-port normal form; the current ledger treats it as the adapter after projection-first EM Stages 004--020.

## Mathematical Role

- Layer: `equation_anchor`
- Type: `equation`
- Status: `exact_within_reduced_bundle`
- Parent node: [[MATH_MAXWELL_MIXED_KERNEL]]

## Equation

$$
N_l(0)=[Ω_A² g_W + R g_A]²/[Ω_A²Ω_W²-R²]² ≥ 0
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- none

### Related math nodes
- [[MATH_MAXWELL_MIXED_KERNEL]]

### Related equations
- [[EQ_COMPACT_L2_FINGERPRINT]]

### Related claims
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MATH_MAXWELL_MIXED_KERNEL]] | Equation anchor belongs to or formalizes this graph node. |
| `MULTIPLIES` | [[EQ_COMPACT_L2_FINGERPRINT]] | Mixed transfer carries the outgoing fingerprint to the wall branch. |
| `SUPPORTS_CLAIM` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Equation anchor supports this named claim. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `CONTAINS_EQUATION` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact contains or supports EQ_MAXWELL_MIXED_TRANSFER. |
| `CONTAINS_EQUATION` | [[FILE_PDE_AUDIT]] | Source artifact contains or supports EQ_MAXWELL_MIXED_TRANSFER. |

## Source Anchors

### Source anchor notes
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]

### Source files
- `research/pde_ledger/paper/stages/stage_021.tex`
- `notes/moving_throat_pde_program_compact.md`
- `moving_throat_pde_stage021_reduced_one_port_normal_form.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/pde_ledger/paper/stages/stage_021.tex`
- Line: `43`
- Heading: Outgoing transfer.
- Nearest label: `eq:app-stage021-transfer-factor` at line `47`
- Match basis: `semantic_heading_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
