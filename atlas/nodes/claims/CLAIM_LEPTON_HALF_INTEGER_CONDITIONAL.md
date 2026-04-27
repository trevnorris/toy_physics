---
id: CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
title: Lepton half-integer mixed-sector route is conditional
type: claim
layer: claim_theorem
status: conditional_open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: A same-charge mixed-sector Berry rotor can force half-integer quantization only if the selective tau-subbundle, central sign holonomy, and autonomous self-reproducing soliton cl...
future_paper_needed: false
source_links:
- '[[FILE_LEPTON_WORK]]'
- '[[SEC_LEPTON_MIXED_ROTOR]]'
- '[[SEC_LEPTON_VERDICT]]'
physical_ids:
- PHYS_CHARGE_BRANCH
- PHYS_MIXED_EM_CORE
equation_ids:
- EQ_LEPTON_HALF_INTEGER
claim_ids:
- CLAIM_MIXED_SECTOR_MICROSCOPIC
open_gate_ids:
- OPEN_LEPTON_SPIN_DISCRETIZER
status_firewall_ids:
- FIREWALL_LEPTON_CONDITIONAL
source_ids:
- FILE_LEPTON_WORK
- SEC_LEPTON_MIXED_ROTOR
- SEC_LEPTON_VERDICT
outgoing_edges:
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: DEPENDS_ON
  status: active
  note: Claim-level dependency added in v0.4.
- target: LEPTON_HALF_INTEGER_CONDITIONAL
  relation: FEEDS_OR_STATUS_OF
  status: conditional_open
  note: Claim feeds this downstream object, output, or open gate.
- target: OPEN_LEPTON_SPIN_DISCRETIZER
  relation: FEEDS_OR_STATUS_OF
  status: conditional_open
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_LEPTON_MIXED_ROTOR
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mixed-sector Berry rotor.
- source: SEC_LEPTON_VERDICT
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Conditional lepton verdict.
- source: BACKLINK_LEPTON_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: conditional_open
  note: 'Claim class: conditional_theorem'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_CHARGE_BRANCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_open
  note: Physical ontology object grounded by this claim.
- source: PHYS_MIXED_EM_CORE
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_open
  note: Physical ontology object grounded by this claim.
- source: FILE_LEPTON_WORK
  relation: OWNS_OR_ANCHORS_CLAIM
  status: conditional_open
  note: Source artifact anchors this claim.
- source: FIREWALL_LEPTON_CONDITIONAL
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: NEG_QUERY_LEPTON_HALF_INTEGER_PROVED
  relation: STARTS_AT
  status: v07
  note: Negative query starts from CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL.
- source: EQ_LEPTON_HALF_INTEGER
  relation: SUPPORTS_CLAIM
  status: conditional_open
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/conditional_open
- topic/charge
- topic/lepton
- topic/maxwell
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Lepton half-integer mixed-sector route is conditional

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL`  
> **Status:** `conditional_open`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

A same-charge mixed-sector Berry rotor can force half-integer quantization only if the selective tau-subbundle, central sign holonomy, and autonomous self-reproducing soliton closure hold.

## Claim

A same-charge mixed-sector Berry rotor can force half-integer quantization only if the selective tau-subbundle, central sign holonomy, and autonomous self-reproducing soliton closure hold.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

A same-charge mixed-sector Berry rotor can force half-integer quantization only if the selective tau-subbundle, central sign holonomy, and autonomous self-reproducing soliton closure hold.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `conditional_open`
- Outputs: `LEPTON_HALF_INTEGER_CONDITIONAL`, `OPEN_LEPTON_SPIN_DISCRETIZER`

## Atlas Links

### Related physical nodes
- [[PHYS_CHARGE_BRANCH]]
- [[PHYS_MIXED_EM_CORE]]

### Related math nodes
- none

### Related equations
- [[EQ_LEPTON_HALF_INTEGER]]

### Related claims
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]

### Open gates
- [[OPEN_LEPTON_SPIN_DISCRETIZER]]

### Status firewalls
- [[FIREWALL_LEPTON_CONDITIONAL]]

### Source anchors
- [[FILE_LEPTON_WORK]]
- [[SEC_LEPTON_MIXED_ROTOR]]
- [[SEC_LEPTON_VERDICT]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DEPENDS_ON` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Claim-level dependency added in v0.4. |
| `FEEDS_OR_STATUS_OF` | [[LEPTON_HALF_INTEGER_CONDITIONAL]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_LEPTON_SPIN_DISCRETIZER]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_LEPTON_MIXED_ROTOR]] | Mixed-sector Berry rotor. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_LEPTON_VERDICT]] | Conditional lepton verdict. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_WORK]] | Paper backlink block references CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: conditional_theorem |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_CHARGE_BRANCH]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MIXED_EM_CORE]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_LEPTON_WORK]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_LEPTON_CONDITIONAL]] | Firewall preserves this correct status boundary. |
| `STARTS_AT` | [[NEG_QUERY_LEPTON_HALF_INTEGER_PROVED]] | Negative query starts from CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL. |
| `SUPPORTS_CLAIM` | [[EQ_LEPTON_HALF_INTEGER]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_LEPTON_WORK]]
- [[SEC_LEPTON_MIXED_ROTOR]]
- [[SEC_LEPTON_VERDICT]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
