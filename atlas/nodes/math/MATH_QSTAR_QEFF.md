---
id: MATH_QSTAR_QEFF
title: Corrected charge dictionary
type: ontology_dictionary
layer: math_object
status: exact_bookkeeping
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Separates branch sign, microscopic charge, and brane-observed localized charge.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- research/4d_em_fields/paper/4d_em_fields.tex
- notes/pde_audit_full.md
- 4d_summary.md
- 4d_em_fields_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- 4d_em_fields_summary.md
- pde_audit_full.md:V2-30
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_EM_FIELDS]]'
physical_ids:
- PHYS_CHARGE_BRANCH
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
- PHYS_REG_CHARGE_EM
equation_ids:
- EQ_QEFF_NORMALIZATION
claim_ids:
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
status_firewall_ids:
- FIREWALL_CHARGE_NOT_CIRCULATION
- FIREWALL_QEFF_THICKNESS_NOT_BREATHING
source_ids:
- FILE_4D_PARENT
- FILE_EM_FIELDS
outgoing_edges:
- target: PHYS_MAGNETIC_VORTICAL_CIRCULATION
  relation: DISTINGUISHES_FROM
  status: firewall
  note: Charge dictionary must not be conflated with circulation.
incoming_edges:
- source: EQ_QEFF_NORMALIZATION
  relation: ANCHORS
  status: equation_anchor
  note: Equation anchor belongs to or formalizes this graph node.
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references MATH_QSTAR_QEFF.
- source: FILE_4D_PARENT
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: FILE_EM_FIELDS
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: FEEDS_OR_STATUS_OF
  status: exact_bookkeeping_firewall
  note: Claim feeds this downstream object, output, or open gate.
- source: ATOM_HYDROGEN_ZERO_MODE
  relation: IMPORTS
  status: reduced
  note: Hydrogenic charge scale uses corrected thickness-controlled q_eff.
- source: PHYS_REG_CHARGE_EM
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FIREWALL_CHARGE_NOT_CIRCULATION
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_QEFF_THICKNESS_NOT_BREATHING
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: PHYS_CHARGE_BRANCH
  relation: REPRESENTED_BY
  status: exact
  note: Charge branch represented by eta_Q,q_*,q_eff dictionary.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_bookkeeping
- topic/charge
- topic/moving_throat
- type/ontology_dictionary
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Corrected charge dictionary

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_QSTAR_QEFF`  
> **Status:** `exact_bookkeeping`  
> **Layer:** `math_object`  
> **Type:** `ontology_dictionary`

## Summary

Separates branch sign, microscopic charge, and brane-observed localized charge.

## Physical Meaning

Separates branch sign, microscopic charge, and brane-observed localized charge.

## Mathematical Role

- Layer: `math_object`
- Type: `ontology_dictionary`
- Status: `exact_bookkeeping`

## Equation

$$
eta_Q=±1
$$

$$
q_*=eta_Q e_*
$$

$$
q_eff=q_*/sqrt(Z_int)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_CHARGE_BRANCH]]
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]
- [[PHYS_REG_CHARGE_EM]]

### Related math nodes
- none

### Related equations
- [[EQ_QEFF_NORMALIZATION]]

### Related claims
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_CHARGE_NOT_CIRCULATION]]
- [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]]

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `DISTINGUISHES_FROM` | [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]] | Charge dictionary must not be conflated with circulation. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[EQ_QEFF_NORMALIZATION]] | Equation anchor belongs to or formalizes this graph node. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references MATH_QSTAR_QEFF. |
| `DOCUMENTS` | [[FILE_4D_PARENT]] | File anchor documents this node. |
| `DOCUMENTS` | [[FILE_EM_FIELDS]] | File anchor documents this node. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Claim feeds this downstream object, output, or open gate. |
| `IMPORTS` | [[ATOM_HYDROGEN_ZERO_MODE]] | Hydrogenic charge scale uses corrected thickness-controlled q_eff. |
| `LINKS_TO` | [[PHYS_REG_CHARGE_EM]] | Physical register entry links to graph object. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_CHARGE_NOT_CIRCULATION]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]] | Firewall preserves this correct status boundary. |
| `REPRESENTED_BY` | [[PHYS_CHARGE_BRANCH]] | Charge branch represented by eta_Q,q_*,q_eff dictionary. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_EM_FIELDS]]

### Source files
- `research/4d/paper/4d.tex`
- `research/4d_em_fields/paper/4d_em_fields.tex`
- `notes/pde_audit_full.md`
- `4d_summary.md`
- `4d_em_fields_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
