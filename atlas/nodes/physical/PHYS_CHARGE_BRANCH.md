---
id: PHYS_CHARGE_BRANCH
title: Electric charge branch
type: charge_ontology
layer: physical_ontology
status: exact_bookkeeping
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Electric sign is an oriented puncture/branch label; observed charge is localization-dressed.
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
physical_ids:
- MT_V2_30_EM_ONTOLOGY
- PHYS_REG_CHARGE_EM
math_ids:
- MATH_FLUXOID
- MATH_QSTAR_QEFF
claim_ids:
- CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
- CLAIM_CHARGE_ONTOLOGY_FIREWALL
- CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
status_firewall_ids:
- FIREWALL_CHARGE_NOT_CIRCULATION
outgoing_edges:
- target: CLAIM_ATOMIC_HYDROGEN_ZERO_MODE
  relation: GROUNDS_PHYSICAL_MEANING
  status: reduced_sector_consequence
  note: Physical ontology object grounded by this claim.
- target: CLAIM_CHARGE_ONTOLOGY_FIREWALL
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_bookkeeping_firewall
  note: Physical ontology object grounded by this claim.
- target: CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_open
  note: Physical ontology object grounded by this claim.
- target: MATH_QSTAR_QEFF
  relation: REPRESENTED_BY
  status: exact
  note: Charge branch represented by eta_Q,q_*,q_eff dictionary.
incoming_edges:
- source: MT_V2_30_EM_ONTOLOGY
  relation: ANCHORS
  status: paper-facing
  note: EM ontology anchors charge dictionary.
- source: PHYS_REG_CHARGE_EM
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FIREWALL_CHARGE_NOT_CIRCULATION
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: MATH_FLUXOID
  relation: SHOULD_NOT_IDENTIFY_WITH
  status: firewall
  note: Circulation is magnetic/vortical, not electric charge.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_bookkeeping
- topic/charge
- topic/moving_throat
- type/charge_ontology
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Electric charge branch

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_CHARGE_BRANCH`  
> **Status:** `exact_bookkeeping`  
> **Layer:** `physical_ontology`  
> **Type:** `charge_ontology`

## Summary

Electric sign is an oriented puncture/branch label; observed charge is localization-dressed.

## Physical Meaning

Electric sign is an oriented puncture/branch label; observed charge is localization-dressed.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `charge_ontology`
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
- [[MT_V2_30_EM_ONTOLOGY]]
- [[PHYS_REG_CHARGE_EM]]

### Related math nodes
- [[MATH_FLUXOID]]
- [[MATH_QSTAR_QEFF]]

### Related equations
- none

### Related claims
- [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]]
- [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]]
- [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_CHARGE_NOT_CIRCULATION]]

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_ATOMIC_HYDROGEN_ZERO_MODE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_CHARGE_ONTOLOGY_FIREWALL]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]] | Physical ontology object grounded by this claim. |
| `REPRESENTED_BY` | [[MATH_QSTAR_QEFF]] | Charge branch represented by eta_Q,q_*,q_eff dictionary. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MT_V2_30_EM_ONTOLOGY]] | EM ontology anchors charge dictionary. |
| `LINKS_TO` | [[PHYS_REG_CHARGE_EM]] | Physical register entry links to graph object. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_CHARGE_NOT_CIRCULATION]] | Firewall preserves this correct status boundary. |
| `SHOULD_NOT_IDENTIFY_WITH` | [[MATH_FLUXOID]] | Circulation is magnetic/vortical, not electric charge. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

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
