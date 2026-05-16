---
id: CLAIM_CHARGE_ONTOLOGY_FIREWALL
title: Corrected electric-charge ontology firewall
type: claim
layer: claim_theorem
status: exact_bookkeeping_firewall
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Electric sign is eta_Q, microscopic branch charge is q_star, brane charge is q_eff=q_star/sqrt(Z_int), and circulation is magnetic/vortical rather than electric charge.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
source_links:
- '[[FILE_1PN_FULL]]'
- '[[FILE_2PN_FULL]]'
- '[[FILE_4D_PARENT]]'
- '[[FILE_ATOM_WORK]]'
- '[[FILE_EM_FIELDS]]'
- '[[FILE_LEPTON_WORK]]'
- '[[SEC_1PN_BRIDGE_DICTIONARY]]'
- '[[SEC_4D_CHARGE_BOOKKEEPING]]'
- '[[SEC_ATOM_FOUNDATIONS]]'
- '[[SEC_EM_CHARGE_ONTOLOGY]]'
- '[[SEC_EM_QEFF]]'
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 250
  heading_level: paragraph
  heading: Minimal coupling and covariant derivatives.
  nearest_label: null
  nearby_labels: []
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 1.0
  confidence: medium
  source_anchor_node: SEC_4D_CHARGE_BOOKKEEPING
physical_ids:
- PHYS_CHARGE_BRANCH
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
- PHYS_REG_CHARGE_EM
math_ids:
- MATH_QSTAR_QEFF
equation_ids:
- EQ_QEFF_NORMALIZATION
claim_ids:
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
status_firewall_ids:
- FIREWALL_CHARGE_NOT_CIRCULATION
- FIREWALL_KAPPARHO_NOT_CHARGE
- FIREWALL_QEFF_THICKNESS_NOT_BREATHING
source_ids:
- FILE_1PN_FULL
- FILE_2PN_FULL
- FILE_4D_PARENT
- FILE_ATOM_WORK
- FILE_EM_FIELDS
- FILE_LEPTON_WORK
- SEC_1PN_BRIDGE_DICTIONARY
- SEC_4D_CHARGE_BOOKKEEPING
- SEC_ATOM_FOUNDATIONS
- SEC_EM_CHARGE_ONTOLOGY
- SEC_EM_QEFF
outgoing_edges:
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: CONSTRAINS_NOTATION_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- target: ATOM_HYDROGEN_ZERO_MODE
  relation: FEEDS_OR_STATUS_OF
  status: exact_bookkeeping_firewall
  note: Claim feeds this downstream object, output, or open gate.
- target: LEPTON_MIXED_BERRY_ROTOR
  relation: FEEDS_OR_STATUS_OF
  status: exact_bookkeeping_firewall
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_QSTAR_QEFF
  relation: FEEDS_OR_STATUS_OF
  status: exact_bookkeeping_firewall
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_1PN_BRIDGE_DICTIONARY
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_CHARGE_ONTOLOGY_FIREWALL
- source: SEC_4D_CHARGE_BOOKKEEPING
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: q-star/effective charge notation.
- source: SEC_ATOM_FOUNDATIONS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Atomic/lepton foundations.
- source: SEC_EM_CHARGE_ONTOLOGY
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Charge sign and q_eff update.
- source: SEC_EM_QEFF
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Thickness-controlled charge.
- source: BACKLINK_1PN_BRIDGE
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL.
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL.
- source: BACKLINK_ATOM_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL.
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL.
- source: BACKLINK_LEPTON_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_bookkeeping_firewall
  note: 'Claim class: exact'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_CHARGE_BRANCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_bookkeeping_firewall
  note: Physical ontology object grounded by this claim.
- source: PHYS_MAGNETIC_VORTICAL_CIRCULATION
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_bookkeeping_firewall
  note: Physical ontology object grounded by this claim.
- source: PHYS_REG_CHARGE_EM
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FILE_1PN_FULL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_bookkeeping_firewall
  note: Source artifact anchors this claim.
- source: FILE_2PN_FULL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_bookkeeping_firewall
  note: Source artifact anchors this claim.
- source: FILE_4D_PARENT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_bookkeeping_firewall
  note: Source artifact anchors this claim.
- source: FILE_ATOM_WORK
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_bookkeeping_firewall
  note: Source artifact anchors this claim.
- source: FILE_EM_FIELDS
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_bookkeeping_firewall
  note: Source artifact anchors this claim.
- source: FILE_LEPTON_WORK
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_bookkeeping_firewall
  note: Source artifact anchors this claim.
- source: FIREWALL_CHARGE_NOT_CIRCULATION
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_KAPPARHO_NOT_CHARGE
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: FIREWALL_QEFF_THICKNESS_NOT_BREATHING
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: EQ_QEFF_NORMALIZATION
  relation: SUPPORTS_CLAIM
  status: exact_bookkeeping_firewall
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_bookkeeping_firewall
- topic/atom
- topic/charge
- topic/lepton
- topic/maxwell
- topic/moving_throat
- topic/pn_chain
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Corrected electric-charge ontology firewall

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_CHARGE_ONTOLOGY_FIREWALL`
> **Status:** `exact_bookkeeping_firewall`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

Electric sign is eta_Q, microscopic branch charge is q_star, brane charge is q_eff=q_star/sqrt(Z_int), and circulation is magnetic/vortical rather than electric charge.

## Claim

Electric sign is eta_Q, microscopic branch charge is q_star, brane charge is q_eff=q_star/sqrt(Z_int), and circulation is magnetic/vortical rather than electric charge.

## Physical Meaning

Electric sign is eta_Q, microscopic branch charge is q_star, brane charge is q_eff=q_star/sqrt(Z_int), and circulation is magnetic/vortical rather than electric charge.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_bookkeeping_firewall`
- Outputs: `MATH_QSTAR_QEFF`, `ATOM_HYDROGEN_ZERO_MODE`, `LEPTON_MIXED_BERRY_ROTOR`

## Atlas Links

### Related physical nodes
- [[PHYS_CHARGE_BRANCH]]
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]
- [[PHYS_REG_CHARGE_EM]]

### Related math nodes
- [[MATH_QSTAR_QEFF]]

### Related equations
- [[EQ_QEFF_NORMALIZATION]]

### Related claims
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_CHARGE_NOT_CIRCULATION]]
- [[FIREWALL_KAPPARHO_NOT_CHARGE]]
- [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]]

### Source anchors
- [[FILE_1PN_FULL]]
- [[FILE_2PN_FULL]]
- [[FILE_4D_PARENT]]
- [[FILE_ATOM_WORK]]
- [[FILE_EM_FIELDS]]
- [[FILE_LEPTON_WORK]]
- [[SEC_1PN_BRIDGE_DICTIONARY]]
- [[SEC_4D_CHARGE_BOOKKEEPING]]
- [[SEC_ATOM_FOUNDATIONS]]
- [[SEC_EM_CHARGE_ONTOLOGY]]
- [[SEC_EM_QEFF]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `CONSTRAINS_NOTATION_FOR` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Claim-level dependency added in v0.4. |
| `FEEDS_OR_STATUS_OF` | [[ATOM_HYDROGEN_ZERO_MODE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[LEPTON_MIXED_BERRY_ROTOR]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_QSTAR_QEFF]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_1PN_BRIDGE_DICTIONARY]] | v0.6 section anchor for CLAIM_CHARGE_ONTOLOGY_FIREWALL |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_CHARGE_BOOKKEEPING]] | q-star/effective charge notation. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_ATOM_FOUNDATIONS]] | Atomic/lepton foundations. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_CHARGE_ONTOLOGY]] | Charge sign and q_eff update. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_QEFF]] | Thickness-controlled charge. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_1PN_BRIDGE]] | Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_ATOM_WORK]] | Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_WORK]] | Paper backlink block references CLAIM_CHARGE_ONTOLOGY_FIREWALL. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_CHARGE_BRANCH]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]] | Physical ontology object grounded by this claim. |
| `LINKS_TO` | [[PHYS_REG_CHARGE_EM]] | Physical register entry links to graph object. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_1PN_FULL]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_2PN_FULL]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_4D_PARENT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_ATOM_WORK]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_EM_FIELDS]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_LEPTON_WORK]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_CHARGE_NOT_CIRCULATION]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_KAPPARHO_NOT_CHARGE]] | Firewall preserves this correct status boundary. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_QEFF_THICKNESS_NOT_BREATHING]] | Firewall preserves this correct status boundary. |
| `SUPPORTS_CLAIM` | [[EQ_QEFF_NORMALIZATION]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_1PN_FULL]]
- [[FILE_2PN_FULL]]
- [[FILE_4D_PARENT]]
- [[FILE_ATOM_WORK]]
- [[FILE_EM_FIELDS]]
- [[FILE_LEPTON_WORK]]
- [[SEC_1PN_BRIDGE_DICTIONARY]]
- [[SEC_4D_CHARGE_BOOKKEEPING]]
- [[SEC_ATOM_FOUNDATIONS]]
- [[SEC_EM_CHARGE_ONTOLOGY]]
- [[SEC_EM_QEFF]]

### Source files
- `research/4d/paper/4d.tex`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `250`
- Heading: Minimal coupling and covariant derivatives.
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
