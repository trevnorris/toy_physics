---
id: CLAIM_MIXED_SECTOR_MICROSCOPIC
title: Mixed EM sector is microscopic ontology, not gauge artifact
type: claim
layer: claim_theorem
status: exact_gauge_invariant_with_reduced_uses
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: A_w, F_mu_w, J^w, E_w, and C_a are retained in the parent/microscopic ontology and are the natural carriers for leakage, beyond-MHD structure, and outgoing quadrupole transfer.
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
source_links:
- '[[FILE_EM_FIELDS]]'
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[FILE_PDE_AUDIT]]'
- '[[FILE_PLASMA]]'
- '[[SEC_EM_DEVIATIONS]]'
- '[[SEC_LEPTON_MIXED_ROTOR]]'
- '[[SEC_MT_MIXED_INVARIANTS]]'
- '[[SEC_PDE_EM_STATUS]]'
- '[[SEC_PDE_MIXED_KERNEL]]'
- '[[SEC_PLASMA_AW_ROLE]]'
- '[[SEC_PLASMA_MIXED_FIELDS]]'
tex_anchor:
  file: research/4d_em_fields/paper/4d_em_fields.tex
  line: 1348
  heading_level: subsection
  heading: When the 3+1 Maxwell limit holds
  nearest_label: null
  nearby_labels: []
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 0.436
  confidence: low
  source_anchor_node: SEC_EM_DEVIATIONS
physical_ids:
- PHYS_MIXED_EM_CORE
- PHYS_OPEN_CONDUIT
- PHYS_OUTGOING_QUADRUPOLE_PORT
math_ids:
- MATH_MAXWELL_MIXED_KERNEL
- MATH_MIXED_FIELDS_EW_CA
equation_ids:
- EQ_MAXWELL_MIXED_TRANSFER
claim_ids:
- CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
- CLAIM_MIXED_RECIRCULATION_OPEN
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
status_firewall_ids:
- FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE
source_ids:
- FILE_EM_FIELDS
- FILE_MOVING_THROAT_COMPACT
- FILE_PDE_AUDIT
- FILE_PLASMA
- SEC_EM_DEVIATIONS
- SEC_LEPTON_MIXED_ROTOR
- SEC_MT_MIXED_INVARIANTS
- SEC_PDE_EM_STATUS
- SEC_PDE_MIXED_KERNEL
- SEC_PLASMA_AW_ROLE
- SEC_PLASMA_MIXED_FIELDS
outgoing_edges:
- target: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: ENABLES
  status: active
  note: Claim-level dependency added in v0.4.
- target: MATH_MAXWELL_MIXED_KERNEL
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_MIXED_FIELDS_EW_CA
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
- target: MT_STAGE004_020_PROJECTED_MAXWELL
  relation: FEEDS_OR_STATUS_OF
  status: exact_gauge_invariant_with_reduced_uses
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_EM_DEVIATIONS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mixed/KK deviations from 3+1 Maxwell.
- source: SEC_LEPTON_MIXED_ROTOR
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mixed-sector Berry rotor.
- source: SEC_MT_MIXED_INVARIANTS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mixed fields and gauge invariance.
- source: SEC_PDE_EM_STATUS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: EM ontology and recirculation/plumbing status.
- source: SEC_PDE_MIXED_KERNEL
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Mixed kernel audit.
- source: SEC_PLASMA_AW_ROLE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: A_w role in beyond-MHD channels.
- source: SEC_PLASMA_MIXED_FIELDS
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: E_w and C_a definitions.
- source: BACKLINK_EM_FIELDS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC.
- source: BACKLINK_LEPTON_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC.
- source: BACKLINK_MOVING_THROAT_COMPACT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC.
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_gauge_invariant_with_reduced_uses
  note: 'Claim class: exact_plus_reduced_use'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
  relation: DEPENDS_ON
  status: active
  note: Claim-level dependency added in v0.4.
- source: PHYS_MIXED_EM_CORE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_gauge_invariant_with_reduced_uses
  note: Physical ontology object grounded by this claim.
- source: PHYS_OPEN_CONDUIT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_gauge_invariant_with_reduced_uses
  note: Physical ontology object grounded by this claim.
- source: PHYS_OUTGOING_QUADRUPOLE_PORT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_gauge_invariant_with_reduced_uses
  note: Physical ontology object grounded by this claim.
- source: FILE_EM_FIELDS
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Source artifact anchors this claim.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Source artifact anchors this claim.
- source: FILE_PDE_AUDIT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Source artifact anchors this claim.
- source: FILE_PLASMA
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Source artifact anchors this claim.
- source: FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: REFINES_OPEN_EM_REQUIREMENT
  status: active
  note: Claim-level dependency added in v0.4.
- source: QUERY_MIXED_SECTOR
  relation: STARTS_AT
  status: v06
  note: Query validation start node.
- source: EQ_MAXWELL_MIXED_TRANSFER
  relation: SUPPORTS_CLAIM
  status: exact_gauge_invariant_with_reduced_uses
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_gauge_invariant_with_reduced_uses
- topic/maxwell
- topic/moving_throat
- topic/quadrupole
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Mixed EM sector is microscopic ontology, not gauge artifact

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_MIXED_SECTOR_MICROSCOPIC`
> **Status:** `exact_gauge_invariant_with_reduced_uses`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

A_w, F_mu_w, J^w, E_w, and C_a are retained in the parent/microscopic ontology and are the natural carriers for leakage, beyond-MHD structure, and outgoing quadrupole transfer.

## Claim

A_w, F_mu_w, J^w, E_w, and C_a are retained in the parent/microscopic ontology and are the natural carriers for leakage, beyond-MHD structure, and outgoing quadrupole transfer.

## Physical Meaning

A_w, F_mu_w, J^w, E_w, and C_a are retained in the parent/microscopic ontology and are the natural carriers for leakage, beyond-MHD structure, and outgoing quadrupole transfer.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_gauge_invariant_with_reduced_uses`
- Outputs: `MATH_MIXED_FIELDS_EW_CA`, `MATH_MAXWELL_MIXED_KERNEL`, `MT_STAGE004_020_PROJECTED_MAXWELL`

## Atlas Links

### Related physical nodes
- [[PHYS_MIXED_EM_CORE]]
- [[PHYS_OPEN_CONDUIT]]
- [[PHYS_OUTGOING_QUADRUPOLE_PORT]]

### Related math nodes
- [[MATH_MAXWELL_MIXED_KERNEL]]
- [[MATH_MIXED_FIELDS_EW_CA]]

### Related equations
- [[EQ_MAXWELL_MIXED_TRANSFER]]

### Related claims
- [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE]]

### Source anchors
- [[FILE_EM_FIELDS]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[FILE_PLASMA]]
- [[SEC_EM_DEVIATIONS]]
- [[SEC_LEPTON_MIXED_ROTOR]]
- [[SEC_MT_MIXED_INVARIANTS]]
- [[SEC_PDE_EM_STATUS]]
- [[SEC_PDE_MIXED_KERNEL]]
- [[SEC_PLASMA_AW_ROLE]]
- [[SEC_PLASMA_MIXED_FIELDS]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `ENABLES` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Claim-level dependency added in v0.4. |
| `FEEDS_OR_STATUS_OF` | [[MATH_MAXWELL_MIXED_KERNEL]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_MIXED_FIELDS_EW_CA]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MT_STAGE004_020_PROJECTED_MAXWELL]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_EM_DEVIATIONS]] | Mixed/KK deviations from 3+1 Maxwell. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_LEPTON_MIXED_ROTOR]] | Mixed-sector Berry rotor. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_MIXED_INVARIANTS]] | Mixed fields and gauge invariance. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_EM_STATUS]] | EM ontology and recirculation/plumbing status. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PDE_MIXED_KERNEL]] | Mixed kernel audit. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PLASMA_AW_ROLE]] | A_w role in beyond-MHD channels. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PLASMA_MIXED_FIELDS]] | E_w and C_a definitions. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_EM_FIELDS]] | Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_WORK]] | Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_MOVING_THROAT_COMPACT]] | Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references CLAIM_MIXED_SECTOR_MICROSCOPIC. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_plus_reduced_use |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `DEPENDS_ON` | [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]] | Claim-level dependency added in v0.4. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_MIXED_EM_CORE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_OPEN_CONDUIT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_OUTGOING_QUADRUPOLE_PORT]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_EM_FIELDS]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PDE_AUDIT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PLASMA]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE]] | Firewall preserves this correct status boundary. |
| `REFINES_OPEN_EM_REQUIREMENT` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Claim-level dependency added in v0.4. |
| `STARTS_AT` | [[QUERY_MIXED_SECTOR]] | Query validation start node. |
| `SUPPORTS_CLAIM` | [[EQ_MAXWELL_MIXED_TRANSFER]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_EM_FIELDS]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[FILE_PDE_AUDIT]]
- [[FILE_PLASMA]]
- [[SEC_EM_DEVIATIONS]]
- [[SEC_LEPTON_MIXED_ROTOR]]
- [[SEC_MT_MIXED_INVARIANTS]]
- [[SEC_PDE_EM_STATUS]]
- [[SEC_PDE_MIXED_KERNEL]]
- [[SEC_PLASMA_AW_ROLE]]
- [[SEC_PLASMA_MIXED_FIELDS]]

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`

### TeX anchor
- File: `research/4d_em_fields/paper/4d_em_fields.tex`
- Line: `1348`
- Heading: When the 3+1 Maxwell limit holds
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `low`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
