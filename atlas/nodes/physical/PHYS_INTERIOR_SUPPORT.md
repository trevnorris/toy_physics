---
id: PHYS_INTERIOR_SUPPORT
title: Interior throat/support region
type: geometry_feature
layer: physical_ontology
status: physical_ontology
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Axial/cavity support structure extending into w and carrying wall, BdG, Maxwell/mixed, and outgoing-port degrees of freedom.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/moving_throat_notes_full.md
- pde_audit_full.md
- moving_throat_output_full.md
legacy_sources:
- pde_audit_full.md:V2-28
- moving_throat_output_full.md
physical_ids:
- PHYS_BRANE_BULK_THROAT_DEFECT
- PHYS_SUPERFLUID_INTAKE_OUTPUT
claim_ids:
- CLAIM_2PN_ADM_WITHIN_HIERARCHY
- CLAIM_LEPTON_MASS_REDUCED_LEDGER
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE3_BDG_SCHUR
- CLAIM_STAGE8_TO_14_SELECTED_BRANCH
outgoing_edges:
- target: READOUT_D0_C_P0_N2_N4
  relation: AFFECTS
  status: reduced/open
  note: BdG/support spectra contribute B_n moments.
- target: CLAIM_2PN_ADM_WITHIN_HIERARCHY
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_assembly_within_closure
  note: Physical ontology object grounded by this claim.
- target: CLAIM_LEPTON_MASS_REDUCED_LEDGER
  relation: GROUNDS_PHYSICAL_MEANING
  status: reduced_closure_theorem_open_scale
  note: Physical ontology object grounded by this claim.
- target: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: GROUNDS_PHYSICAL_MEANING
  status: paper_facing_ontology_discipline
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE3_BDG_SCHUR
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_stable_modes
  note: Physical ontology object grounded by this claim.
- target: CLAIM_STAGE8_TO_14_SELECTED_BRANCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_selected_reduced_branch
  note: Physical ontology object grounded by this claim.
- target: PHYS_SUPERFLUID_INTAKE_OUTPUT
  relation: SUPPORTS
  status: physical
  note: Interior/open conduit hosts intake/output and support channels.
incoming_edges:
- source: BACKLINK_1PN_BRIDGE
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_INTERIOR_SUPPORT.
- source: BACKLINK_LEPTON_MASS
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_INTERIOR_SUPPORT.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: HAS_PART
  status: ontology
  note: Interior carries axial/cavity support.
- source: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: HAS_PART
  status: physical
  note: Interior support/cavity region extends into bulk.
- source: LEPTON_MASS_REDUCED_LEDGER
  relation: INTERPRETS
  status: reduced
  note: Mass ledger is internal support/rest-energy branch data, not primitive mass.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/physical_ontology
- topic/maxwell
- topic/moving_throat
- type/geometry_feature
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Interior throat/support region

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_INTERIOR_SUPPORT`  
> **Status:** `physical_ontology`  
> **Layer:** `physical_ontology`  
> **Type:** `geometry_feature`

## Summary

Axial/cavity support structure extending into w and carrying wall, BdG, Maxwell/mixed, and outgoing-port degrees of freedom.

## Physical Meaning

Axial/cavity support structure extending into w and carrying wall, BdG, Maxwell/mixed, and outgoing-port degrees of freedom.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `geometry_feature`
- Status: `physical_ontology`

## Equation

$$
0<=w<=L
$$

$$
support modes
$$

$$
D/N ladder
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]
- [[PHYS_SUPERFLUID_INTAKE_OUTPUT]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]]
- [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE3_BDG_SCHUR]]
- [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `AFFECTS` | [[READOUT_D0_C_P0_N2_N4]] | BdG/support spectra contribute B_n moments. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_2PN_ADM_WITHIN_HIERARCHY]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE3_BDG_SCHUR]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE8_TO_14_SELECTED_BRANCH]] | Physical ontology object grounded by this claim. |
| `SUPPORTS` | [[PHYS_SUPERFLUID_INTAKE_OUTPUT]] | Interior/open conduit hosts intake/output and support channels. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_1PN_BRIDGE]] | Paper backlink block references PHYS_INTERIOR_SUPPORT. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_MASS]] | Paper backlink block references PHYS_INTERIOR_SUPPORT. |
| `HAS_PART` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Interior carries axial/cavity support. |
| `HAS_PART` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Interior support/cavity region extends into bulk. |
| `INTERPRETS` | [[LEPTON_MASS_REDUCED_LEDGER]] | Mass ledger is internal support/rest-energy branch data, not primitive mass. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/moving_throat_notes_full.md`
- `pde_audit_full.md`
- `moving_throat_output_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
