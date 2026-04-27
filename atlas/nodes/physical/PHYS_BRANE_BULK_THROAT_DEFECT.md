---
id: PHYS_BRANE_BULK_THROAT_DEFECT
title: Finite brane-bulk throat defect
type: defect
layer: physical_ontology
status: physical_ontology
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Particle-like object is a finite-radius puncture/open conduit through the brane into a bulk throat, not a dimple or capped pocket.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- notes/atom_work.md
- pde_audit_full.md
- atom_work.md
legacy_sources:
- pde_audit_full.md:V2-28
- atom_work.md
physical_ids:
- MT_V2_28_ONTOLOGY_CHECKLIST
- PHYS_INTERIOR_SUPPORT
- PHYS_MOUTH_CROSS_SECTION
- PHYS_OPEN_CONDUIT
- PHYS_OPEN_FINITE_EXIT
- PHYS_REG_PARTICLE_THROAT
math_ids:
- MATH_SIGMA_R_FIELD
claim_ids:
- CLAIM_1PN_EIH_WITHIN_HIERARCHY
- CLAIM_LEPTON_MASS_REDUCED_LEDGER
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE1_GEOMETRY_LIFT
outgoing_edges:
- target: CLAIM_1PN_EIH_WITHIN_HIERARCHY
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
- target: CLAIM_STAGE1_GEOMETRY_LIFT
  relation: GROUNDS_PHYSICAL_MEANING
  status: effective_geometry_lift
  note: Physical ontology object grounded by this claim.
- target: PHYS_OPEN_CONDUIT
  relation: HAS_BOUNDARY_CLASS
  status: effective branch
  note: Physical branch is open finite-radius conduit.
- target: PHYS_OPEN_FINITE_EXIT
  relation: HAS_BOUNDARY_FEATURE
  status: physical
  note: Open finite exit is part of the updated physical picture.
- target: PHYS_INTERIOR_SUPPORT
  relation: HAS_PART
  status: ontology
  note: Interior carries axial/cavity support.
- target: PHYS_INTERIOR_SUPPORT
  relation: HAS_PART
  status: physical
  note: Interior support/cavity region extends into bulk.
- target: PHYS_MOUTH_CROSS_SECTION
  relation: HAS_PART
  status: ontology
  note: Mouth is brane entrance cross-section.
- target: PHYS_MOUTH_CROSS_SECTION
  relation: HAS_PART
  status: physical
  note: Mouth is the brane-side cross-section.
- target: MATH_SIGMA_R_FIELD
  relation: REPRESENTED_BY
  status: effective closure
  note: Finite throat geometry represented by Sigma=r-R.
incoming_edges:
- source: MT_V2_28_ONTOLOGY_CHECKLIST
  relation: ANCHORS
  status: paper-facing
  note: Physical picture anchors the atlas ontology layer.
- source: ATLAS_RULE_PHYS_TO_MATH
  relation: APPLIES_TO
  status: active
  note: Start with physical defect, then map to Sigma/R, wall action, response readouts.
- source: BACKLINK_1PN_BRIDGE
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_BRANE_BULK_THROAT_DEFECT.
- source: BACKLINK_1PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_BRANE_BULK_THROAT_DEFECT.
- source: PHYS_REG_PARTICLE_THROAT
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: READOUT_D0_C_P0_N2_N4
  relation: NOT_EQUAL_TO
  status: firewall
  note: Readout packet is not a full physical description of the throat.
- source: MT_V2_28_ONTOLOGY_CHECKLIST
  relation: PATCHES_LANGUAGE
  status: paper_facing
  note: Paper-facing physical picture must be used consistently.
- source: NEG_QUERY_QEFF_FROM_THROAT_GEOMETRY
  relation: STARTS_AT
  status: v07
  note: Negative query starts from PHYS_BRANE_BULK_THROAT_DEFECT.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/physical_ontology
- topic/atom
- topic/moving_throat
- type/defect
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Finite brane-bulk throat defect

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_BRANE_BULK_THROAT_DEFECT`  
> **Status:** `physical_ontology`  
> **Layer:** `physical_ontology`  
> **Type:** `defect`

## Summary

Particle-like object is a finite-radius puncture/open conduit through the brane into a bulk throat, not a dimple or capped pocket.

## Physical Meaning

Particle-like object is a finite-radius puncture/open conduit through the brane into a bulk throat, not a dimple or capped pocket.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `defect`
- Status: `physical_ontology`

## Equation

$$
Sigma=0
$$

$$
r=R(Omega,w,t)
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[MT_V2_28_ONTOLOGY_CHECKLIST]]
- [[PHYS_INTERIOR_SUPPORT]]
- [[PHYS_MOUTH_CROSS_SECTION]]
- [[PHYS_OPEN_CONDUIT]]
- [[PHYS_OPEN_FINITE_EXIT]]
- [[PHYS_REG_PARTICLE_THROAT]]

### Related math nodes
- [[MATH_SIGMA_R_FIELD]]

### Related equations
- none

### Related claims
- [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]]
- [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]]
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE1_GEOMETRY_LIFT]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_LEPTON_MASS_REDUCED_LEDGER]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_STAGE1_GEOMETRY_LIFT]] | Physical ontology object grounded by this claim. |
| `HAS_BOUNDARY_CLASS` | [[PHYS_OPEN_CONDUIT]] | Physical branch is open finite-radius conduit. |
| `HAS_BOUNDARY_FEATURE` | [[PHYS_OPEN_FINITE_EXIT]] | Open finite exit is part of the updated physical picture. |
| `HAS_PART` | [[PHYS_INTERIOR_SUPPORT]] | Interior carries axial/cavity support. |
| `HAS_PART` | [[PHYS_INTERIOR_SUPPORT]] | Interior support/cavity region extends into bulk. |
| `HAS_PART` | [[PHYS_MOUTH_CROSS_SECTION]] | Mouth is brane entrance cross-section. |
| `HAS_PART` | [[PHYS_MOUTH_CROSS_SECTION]] | Mouth is the brane-side cross-section. |
| `REPRESENTED_BY` | [[MATH_SIGMA_R_FIELD]] | Finite throat geometry represented by Sigma=r-R. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS` | [[MT_V2_28_ONTOLOGY_CHECKLIST]] | Physical picture anchors the atlas ontology layer. |
| `APPLIES_TO` | [[ATLAS_RULE_PHYS_TO_MATH]] | Start with physical defect, then map to Sigma/R, wall action, response readouts. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_1PN_BRIDGE]] | Paper backlink block references PHYS_BRANE_BULK_THROAT_DEFECT. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_1PN_FULL]] | Paper backlink block references PHYS_BRANE_BULK_THROAT_DEFECT. |
| `LINKS_TO` | [[PHYS_REG_PARTICLE_THROAT]] | Physical register entry links to graph object. |
| `NOT_EQUAL_TO` | [[READOUT_D0_C_P0_N2_N4]] | Readout packet is not a full physical description of the throat. |
| `PATCHES_LANGUAGE` | [[MT_V2_28_ONTOLOGY_CHECKLIST]] | Paper-facing physical picture must be used consistently. |
| `STARTS_AT` | [[NEG_QUERY_QEFF_FROM_THROAT_GEOMETRY]] | Negative query starts from PHYS_BRANE_BULK_THROAT_DEFECT. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `notes/atom_work.md`
- `pde_audit_full.md`
- `atom_work.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
