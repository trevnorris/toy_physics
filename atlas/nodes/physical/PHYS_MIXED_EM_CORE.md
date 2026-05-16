---
id: PHYS_MIXED_EM_CORE
title: Mixed brane-bulk EM core
type: field_sector
layer: physical_ontology
status: exact_parent_observables_reduced_open
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Hidden brane-bulk EM transport channels retained microscopically outside strict zero-mode brane reduction.
future_paper_needed: false
source_files:
- research/4d_em_fields/paper/4d_em_fields.tex
- research/4d_plasma/paper/4d_plasma.tex
- notes/pde_audit_full.md
- 4d_em_fields_summary.md
- 4d_plasma_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_em_fields_summary.md
- 4d_plasma_summary.md
- pde_audit_full.md:V2-30
source_links:
- '[[FILE_PLASMA]]'
physical_ids:
- MT_V2_30_EM_ONTOLOGY
- PHYS_REG_CHARGE_EM
math_ids:
- MATH_MIXED_FIELDS_EW_CA
- MATH_ZERO_MODE_MAXWELL
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
- CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
- CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
- CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
- CLAIM_MIXED_RECIRCULATION_OPEN
- CLAIM_MIXED_SECTOR_MICROSCOPIC
- CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
- CLAIM_ZERO_MODE_MAXWELL_REDUCTION
open_gate_ids:
- OPEN_LEPTON_SPIN_DISCRETIZER
status_firewall_ids:
- FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE
source_ids:
- FILE_PLASMA
outgoing_edges:
- target: READOUT_D0_C_P0_N2_N4
  relation: AFFECTS
  status: reduced/open
  note: Mixed spectra contribute Z_n,N_n moments.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_reduced_residual
  note: Physical ontology object grounded by this claim.
- target: CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL
  relation: GROUNDS_PHYSICAL_MEANING
  status: conditional_open
  note: Physical ontology object grounded by this claim.
- target: CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH
  relation: GROUNDS_PHYSICAL_MEANING
  status: safe_interpretation_or_structural_patch
  note: Physical ontology object grounded by this claim.
- target: CLAIM_MIXED_RECIRCULATION_OPEN
  relation: GROUNDS_PHYSICAL_MEANING
  status: open
  note: Physical ontology object grounded by this claim.
- target: CLAIM_MIXED_SECTOR_MICROSCOPIC
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_gauge_invariant_with_reduced_uses
  note: Physical ontology object grounded by this claim.
- target: CLAIM_PROJECTED_EM_OUTGOING_BRIDGE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_mixed_kernel
  note: Physical ontology object grounded by this claim.
- target: CLAIM_ZERO_MODE_MAXWELL_REDUCTION
  relation: GROUNDS_PHYSICAL_MEANING
  status: controlled_reduction
  note: Physical ontology object grounded by this claim.
- target: MATH_MIXED_FIELDS_EW_CA
  relation: REPRESENTED_BY
  status: exact
  note: Mixed physical channels represented by gauge-invariant E_w,C_a and A_w/F_mu_w/J^w.
incoming_edges:
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_MIXED_EM_CORE.
- source: BACKLINK_LEPTON_WORK
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_MIXED_EM_CORE.
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references PHYS_MIXED_EM_CORE.
- source: CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
  relation: DEPENDS_ON
  status: conditional_open_plumbing
  note: Plumbing condition lives in the mixed brane-bulk EM channels.
- source: FILE_PLASMA
  relation: DOCUMENTS
  status: source_anchor
  note: File anchor documents this node.
- source: PHYS_REG_CHARGE_EM
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: OPEN_LEPTON_SPIN_DISCRETIZER
  relation: MAY_DEPEND_ON
  status: conditional
  note: Lepton same-charge route likely depends on mixed-sector internal structure.
- source: MT_V2_30_EM_ONTOLOGY
  relation: PATCHES_LANGUAGE
  status: paper_facing
  note: EM ontology/status keeps mixed channels live while respecting zero-mode reduction.
- source: FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: LEPTON_MIXED_BERRY_ROTOR
  relation: REQUIRES
  status: conditional
  note: Same-charge Berry rotor corridor requires mixed core channels to remain live.
- source: MATH_ZERO_MODE_MAXWELL
  relation: SUPPRESSES_BUT_DOES_NOT_REMOVE
  status: firewall
  note: Mixed channels suppressed only in far-field brane limit.
tags:
- atlas/node
- atlas/physical
- layer/physical_ontology
- status/exact_parent_observables_reduced_open
- topic/maxwell
- topic/moving_throat
- topic/projection
- type/field_sector
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Mixed brane-bulk EM core

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `PHYS_MIXED_EM_CORE`
> **Status:** `exact_parent_observables_reduced_open`
> **Layer:** `physical_ontology`
> **Type:** `field_sector`

## Summary

Hidden brane-bulk EM transport channels retained microscopically outside strict zero-mode brane reduction.

## Physical Meaning

Hidden brane-bulk EM transport channels retained microscopically outside strict zero-mode brane reduction.

## Mathematical Role

- Layer: `physical_ontology`
- Type: `field_sector`
- Status: `exact_parent_observables_reduced_open`

## Equation

$$
A_w
$$

$$
J^w
$$

$$
F_{mu w}
$$

$$
E_w
$$

$$
C_a
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[MT_V2_30_EM_ONTOLOGY]]
- [[PHYS_REG_CHARGE_EM]]

### Related math nodes
- [[MATH_MIXED_FIELDS_EW_CA]]
- [[MATH_ZERO_MODE_MAXWELL]]

### Related equations
- none

### Related claims
- [[CLAIM_G2_COMMON_QUOTIENT]]
- [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]]
- [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]]
- [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]]
- [[CLAIM_MIXED_RECIRCULATION_OPEN]]
- [[CLAIM_MIXED_SECTOR_MICROSCOPIC]]
- [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]]
- [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]]

### Open gates
- [[OPEN_LEPTON_SPIN_DISCRETIZER]]

### Status firewalls
- [[FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE]]

### Source anchors
- [[FILE_PLASMA]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `AFFECTS` | [[READOUT_D0_C_P0_N2_N4]] | Mixed spectra contribute Z_n,N_n moments. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_G2_COMMON_QUOTIENT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_LEPTON_HALF_INTEGER_CONDITIONAL]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MAXWELL_GAUGE_LOCALIZATION_PATCH]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MIXED_RECIRCULATION_OPEN]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_MIXED_SECTOR_MICROSCOPIC]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_PROJECTED_EM_OUTGOING_BRIDGE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[CLAIM_ZERO_MODE_MAXWELL_REDUCTION]] | Physical ontology object grounded by this claim. |
| `REPRESENTED_BY` | [[MATH_MIXED_FIELDS_EW_CA]] | Mixed physical channels represented by gauge-invariant E_w,C_a and A_w/F_mu_w/J^w. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references PHYS_MIXED_EM_CORE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_LEPTON_WORK]] | Paper backlink block references PHYS_MIXED_EM_CORE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references PHYS_MIXED_EM_CORE. |
| `DEPENDS_ON` | [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]] | Plumbing condition lives in the mixed brane-bulk EM channels. |
| `DOCUMENTS` | [[FILE_PLASMA]] | File anchor documents this node. |
| `LINKS_TO` | [[PHYS_REG_CHARGE_EM]] | Physical register entry links to graph object. |
| `MAY_DEPEND_ON` | [[OPEN_LEPTON_SPIN_DISCRETIZER]] | Lepton same-charge route likely depends on mixed-sector internal structure. |
| `PATCHES_LANGUAGE` | [[MT_V2_30_EM_ONTOLOGY]] | EM ontology/status keeps mixed channels live while respecting zero-mode reduction. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_ZERO_MODE_NOT_MIXED_ERASURE]] | Firewall preserves this correct status boundary. |
| `REQUIRES` | [[LEPTON_MIXED_BERRY_ROTOR]] | Same-charge Berry rotor corridor requires mixed core channels to remain live. |
| `SUPPRESSES_BUT_DOES_NOT_REMOVE` | [[MATH_ZERO_MODE_MAXWELL]] | Mixed channels suppressed only in far-field brane limit. |

## Source Anchors

### Source anchor notes
- [[FILE_PLASMA]]

### Source files
- `research/4d_em_fields/paper/4d_em_fields.tex`
- `research/4d_plasma/paper/4d_plasma.tex`
- `notes/pde_audit_full.md`
- `4d_em_fields_summary.md`
- `4d_plasma_summary.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
