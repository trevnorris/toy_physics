---
id: READOUT_D0_C_P0_N2_N4
title: Reduced response readout packet
type: readout_packet
layer: status_audit
status: compressed_target_readouts
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Low-order projected response moments used to test branch; not full physical description.
future_paper_needed: false
source_files:
- notes/pde_audit_full.md
- pde_audit_full.md
legacy_sources:
- pde_audit_full.md:V2-26
- pde_audit_full.md:V2-28
physical_ids:
- MT_V2_28_ONTOLOGY_CHECKLIST
- PHYS_BRANE_BULK_THROAT_DEFECT
- PHYS_INTERIOR_SUPPORT
- PHYS_MATERIAL_CLOSURE
- PHYS_MIXED_EM_CORE
- PHYS_RESPONSE_READOUTS
claim_ids:
- CLAIM_RESPONSE_READOUT_DISCIPLINE
- CLAIM_STAGE5_P0_TARGET
open_gate_ids:
- TARGET_PACKET_A
outgoing_edges:
- target: MT_V2_28_ONTOLOGY_CHECKLIST
  relation: CONSTRAINED_BY
  status: status
  note: Checklist warns readouts are not physical ontology variables.
- target: PHYS_BRANE_BULK_THROAT_DEFECT
  relation: NOT_EQUAL_TO
  status: firewall
  note: Readout packet is not a full physical description of the throat.
incoming_edges:
- source: PHYS_INTERIOR_SUPPORT
  relation: AFFECTS
  status: reduced/open
  note: BdG/support spectra contribute B_n moments.
- source: PHYS_MATERIAL_CLOSURE
  relation: AFFECTS
  status: open
  note: Material sector can shift response readouts through density/speed/flux.
- source: PHYS_MIXED_EM_CORE
  relation: AFFECTS
  status: reduced/open
  note: Mixed spectra contribute Z_n,N_n moments.
- source: MT_STAGE6_FULL_GROUPED_BUNDLE
  relation: COMPRESSES_TO
  status: target packet
  note: Full bundle exports compressed audit readouts.
- source: CLAIM_RESPONSE_READOUT_DISCIPLINE
  relation: FEEDS_OR_STATUS_OF
  status: paper_facing_ontology_discipline
  note: Claim feeds this downstream object, output, or open gate.
- source: CLAIM_STAGE5_P0_TARGET
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_grouped_bridge
  note: Claim feeds this downstream object, output, or open gate.
- source: PHYS_RESPONSE_READOUTS
  relation: INTERPRETS
  status: status
  note: Low-order readouts are compressed projected quantities, not the object.
- source: NEG_QUERY_READOUTS_ARE_PHYSICAL_THROAT
  relation: STARTS_AT
  status: v07
  note: Negative query starts from READOUT_D0_C_P0_N2_N4.
- source: TARGET_PACKET_A
  relation: USES
  status: target
  note: Packet A is evaluated through reduced response readouts.
tags:
- atlas/audits
- atlas/node
- layer/status_audit
- status/compressed_target_readouts
- topic/moving_throat
- topic/quadrupole
- type/readout_packet
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Reduced response readout packet

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `READOUT_D0_C_P0_N2_N4`  
> **Status:** `compressed_target_readouts`  
> **Layer:** `status_audit`  
> **Type:** `readout_packet`

## Summary

Low-order projected response moments used to test branch; not full physical description.

## Physical Meaning

Low-order projected response moments used to test branch; not full physical description.

## Mathematical Role

- Layer: `status_audit`
- Type: `readout_packet`
- Status: `compressed_target_readouts`

## Equation

$$
D0
$$

$$
C
$$

$$
P0
$$

$$
N2
$$

$$
N4
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[MT_V2_28_ONTOLOGY_CHECKLIST]]
- [[PHYS_BRANE_BULK_THROAT_DEFECT]]
- [[PHYS_INTERIOR_SUPPORT]]
- [[PHYS_MATERIAL_CLOSURE]]
- [[PHYS_MIXED_EM_CORE]]
- [[PHYS_RESPONSE_READOUTS]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_RESPONSE_READOUT_DISCIPLINE]]
- [[CLAIM_STAGE5_P0_TARGET]]

### Open gates
- [[TARGET_PACKET_A]]

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `CONSTRAINED_BY` | [[MT_V2_28_ONTOLOGY_CHECKLIST]] | Checklist warns readouts are not physical ontology variables. |
| `NOT_EQUAL_TO` | [[PHYS_BRANE_BULK_THROAT_DEFECT]] | Readout packet is not a full physical description of the throat. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `AFFECTS` | [[PHYS_INTERIOR_SUPPORT]] | BdG/support spectra contribute B_n moments. |
| `AFFECTS` | [[PHYS_MATERIAL_CLOSURE]] | Material sector can shift response readouts through density/speed/flux. |
| `AFFECTS` | [[PHYS_MIXED_EM_CORE]] | Mixed spectra contribute Z_n,N_n moments. |
| `COMPRESSES_TO` | [[MT_STAGE6_FULL_GROUPED_BUNDLE]] | Full bundle exports compressed audit readouts. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_RESPONSE_READOUT_DISCIPLINE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[CLAIM_STAGE5_P0_TARGET]] | Claim feeds this downstream object, output, or open gate. |
| `INTERPRETS` | [[PHYS_RESPONSE_READOUTS]] | Low-order readouts are compressed projected quantities, not the object. |
| `STARTS_AT` | [[NEG_QUERY_READOUTS_ARE_PHYSICAL_THROAT]] | Negative query starts from READOUT_D0_C_P0_N2_N4. |
| `USES` | [[TARGET_PACKET_A]] | Packet A is evaluated through reduced response readouts. |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `notes/pde_audit_full.md`
- `pde_audit_full.md`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
