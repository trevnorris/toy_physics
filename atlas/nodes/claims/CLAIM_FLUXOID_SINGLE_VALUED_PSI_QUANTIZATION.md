---
id: CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION
title: Fluxoid integer winding from single-valued psi
type: claim
layer: claim_theorem
status: exact_identity_audit
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: The circulation package derives Delta_theta in 2*pi*Z from psi(2*pi)=psi(0), verifies single-valued chi contributes zero loop winding, and preserves eta_Q as a separate electric...
future_paper_needed: false
source_links:
- '[[FILE_CIRCULATION_PACKAGE]]'
physical_ids:
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
math_ids:
- MATH_FLUXOID
status_firewall_ids:
- FIREWALL_CHARGE_NOT_CIRCULATION
source_ids:
- FILE_CIRCULATION_PACKAGE
outgoing_edges:
- target: FIREWALL_CHARGE_NOT_CIRCULATION
  relation: REINFORCES_FIREWALL
  status: exact_identity_audit
  note: Audit keeps eta_Q/q_* bookkeeping separate from the circulation integer.
- target: MATH_FLUXOID
  relation: SUPPORTS
  status: exact_identity_audit
  note: Circulation package derives integer phase winding from single-valued psi and supports the fluxoid/circulation law.
incoming_edges:
- source: FILE_CIRCULATION_PACKAGE
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_identity_audit
  note: Circulation package Step 01 anchors this quantized phase-winding audit.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_identity_audit
- topic/charge
- topic/moving_throat
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Fluxoid integer winding from single-valued psi

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION`
> **Status:** `exact_identity_audit`
> **Layer:** `claim_theorem`
> **Type:** `claim`

## Summary

The circulation package derives Delta_theta in 2*pi*Z from psi(2*pi)=psi(0), verifies single-valued chi contributes zero loop winding, and preserves eta_Q as a separate electric branch label.

## Claim

The circulation package derives Delta_theta in 2*pi*Z from psi(2*pi)=psi(0), verifies single-valued chi contributes zero loop winding, and preserves eta_Q as a separate electric branch label.

## Physical Meaning

The circulation package derives Delta_theta in 2*pi*Z from psi(2*pi)=psi(0), verifies single-valued chi contributes zero loop winding, and preserves eta_Q as a separate electric branch label.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_identity_audit`
- Outputs: `MATH_FLUXOID`, `FIREWALL_CHARGE_NOT_CIRCULATION`

## Atlas Links

### Related physical nodes
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]

### Related math nodes
- [[MATH_FLUXOID]]

### Related equations
- none

### Related claims
- none

### Open gates
- none

### Status firewalls
- [[FIREWALL_CHARGE_NOT_CIRCULATION]]

### Source anchors
- [[FILE_CIRCULATION_PACKAGE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `REINFORCES_FIREWALL` | [[FIREWALL_CHARGE_NOT_CIRCULATION]] | Audit keeps eta_Q/q_* bookkeeping separate from the circulation integer. |
| `SUPPORTS` | [[MATH_FLUXOID]] | Circulation package derives integer phase winding from single-valued psi and supports the fluxoid/circulati... |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_CIRCULATION_PACKAGE]] | Circulation package Step 01 anchors this quantized phase-winding audit. |

## Source Anchors

### Source anchor notes
- [[FILE_CIRCULATION_PACKAGE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
