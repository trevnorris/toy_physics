---
id: MATH_FLUXOID
title: Fluxoid/circulation law
type: topological_identity
layer: math_object
status: exact_identity_open_plumbing
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-05-16T02:18:48Z'
summary_short: Topological magnetic/vortical class; not electric charge definition.
future_paper_needed: false
source_files:
- research/4d/paper/4d.tex
- notes/pde_audit_full.md
- notes/circulation/step_01_fluxoid_firewall.md
- notes/circulation/step_01_fluxoid_firewall_sympy.py
- 4d_summary.md
- pde_audit_full.md
legacy_sources:
- 4d_summary.md
- pde_audit_full.md:V2-30
tex_anchor:
  file: research/4d/paper/4d.tex
  line: 725
  heading_level: subsection
  heading: Quantized phase winding and fluxoid quantization
  nearest_label:
    name: eq:circulation-fluxoid
    line: 725
  nearby_labels:
  - name: eq:circulation-fluxoid
    line: 725
  match_basis: semantic_label_match
  match_score: 0.647
  confidence: medium
physical_ids:
- PHYS_CHARGE_BRANCH
- PHYS_MAGNETIC_VORTICAL_CIRCULATION
claim_ids:
- CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION
- CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
- CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
outgoing_edges:
- target: PHYS_CHARGE_BRANCH
  relation: SHOULD_NOT_IDENTIFY_WITH
  status: firewall
  note: Circulation is magnetic/vortical, not electric charge.
incoming_edges:
- source: CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL
  relation: IMPORTS
  status: conditional_open_plumbing
  note: Effective current map starts from the fluxoid/circulation integer but does not derive from it alone.
- source: CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID
  relation: LIMITS_SCOPE_OF
  status: exact_negative_within_3d_closure_analysis
  note: Fluxoid holonomy alone is not a 3D radial-force law.
- source: PHYS_MAGNETIC_VORTICAL_CIRCULATION
  relation: REPRESENTED_BY
  status: exact/open plumbing
  note: Fluxoid/circulation represented by quantized holonomy law.
- source: CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION
  relation: SUPPORTS
  status: exact_identity_audit
  note: Circulation package derives integer phase winding from single-valued psi and supports the fluxoid/circulation law.
tags:
- atlas/math
- atlas/node
- layer/math_object
- status/exact_identity_open_plumbing
- topic/charge
- topic/moving_throat
- type/topological_identity
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Fluxoid/circulation law

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `MATH_FLUXOID`
> **Status:** `exact_identity_open_plumbing`
> **Layer:** `math_object`
> **Type:** `topological_identity`

## Summary

Topological magnetic/vortical class; not electric charge definition.

## Physical Meaning

Topological magnetic/vortical class; not electric charge definition.

## Mathematical Role

- Layer: `math_object`
- Type: `topological_identity`
- Status: `exact_identity_open_plumbing`

## Equation

$$
∮(∂θ-q_*A/ℏ)·dl=2πn
$$

## Variable Dictionary

The graph currently records the equation text but not a full variable dictionary for this generated note. Add dictionary detail in the graph source before regenerating if this equation becomes a reusable paper-facing object.

## Atlas Links

### Related physical nodes
- [[PHYS_CHARGE_BRANCH]]
- [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]]

### Related math nodes
- none

### Related equations
- none

### Related claims
- [[CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION]]
- [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]]
- [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]]

### Open gates
- none

### Status firewalls
- none

### Source anchors
- none

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `SHOULD_NOT_IDENTIFY_WITH` | [[PHYS_CHARGE_BRANCH]] | Circulation is magnetic/vortical, not electric charge. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `IMPORTS` | [[CLAIM_MIXED_CIRCULATION_PLUMBING_CONDITIONAL]] | Effective current map starts from the fluxoid/circulation integer but does not derive from it alone. |
| `LIMITS_SCOPE_OF` | [[CLAIM_NO_UNIVERSAL_FORCE_FROM_FLUXOID]] | Fluxoid holonomy alone is not a 3D radial-force law. |
| `REPRESENTED_BY` | [[PHYS_MAGNETIC_VORTICAL_CIRCULATION]] | Fluxoid/circulation represented by quantized holonomy law. |
| `SUPPORTS` | [[CLAIM_FLUXOID_SINGLE_VALUED_PSI_QUANTIZATION]] | Circulation package derives integer phase winding from single-valued psi and supports the fluxoid/circulati... |

## Source Anchors

### Source anchor notes
- No source anchor note recorded.

### Source files
- `research/4d/paper/4d.tex`
- `notes/pde_audit_full.md`
- `notes/circulation/step_01_fluxoid_firewall.md`
- `notes/circulation/step_01_fluxoid_firewall_sympy.py`
- `4d_summary.md`
- `pde_audit_full.md`

### TeX anchor
- File: `research/4d/paper/4d.tex`
- Line: `725`
- Heading: Quantized phase winding and fluxoid quantization
- Nearest label: `eq:circulation-fluxoid` at line `725`
- Match basis: `semantic_label_match`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
