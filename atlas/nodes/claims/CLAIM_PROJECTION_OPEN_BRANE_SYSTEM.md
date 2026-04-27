---
id: CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
title: Projection makes the brane an open subsystem
type: claim
layer: claim_theorem
status: exact_projection_plus_controlled_hook
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: Projection with W(w) gives exact projected continuity with leakage; the longitudinal identity becomes the Poisson hook only in a controlled quasi-static/longitudinal regime.
future_paper_needed: false
source_files:
- research/4d_plasma/paper/4d_plasma.tex
source_links:
- '[[FILE_4D_PARENT]]'
- '[[FILE_PLASMA]]'
- '[[SEC_1PN_BRIDGE_DICTIONARY]]'
- '[[SEC_4D_POISSON_HOOK]]'
- '[[SEC_4D_PROJECTED_CONTINUITY]]'
- '[[SEC_4D_PROJECTION_REDUCTION]]'
- '[[SEC_MT_PROJECTION_ZERO_MODE]]'
- '[[SEC_PLASMA_LEAKAGE]]'
- '[[SEC_PLASMA_Z_VS_W]]'
tex_anchor:
  file: research/4d_plasma/paper/4d_plasma.tex
  line: 252
  heading_level: subsection
  heading: Localization and Projection Operators
  nearest_label:
    name: sec:projection_defs
    line: 253
  nearby_labels:
  - name: sec:projection_defs
    line: 253
  match_basis: graph_edge:ANCHORS_CLAIM_SECTION
  match_score: 0.722
  confidence: medium
  source_anchor_node: SEC_PLASMA_Z_VS_W
physical_ids:
- PHYS_BRANE_OBSERVER
- PHYS_OPEN_CONDUIT
- PHYS_REG_OBSERVER_PROJECTION
- PHYS_SUPERFLUID_INTAKE_OUTPUT
math_ids:
- MATH_LONGITUDINAL_IDENTITY
- MATH_POISSON_HOOK
- MATH_PROJECTED_CONTINUITY
equation_ids:
- EQ_LONGITUDINAL_IDENTITY
- EQ_PROJECTED_CONTINUITY
claim_ids:
- CLAIM_1PN_EIH_WITHIN_HIERARCHY
- CLAIM_PARENT_ACTION_CURRENT_EXACT
status_firewall_ids:
- FIREWALL_PROJECTION_NOT_REDUCTION
source_ids:
- FILE_4D_PARENT
- FILE_PLASMA
- SEC_1PN_BRIDGE_DICTIONARY
- SEC_4D_POISSON_HOOK
- SEC_4D_PROJECTED_CONTINUITY
- SEC_4D_PROJECTION_REDUCTION
- SEC_MT_PROJECTION_ZERO_MODE
- SEC_PLASMA_LEAKAGE
- SEC_PLASMA_Z_VS_W
outgoing_edges:
- target: CLAIM_1PN_EIH_WITHIN_HIERARCHY
  relation: BACKBONE_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- target: MATH_LONGITUDINAL_IDENTITY
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_POISSON_HOOK
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
- target: MATH_PROJECTED_CONTINUITY
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
- target: PN_0_NEWTONIAN_HOOK
  relation: FEEDS_OR_STATUS_OF
  status: exact_projection_plus_controlled_hook
  note: Claim feeds this downstream object, output, or open gate.
incoming_edges:
- source: SEC_1PN_BRIDGE_DICTIONARY
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_PROJECTION_OPEN_BRANE_SYSTEM
- source: SEC_4D_POISSON_HOOK
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Poisson hook from longitudinal identity.
- source: SEC_4D_PROJECTED_CONTINUITY
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Exact leakage identity.
- source: SEC_4D_PROJECTION_REDUCTION
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Projection vs reduction firewall.
- source: SEC_MT_PROJECTION_ZERO_MODE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Projection and zero-mode hooks.
- source: SEC_PLASMA_LEAKAGE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Open-system brane identities.
- source: SEC_PLASMA_Z_VS_W
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Z controls action; W controls observation.
- source: BACKLINK_4D_PARENT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_PROJECTION_OPEN_BRANE_SYSTEM.
- source: BACKLINK_PLASMA
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_PROJECTION_OPEN_BRANE_SYSTEM.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_projection_plus_controlled_hook
  note: 'Claim class: exact_and_controlled_reduction'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_BRANE_OBSERVER
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_projection_plus_controlled_hook
  note: Physical ontology object grounded by this claim.
- source: PHYS_OPEN_CONDUIT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_projection_plus_controlled_hook
  note: Physical ontology object grounded by this claim.
- source: PHYS_SUPERFLUID_INTAKE_OUTPUT
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_projection_plus_controlled_hook
  note: Physical ontology object grounded by this claim.
- source: PHYS_REG_OBSERVER_PROJECTION
  relation: LINKS_TO
  status: active
  note: Physical register entry links to graph object.
- source: FILE_4D_PARENT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_projection_plus_controlled_hook
  note: Source artifact anchors this claim.
- source: FILE_PLASMA
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_projection_plus_controlled_hook
  note: Source artifact anchors this claim.
- source: FIREWALL_PROJECTION_NOT_REDUCTION
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: CLAIM_PARENT_ACTION_CURRENT_EXACT
  relation: SUPPLIES_PARENT_FOR
  status: active
  note: Claim-level dependency added in v0.4.
- source: EQ_LONGITUDINAL_IDENTITY
  relation: SUPPORTS_CLAIM
  status: exact_projection_plus_controlled_hook
  note: Equation anchor supports this named claim.
- source: EQ_PROJECTED_CONTINUITY
  relation: SUPPORTS_CLAIM
  status: exact_projection_plus_controlled_hook
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_projection_plus_controlled_hook
- topic/projection
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# Projection makes the brane an open subsystem

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_PROJECTION_OPEN_BRANE_SYSTEM`  
> **Status:** `exact_projection_plus_controlled_hook`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

Projection with W(w) gives exact projected continuity with leakage; the longitudinal identity becomes the Poisson hook only in a controlled quasi-static/longitudinal regime.

## Claim

Projection with W(w) gives exact projected continuity with leakage; the longitudinal identity becomes the Poisson hook only in a controlled quasi-static/longitudinal regime.

## Physical Meaning

Projection with W(w) gives exact projected continuity with leakage; the longitudinal identity becomes the Poisson hook only in a controlled quasi-static/longitudinal regime.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_projection_plus_controlled_hook`
- Outputs: `MATH_PROJECTED_CONTINUITY`, `MATH_LONGITUDINAL_IDENTITY`, `MATH_POISSON_HOOK`, `PN_0_NEWTONIAN_HOOK`

## Atlas Links

### Related physical nodes
- [[PHYS_BRANE_OBSERVER]]
- [[PHYS_OPEN_CONDUIT]]
- [[PHYS_REG_OBSERVER_PROJECTION]]
- [[PHYS_SUPERFLUID_INTAKE_OUTPUT]]

### Related math nodes
- [[MATH_LONGITUDINAL_IDENTITY]]
- [[MATH_POISSON_HOOK]]
- [[MATH_PROJECTED_CONTINUITY]]

### Related equations
- [[EQ_LONGITUDINAL_IDENTITY]]
- [[EQ_PROJECTED_CONTINUITY]]

### Related claims
- [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]]
- [[CLAIM_PARENT_ACTION_CURRENT_EXACT]]

### Open gates
- none

### Status firewalls
- [[FIREWALL_PROJECTION_NOT_REDUCTION]]

### Source anchors
- [[FILE_4D_PARENT]]
- [[FILE_PLASMA]]
- [[SEC_1PN_BRIDGE_DICTIONARY]]
- [[SEC_4D_POISSON_HOOK]]
- [[SEC_4D_PROJECTED_CONTINUITY]]
- [[SEC_4D_PROJECTION_REDUCTION]]
- [[SEC_MT_PROJECTION_ZERO_MODE]]
- [[SEC_PLASMA_LEAKAGE]]
- [[SEC_PLASMA_Z_VS_W]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `BACKBONE_FOR` | [[CLAIM_1PN_EIH_WITHIN_HIERARCHY]] | Claim-level dependency added in v0.4. |
| `FEEDS_OR_STATUS_OF` | [[MATH_LONGITUDINAL_IDENTITY]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_POISSON_HOOK]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[MATH_PROJECTED_CONTINUITY]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[PN_0_NEWTONIAN_HOOK]] | Claim feeds this downstream object, output, or open gate. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_1PN_BRIDGE_DICTIONARY]] | v0.6 section anchor for CLAIM_PROJECTION_OPEN_BRANE_SYSTEM |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_POISSON_HOOK]] | Poisson hook from longitudinal identity. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_PROJECTED_CONTINUITY]] | Exact leakage identity. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_4D_PROJECTION_REDUCTION]] | Projection vs reduction firewall. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_MT_PROJECTION_ZERO_MODE]] | Projection and zero-mode hooks. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PLASMA_LEAKAGE]] | Open-system brane identities. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_PLASMA_Z_VS_W]] | Z controls action; W controls observation. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_4D_PARENT]] | Paper backlink block references CLAIM_PROJECTION_OPEN_BRANE_SYSTEM. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_PLASMA]] | Paper backlink block references CLAIM_PROJECTION_OPEN_BRANE_SYSTEM. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_and_controlled_reduction |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_BRANE_OBSERVER]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_OPEN_CONDUIT]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_SUPERFLUID_INTAKE_OUTPUT]] | Physical ontology object grounded by this claim. |
| `LINKS_TO` | [[PHYS_REG_OBSERVER_PROJECTION]] | Physical register entry links to graph object. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_4D_PARENT]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_PLASMA]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_PROJECTION_NOT_REDUCTION]] | Firewall preserves this correct status boundary. |
| `SUPPLIES_PARENT_FOR` | [[CLAIM_PARENT_ACTION_CURRENT_EXACT]] | Claim-level dependency added in v0.4. |
| `SUPPORTS_CLAIM` | [[EQ_LONGITUDINAL_IDENTITY]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_PROJECTED_CONTINUITY]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_4D_PARENT]]
- [[FILE_PLASMA]]
- [[SEC_1PN_BRIDGE_DICTIONARY]]
- [[SEC_4D_POISSON_HOOK]]
- [[SEC_4D_PROJECTED_CONTINUITY]]
- [[SEC_4D_PROJECTION_REDUCTION]]
- [[SEC_MT_PROJECTION_ZERO_MODE]]
- [[SEC_PLASMA_LEAKAGE]]
- [[SEC_PLASMA_Z_VS_W]]

### Source files
- `research/4d_plasma/paper/4d_plasma.tex`

### TeX anchor
- File: `research/4d_plasma/paper/4d_plasma.tex`
- Line: `252`
- Heading: Localization and Projection Operators
- Nearest label: `sec:projection_defs` at line `253`
- Match basis: `graph_edge:ANCHORS_CLAIM_SECTION`
- Confidence: `medium`

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
