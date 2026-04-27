---
id: CLAIM_5PN_FULL_BUNDLE_SURFACE
title: 5PN full-bundle surface and weak-axisymmetric continuation
type: claim
layer: claim_theorem
status: exact_within_reduced_bundle_open_branch
atlas_version: obsidian-v0.1
source_graph_version: v0.8-codex-handoff
source_graph_file: graph/fluid_universe_derivation_atlas_graph.yaml
generated_by: codex
generated: true
last_generated_utc: '2026-04-27T19:11:17Z'
summary_short: The 5PN/moving-throat continuation gives exact isotropic full-bundle target surfaces and weak-axisymmetric quotient/orbit structures, but the actual branch still must land on them.
future_paper_needed: false
source_links:
- '[[FILE_5PN_FULL]]'
- '[[FILE_MOVING_THROAT_COMPACT]]'
- '[[SEC_5PN_STAGE15_XI]]'
- '[[SEC_5PN_STAGE18_SURFACE]]'
- '[[SEC_5PN_STAGE20_EVEN_GATES]]'
- '[[SEC_5PN_SUMMARY_BOTTOM_LINE]]'
- '[[SEC_5PN_SUMMARY_CLAIMS]]'
- '[[SEC_5PN_SUMMARY_MONOMIALS]]'
- '[[SEC_5PN_SUMMARY_SCOPE]]'
physical_ids:
- PHYS_FINITE_MOUTH_SHAPE
- PHYS_RESPONSE_READOUTS
equation_ids:
- EQ_FULL_BUNDLE_TARGET_SURFACE
- EQ_MONOMIAL_QUOTIENT
- EQ_WEAK_AXISYM_SIGNATURE
claim_ids:
- CLAIM_G2_COMMON_QUOTIENT
- CLAIM_STAGE6_FULL_BUNDLE_RATIO
open_gate_ids:
- OPEN_N2_N4_MOMENT_SHAPE
- OPEN_WEAK_AXISYM_ORBIT_LOCK
status_firewall_ids:
- FIREWALL_SIMILARITY_NOT_FULL_5PN
source_ids:
- FILE_5PN_FULL
- FILE_MOVING_THROAT_COMPACT
- SEC_5PN_STAGE15_XI
- SEC_5PN_STAGE18_SURFACE
- SEC_5PN_STAGE20_EVEN_GATES
- SEC_5PN_SUMMARY_BOTTOM_LINE
- SEC_5PN_SUMMARY_CLAIMS
- SEC_5PN_SUMMARY_MONOMIALS
- SEC_5PN_SUMMARY_SCOPE
outgoing_edges:
- target: OPEN_N2_N4_MOMENT_SHAPE
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle_open_branch
  note: Claim feeds this downstream object, output, or open gate.
- target: OPEN_WEAK_AXISYM_ORBIT_LOCK
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle_open_branch
  note: Claim feeds this downstream object, output, or open gate.
- target: PN_5_FULL_BUNDLE_SURFACE
  relation: FEEDS_OR_STATUS_OF
  status: exact_within_reduced_bundle_open_branch
  note: Claim feeds this downstream object, output, or open gate.
- target: CLAIM_STAGE6_FULL_BUNDLE_RATIO
  relation: REFINES
  status: active
  note: Claim-level dependency added in v0.4.
- target: CLAIM_G2_COMMON_QUOTIENT
  relation: SUPPLIES_XI1_LANGUAGE_FOR
  status: active
  note: Claim-level dependency added in v0.4.
incoming_edges:
- source: SEC_5PN_STAGE15_XI
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Xi_load bridge to prefactor slope.
- source: SEC_5PN_STAGE18_SURFACE
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Exact isotropic full-bundle target surface.
- source: SEC_5PN_STAGE20_EVEN_GATES
  relation: ANCHORS_CLAIM_SECTION
  status: v05
  note: Full even-gate solve with reinstated Z-sector.
- source: SEC_5PN_SUMMARY_BOTTOM_LINE
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE
- source: SEC_5PN_SUMMARY_CLAIMS
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE
- source: SEC_5PN_SUMMARY_MONOMIALS
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE
- source: SEC_5PN_SUMMARY_SCOPE
  relation: ANCHORS_CLAIM_SECTION
  status: v06
  note: v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE
- source: BACKLINK_5PN_FULL
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_5PN_FULL_BUNDLE_SURFACE.
- source: BACKLINK_G2_OUTPUT
  relation: BACKLINKS_ATLAS_NODE
  status: v06
  note: Paper backlink block references CLAIM_5PN_FULL_BUNDLE_SURFACE.
- source: STATUS_LADDER_EXACT_TO_OPEN
  relation: CLASSIFIES
  status: exact_within_reduced_bundle_open_branch
  note: 'Claim class: exact_within_closure_plus_open_gate'
- source: VIEW_CLAIM_LAYER
  relation: CONTAINS_CLAIM
  status: active
  note: v0.4 claim/theorem layer item.
- source: PHYS_FINITE_MOUTH_SHAPE
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_bundle_open_branch
  note: Physical ontology object grounded by this claim.
- source: PHYS_RESPONSE_READOUTS
  relation: GROUNDS_PHYSICAL_MEANING
  status: exact_within_reduced_bundle_open_branch
  note: Physical ontology object grounded by this claim.
- source: FILE_5PN_FULL
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Source artifact anchors this claim.
- source: FILE_MOVING_THROAT_COMPACT
  relation: OWNS_OR_ANCHORS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Source artifact anchors this claim.
- source: FIREWALL_SIMILARITY_NOT_FULL_5PN
  relation: PROTECTS_STATUS_OF
  status: v07
  note: Firewall preserves this correct status boundary.
- source: EQ_FULL_BUNDLE_TARGET_SURFACE
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Equation anchor supports this named claim.
- source: EQ_MONOMIAL_QUOTIENT
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Equation anchor supports this named claim.
- source: EQ_WEAK_AXISYM_SIGNATURE
  relation: SUPPORTS_CLAIM
  status: exact_within_reduced_bundle_open_branch
  note: Equation anchor supports this named claim.
tags:
- atlas/claims
- atlas/node
- layer/claim_theorem
- status/exact_within_reduced_bundle_open_branch
- topic/pn_chain
- type/claim
---

<!-- GENERATED FROM graph/fluid_universe_derivation_atlas_graph.yaml; DO NOT EDIT BY HAND. -->

# 5PN full-bundle surface and weak-axisymmetric continuation

> [!warning] Generated note
> This Obsidian note is generated from `graph/fluid_universe_derivation_atlas_graph.yaml`. Do not edit by hand; update the graph and regenerate.

> **Atlas ID:** `CLAIM_5PN_FULL_BUNDLE_SURFACE`  
> **Status:** `exact_within_reduced_bundle_open_branch`  
> **Layer:** `claim_theorem`  
> **Type:** `claim`

## Summary

The 5PN/moving-throat continuation gives exact isotropic full-bundle target surfaces and weak-axisymmetric quotient/orbit structures, but the actual branch still must land on them.

## Claim

The 5PN/moving-throat continuation gives exact isotropic full-bundle target surfaces and weak-axisymmetric quotient/orbit structures, but the actual branch still must land on them.

## What It Does Not Claim

This generated note preserves the graph status. It should not be read as closing an open gate or promoting a conditional theorem.

## Physical Meaning

The 5PN/moving-throat continuation gives exact isotropic full-bundle target surfaces and weak-axisymmetric quotient/orbit structures, but the actual branch still must land on them.

## Mathematical Role

- Layer: `claim_theorem`
- Type: `claim`
- Status: `exact_within_reduced_bundle_open_branch`
- Outputs: `PN_5_FULL_BUNDLE_SURFACE`, `OPEN_N2_N4_MOMENT_SHAPE`, `OPEN_WEAK_AXISYM_ORBIT_LOCK`

## Atlas Links

### Related physical nodes
- [[PHYS_FINITE_MOUTH_SHAPE]]
- [[PHYS_RESPONSE_READOUTS]]

### Related math nodes
- none

### Related equations
- [[EQ_FULL_BUNDLE_TARGET_SURFACE]]
- [[EQ_MONOMIAL_QUOTIENT]]
- [[EQ_WEAK_AXISYM_SIGNATURE]]

### Related claims
- [[CLAIM_G2_COMMON_QUOTIENT]]
- [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]]

### Open gates
- [[OPEN_N2_N4_MOMENT_SHAPE]]
- [[OPEN_WEAK_AXISYM_ORBIT_LOCK]]

### Status firewalls
- [[FIREWALL_SIMILARITY_NOT_FULL_5PN]]

### Source anchors
- [[FILE_5PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[SEC_5PN_STAGE15_XI]]
- [[SEC_5PN_STAGE18_SURFACE]]
- [[SEC_5PN_STAGE20_EVEN_GATES]]
- [[SEC_5PN_SUMMARY_BOTTOM_LINE]]
- [[SEC_5PN_SUMMARY_CLAIMS]]
- [[SEC_5PN_SUMMARY_MONOMIALS]]
- [[SEC_5PN_SUMMARY_SCOPE]]

## Outgoing Edges

| Relation | Node | Note |
|---|---|---|
| `FEEDS_OR_STATUS_OF` | [[OPEN_N2_N4_MOMENT_SHAPE]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[OPEN_WEAK_AXISYM_ORBIT_LOCK]] | Claim feeds this downstream object, output, or open gate. |
| `FEEDS_OR_STATUS_OF` | [[PN_5_FULL_BUNDLE_SURFACE]] | Claim feeds this downstream object, output, or open gate. |
| `REFINES` | [[CLAIM_STAGE6_FULL_BUNDLE_RATIO]] | Claim-level dependency added in v0.4. |
| `SUPPLIES_XI1_LANGUAGE_FOR` | [[CLAIM_G2_COMMON_QUOTIENT]] | Claim-level dependency added in v0.4. |

## Incoming Edges

| Relation | Node | Note |
|---|---|---|
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_STAGE15_XI]] | Xi_load bridge to prefactor slope. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_STAGE18_SURFACE]] | Exact isotropic full-bundle target surface. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_STAGE20_EVEN_GATES]] | Full even-gate solve with reinstated Z-sector. |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_BOTTOM_LINE]] | v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_CLAIMS]] | v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_MONOMIALS]] | v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE |
| `ANCHORS_CLAIM_SECTION` | [[SEC_5PN_SUMMARY_SCOPE]] | v0.6 section anchor for CLAIM_5PN_FULL_BUNDLE_SURFACE |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_5PN_FULL]] | Paper backlink block references CLAIM_5PN_FULL_BUNDLE_SURFACE. |
| `BACKLINKS_ATLAS_NODE` | [[BACKLINK_G2_OUTPUT]] | Paper backlink block references CLAIM_5PN_FULL_BUNDLE_SURFACE. |
| `CLASSIFIES` | [[STATUS_LADDER_EXACT_TO_OPEN]] | Claim class: exact_within_closure_plus_open_gate |
| `CONTAINS_CLAIM` | [[VIEW_CLAIM_LAYER]] | v0.4 claim/theorem layer item. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_FINITE_MOUTH_SHAPE]] | Physical ontology object grounded by this claim. |
| `GROUNDS_PHYSICAL_MEANING` | [[PHYS_RESPONSE_READOUTS]] | Physical ontology object grounded by this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_5PN_FULL]] | Source artifact anchors this claim. |
| `OWNS_OR_ANCHORS_CLAIM` | [[FILE_MOVING_THROAT_COMPACT]] | Source artifact anchors this claim. |
| `PROTECTS_STATUS_OF` | [[FIREWALL_SIMILARITY_NOT_FULL_5PN]] | Firewall preserves this correct status boundary. |
| `SUPPORTS_CLAIM` | [[EQ_FULL_BUNDLE_TARGET_SURFACE]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_MONOMIAL_QUOTIENT]] | Equation anchor supports this named claim. |
| `SUPPORTS_CLAIM` | [[EQ_WEAK_AXISYM_SIGNATURE]] | Equation anchor supports this named claim. |

## Source Anchors

### Source anchor notes
- [[FILE_5PN_FULL]]
- [[FILE_MOVING_THROAT_COMPACT]]
- [[SEC_5PN_STAGE15_XI]]
- [[SEC_5PN_STAGE18_SURFACE]]
- [[SEC_5PN_STAGE20_EVEN_GATES]]
- [[SEC_5PN_SUMMARY_BOTTOM_LINE]]
- [[SEC_5PN_SUMMARY_CLAIMS]]
- [[SEC_5PN_SUMMARY_MONOMIALS]]
- [[SEC_5PN_SUMMARY_SCOPE]]

### Source files
- No source file path recorded.

## AI Maintenance Notes

- Treat this file as generated read-only presentation material.
- Change graph YAML, TeX papers, or source notes first, then regenerate the Obsidian layer.
- Do not use this note to upgrade statuses or weaken firewalls.
