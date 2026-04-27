# Fluid Universe Derivation Atlas v0.8 — Bidirectional Linking Policy

## Purpose

This policy defines how the atlas and the actual paper drafts should point to one another once Codex has access to the full repository.

## Canonical directionality

The atlas remains the canonical source for graph IDs, claim IDs, equation IDs, open-gate IDs, and firewall IDs.

The papers remain the canonical source for the prose derivations, equation numbering, and final presentation language.

The graph YAML may carry `tex_anchor` hints that point to paper-local lines,
headings, and nearby `\label{...}` entries. These are navigation aids, not new
claims, and should be regenerated if TeX headings or labels move.

The backlink blocks are deliberately small. They should not rewrite derivations; they should make navigation and status checking possible.

## Link contract

Each paper-side block should contain:

```text
Atlas version
Backlink block ID
Purpose
Primary atlas anchors
Source-section anchors
Open gates to preserve, if any
Status note
```

A paper section may have multiple atlas anchors, but every reusable theorem or claim should eventually map to at least one of:

```text
CLAIM_* node
EQ_* node
PHYS_* node
MATH_* node
OPEN_* node, if conditional/open
```

## Insertion rules

1. Insert blocks near theorem summaries, claim taxonomies, status/caveat sections, or transition sections.
2. Avoid placing blocks in the middle of a mathematical derivation unless the paper already has navigation sidebars or boxed metadata.
3. Do not duplicate existing blocks.
4. If a paper has multiple maintained versions, patch only the canonical maintained draft and report the rest.
5. Replace summary-line anchors with actual section anchors if the full paper has different section numbering; prefer existing `tex_anchor` metadata when it is present and semantically correct.
6. If a full paper has a stronger caveat than the atlas block, preserve the stronger caveat.

## Status discipline

Backlinks are navigational. They must not imply upgraded proof status.

Forbidden upgrades include:

- controlled zero-mode Maxwell -> exact full microscopic Maxwell erasure of mixed fields
- effective wall closure -> strict parent-level throat PDE
- angular STF identity -> solved radial/axial normalization
- local 4PN closure -> unconditional full 4PN tail theorem
- hydrogen reduced sector -> full atomic PDE theorem
- lepton conditional holonomy -> completed lepton theorem
- similarity orbit / monomial quotient -> full 5PN closure

## Versioning

When Codex applies blocks to full paper drafts, it may either:

- keep `Atlas version: v0.6` inside the original block to indicate the source block version, or
- update to `Atlas version: v0.8` and add `Backlink source: v0.6 block register`.

The preferred form is:

```text
Atlas version: v0.8
Backlink source: v0.6 block register
```

## Validation

After insertion, run a repository search for:

```text
Backlink block:
CLAIM_
OPEN_
FIREWALL_
```

The application report should show whether each manifest entry is patched, already present, missing, or ambiguous.
