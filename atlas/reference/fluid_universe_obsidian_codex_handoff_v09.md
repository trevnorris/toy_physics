# Fluid Universe Atlas → Obsidian Vault Handoff v0.9

## Executive intent

Build an Obsidian-facing atlas layer for the Fluid Universe derivation graph.

This is **not** a manual note-taking system. The user will treat the vault as read-only. Codex / AI agents are responsible for generating, updating, validating, and synchronizing the atlas data.

The Obsidian layer should make the atlas easier to view, browse, and reason over by combining:

1. structured YAML properties for Dataview / Bases / validation,
2. readable Markdown summaries with MathJax-compatible equations,
3. internal wikilinks for Obsidian graph and local navigation,
4. curated Canvas maps for visual traversal,
5. JSON/YAML graph exports for AI and script workflows,
6. Git version control for all durable artifacts.

The key design rule is:

> Markdown node notes are the human-facing canonical atlas. JSON/YAML graph exports and Canvas files are generated or synchronized products. Canvas is a visualization layer, not the source of truth.

---

## Current atlas inputs

Use the latest existing atlas artifacts already present in the repository or workspace:

- `fluid_universe_derivation_atlas_graph_v08.json`
- `fluid_universe_derivation_atlas_graph_v08.yaml`
- `fluid_universe_derivation_atlas_v08_status_dashboard.md`
- `fluid_universe_derivation_atlas_status_firewall_register_v07.md`
- `fluid_universe_derivation_atlas_negative_query_report_v07.md`
- `fluid_universe_derivation_atlas_query_validation_report_v06.md`
- `fluid_universe_derivation_atlas_physical_register_v04.md`
- `fluid_universe_derivation_atlas_claim_register_v04.md`
- `fluid_universe_derivation_atlas_open_gates_v04.md`
- `fluid_universe_derivation_atlas_paper_backlink_register_v06.md`
- `fluid_universe_derivation_atlas_paper_insertion_manifest_v08.json`
- Codex paper-linking application report, if present locally.

The v0.8 graph currently has 373 nodes and 1325 edges. Preserve stable node IDs exactly.

---

## Repository layout to create

Create this structure under the existing git repository:

```text
atlas/
  README.md
  AGENTS.md

  nodes/
    physical/
    math/
    equations/
    claims/
    open_gates/
    status_firewalls/
    sources/
    derivations/
    audits/
    extensions/
    meta/

  views/
    00_home.md
    01_open_gates.md
    02_claim_register.md
    03_physics_to_math.md
    04_moving_throat_pde.md
    05_pn_chain.md
    06_charge_ontology_firewall.md
    07_status_firewalls.md
    08_coverage_audit.md
    09_source_crosswalk.md
    10_query_validation.md

  canvas/
    00_backbone.canvas
    01_physical_ontology.canvas
    02_moving_throat_pde.canvas
    03_quadrupole_normalization.canvas
    04_pn_chain.canvas
    05_open_gates.canvas
    06_charge_ontology_firewall.canvas
    07_status_firewalls.canvas
    08_lepton_atom_anomaly.canvas

  bases/
    claims.base
    open_gates.base
    physical_objects.base
    status_firewalls.base
    source_crosswalk.base

  exports/
    atlas_graph.json
    atlas_graph.yaml
    atlas_edges.tsv
    atlas_nodes.tsv
    mermaid_backbone.mmd
    mermaid_open_gates.mmd
    validation_report.md

  templates/
    node_template.md
    claim_node_template.md
    physical_node_template.md
    equation_node_template.md
    open_gate_template.md
    status_firewall_template.md
    source_anchor_template.md

  scripts/
    generate_obsidian_nodes.py
    export_graph_from_nodes.py
    generate_canvas_views.py
    validate_atlas.py
    sync_backlinks.py

  .obsidian-suggested/
    README.md
    app.json
    community-plugins.json
```

Do not commit volatile Obsidian workspace state. Add or update `.gitignore` with:

```gitignore
.obsidian/workspace.json
.obsidian/workspace-mobile.json
.obsidian/cache/
.obsidian/plugins/*/data.json
```

It is acceptable to commit a minimal `.obsidian-suggested/` folder, but the atlas must work as plain Markdown even without Obsidian installed.

---

## Canonical note strategy

Each atlas node becomes one Markdown file.

Filename format:

```text
atlas/nodes/<folder>/<NODE_ID>.md
```

Example:

```text
atlas/nodes/claims/CLAIM_PARENT_WALL_STATUS_SPLIT.md
atlas/nodes/physical/PHYS_FINITE_BRANE_BULK_THROAT.md
atlas/nodes/open_gates/OPEN_ACTUAL_BRANCH_EXPORTER.md
```

Do not rename IDs to prettier titles. Stable IDs are critical for AI continuity.

Each note has two layers:

1. **YAML frontmatter**: short, structured, machine-queryable fields.
2. **Markdown body**: readable summaries, MathJax equations, physical interpretation, status explanation, links, and source anchors.

Do not put long prose or display equations in YAML frontmatter. Put those in the body.

---

## Required YAML schema

Use this schema for all nodes. Some fields may be empty arrays, but the keys should exist.

```yaml
---
id: CLAIM_PARENT_WALL_STATUS_SPLIT
title: Parent wall status split
type: claim
layer: status_audit
status: patch_required
claim_strength: audit_patch
atlas_version: v0.9
source_graph_version: v0.8
generated_by: codex
last_generated_utc: 2026-04-25T00:00:00Z
summary_short: Confinement-only parent action supplies a wall force but not an autonomous moving-wall PDE.
source_files:
  - pde_audit_full.md
  - moving_throat_pde_program_compact.md
source_sections:
  - stage_v2_01_parent_wall_action_derivation.md
physical_ids:
  - PHYS_FINITE_BRANE_BULK_THROAT
math_ids:
  - MATH_SIGMA_R_FIELD
  - MATH_S_ETA_QUADRATIC
equation_ids:
  - EQ_CONFINEMENT_VARIATION
  - EQ_WALL_MODAL_PDE
claim_ids: []
open_gate_ids:
  - OPEN_PARENT_PROMOTION_S_SIGMA
status_firewall_ids:
  - FIREWALL_EFFECTIVE_WALL_NOT_PARENT_THEOREM
outgoing_edges:
  - target: OPEN_PARENT_PROMOTION_S_SIGMA
    relation: BLOCKED_BY
  - target: CLAIM_EFFECTIVE_LINEAR_WALL_CLOSURE
    relation: DISTINGUISHES_FROM
incoming_edges: []
tags:
  - atlas/claim
  - status/patch_required
  - layer/status_audit
  - topic/moving_throat
---
```

Field rules:

- `id` must equal the filename stem.
- `status` must belong to the existing atlas status enum or a documented extension.
- `source_files` are repository filenames, not web URLs.
- `*_ids` store stable IDs as strings.
- Actual Obsidian wikilinks should appear in the Markdown body, not only YAML.
- `summary_short` should be one sentence, no equations unless unavoidable.

---

## Required Markdown body structure

Each node note should use this body skeleton:

```markdown
# Parent wall status split

> **Status:** patch_required  
> **Layer:** status_audit  
> **Type:** claim

## Summary

Confinement-only geometry enters the current parent action through
$V_{\rm conf}(\mathbf X;\Sigma)$ and therefore supplies a wall force/source,
but it does not by itself supply an autonomous moving-wall PDE.

## Physical meaning

This node protects the distinction between a finite throat as a physical object
and the currently declared parent action as a confinement-coupled GNLS/Maxwell
system. The throat surface can be used as a coupling argument before it is a
fully dynamical parent field.

## Mathematical role

The confinement variation gives a source term proportional to

$$
\frac{\delta \mathcal L_\psi}{\delta\eta}
= \rho\,\frac{V_{\rm wall}'(\Sigma_0/\ell_c)}{\ell_c},
$$

but it does not produce terms of the form

$$
\eta_{tt},\qquad \partial_w(T_w\eta_w),\qquad \Delta_{S^2}\eta.
$$

Those arise only after adding a throat/wall action such as $S_\eta^{(2)}$ or a
nonlinear parent action $S_\Sigma$.

## Atlas role

- Blocks: [[OPEN_PARENT_PROMOTION_S_SIGMA]]
- Distinguishes from: [[CLAIM_EFFECTIVE_LINEAR_WALL_CLOSURE]]
- Related physical object: [[PHYS_FINITE_BRANE_BULK_THROAT]]
- Related math object: [[MATH_SIGMA_R_FIELD]]

## Status firewall

Do not upgrade the effective wall PDE into a strict parent-level theorem unless
$S_\eta$ or $S_\Sigma$ is promoted into $S_{\rm total}$.

## Source anchors

- `pde_audit_full.md`
- `moving_throat_pde_program_compact.md`

## AI maintenance notes

If the parent action is updated to include $S_\Sigma$, revise this node, the
parent-action claim nodes, and all downstream moving-throat status gates.
```

### Summary section policy

Because the vault is read-only for the user, every important node should include enough body text to be useful without opening the original paper or summary file. Use embedded summaries liberally, but keep them scoped:

- physical nodes: 2–5 paragraphs of physical interpretation,
- claim nodes: 1–3 paragraphs explaining what is claimed and what is not,
- equation nodes: equation, variable dictionary, derivation role,
- open gates: why open, what would close it, downstream consequences,
- status firewalls: invalid inference, corrected inference, affected nodes,
- source nodes: paper/summary purpose, claim taxonomy, linked atlas anchors.

This vault should feel like an atlas, not a spreadsheet.

---

## MathJax / equation policy

Use Obsidian-compatible Markdown math in the body:

```latex
$$
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}}
$$
```

Rules:

1. Put display equations in the body, not YAML.
2. Keep a blank line after closing `$$` before continuing prose.
3. Store equation IDs in YAML as strings, not equation text.
4. Use equation notes for reusable equations.
5. In equation notes, include a short variable dictionary.

Example equation note:

```markdown
# Zero-mode charge normalization

## Equation

$$
q_{\rm eff}=\frac{q_\star}{\sqrt{Z_{\rm int}}},
\qquad
q_\star=\eta_Q e_\star.
$$

## Meaning

The microscopic branch charge $q_\star$ is fixed by the topological electric
orientation $\eta_Q$, while the observable brane coupling is weakened by the
localization thickness $Z_{\rm int}$.
```

---

## Dataview dashboards

Create dashboard notes in `atlas/views/`. Keep queries scoped to the `atlas/` folder.

### Open gates

```dataview
TABLE status, summary_short, source_files
FROM "atlas/nodes/open_gates"
SORT status ASC, file.name ASC
```

### Claims by status

```dataview
TABLE layer, status, claim_strength, summary_short, source_files
FROM "atlas/nodes/claims"
SORT status ASC, layer ASC, file.name ASC
```

### Moving-throat PDE dashboard

```dataview
TABLE type, layer, status, summary_short, open_gate_ids
FROM "atlas/nodes"
WHERE contains(tags, "topic/moving_throat")
SORT layer ASC, status ASC, file.name ASC
```

### Status firewalls

```dataview
TABLE status, summary_short, source_files
FROM "atlas/nodes/status_firewalls"
SORT file.name ASC
```

### Coverage audit

```dataview
TABLE type, layer, status, physical_ids, math_ids, equation_ids, source_files
FROM "atlas/nodes"
WHERE type = "claim" AND (length(physical_ids) = 0 OR length(source_files) = 0)
SORT file.name ASC
```

### Nodes with no outgoing edges

```dataview
TABLE type, layer, status, summary_short
FROM "atlas/nodes"
WHERE length(outgoing_edges) = 0 AND type != "source"
SORT layer ASC, file.name ASC
```

---

## Bases dashboards

Also create `.base` files for users who prefer Obsidian Bases.

Example `atlas/bases/open_gates.base`:

```yaml
filters:
  and:
    - file.inFolder("atlas/nodes/open_gates")
views:
  - type: table
    name: Open gates
    order:
      - file.name
      - status
      - summary_short
      - source_files
```

Example `atlas/bases/claims.base`:

```yaml
filters:
  and:
    - file.inFolder("atlas/nodes/claims")
views:
  - type: table
    name: Claims
    order:
      - file.name
      - layer
      - status
      - claim_strength
      - summary_short
      - source_files
```

---

## Canvas policy

Canvas files should be curated visual maps generated from note files.

Do not create one enormous canvas containing every node.

Create focused canvases:

1. `00_backbone.canvas`: parent action → projection/reduction → PN chain → moving-throat bottleneck.
2. `01_physical_ontology.canvas`: bulk medium → throat defect → brane mouth → interior support → outgoing port.
3. `02_moving_throat_pde.canvas`: Stage 1 geometry lift → BdG → Maxwell/mixed → grouped P2 → selected branch.
4. `03_quadrupole_normalization.canvas`: 2.5PN → grouped P2 → outgoing fingerprint → P0 target → 4PN tail.
5. `04_pn_chain.canvas`: 1PN → 2PN → 2.5PN → 3PN → 4PN → 5PN continuation.
6. `05_open_gates.canvas`: live unresolved theorem gates and what would close them.
7. `06_charge_ontology_firewall.canvas`: eta_Q / q_star / q_eff / circulation / kappa_rho separation.
8. `07_status_firewalls.canvas`: invalid inference → corrected inference maps.
9. `08_lepton_atom_anomaly.canvas`: reduced atomic/lepton/anomaly branches and their status gates.

Use file cards pointing to real node notes. Avoid text-only cards except for headings or grouping labels, because file cards preserve backlinks and graph integrity.

Canvas generation should use a deterministic layout, even if crude. Manual visual refinements are acceptable only if preserved as layout metadata and not used to encode graph truth.

---

## Required validation script behavior

Create `atlas/scripts/validate_atlas.py` with these checks:

1. Every Markdown node filename stem equals its YAML `id`.
2. Every `id` is unique.
3. Every edge target/source refers to an existing node or an allowed external source anchor.
4. Every node has `title`, `type`, `layer`, `status`, `summary_short`, and `source_graph_version`.
5. Every claim node has at least one `source_files` entry.
6. Every open gate has at least one incoming or outgoing edge.
7. Every status firewall has an invalid/corrected inference section in the body.
8. Every equation node has an `## Equation` section.
9. No body contains obsolete charge language such as `q ~ a^2 Gamma` unless explicitly marked as deprecated.
10. No node upgrades a conditional/open status without a corresponding source/status change.
11. All Dataview and Bases files reference valid folders.
12. Canvas file cards point to existing Markdown files.
13. Generated graph export node/edge counts match the Markdown source unless explicitly documented.

The script should write:

```text
atlas/exports/validation_report.md
```

and exit nonzero on hard failures.

---

## AI update protocol

Because the user reads the vault but does not manually maintain it, Codex should follow this update sequence whenever a new paper, summary, or derivation file changes:

1. Read changed source files.
2. Identify affected atlas nodes.
3. Update YAML metadata.
4. Update body summaries and equations.
5. Add or remove node links.
6. Update open gates and status firewalls.
7. Regenerate graph exports.
8. Regenerate or patch Canvas views if the visual story changed.
9. Run validation.
10. Commit changes with a message like:

```text
atlas: update moving-throat PDE gate after parent-action audit
```

Status changes require special care. Never promote:

- effective wall closure → strict parent-level theorem,
- zero-mode Maxwell → mixed-sector erasure,
- local 4PN closure → unconditional full 4PN tail theorem,
- angular source-map closure → radial/axial normalization solved,
- monomial similarity orbit → full 5PN closure,
- reduced atom/lepton/anomaly result → full solved moving-throat theorem.

---

## Suggested initial conversion priority

Do not wait for perfect conversion of all 373 nodes before making the vault useful.

Phase 1 should convert all nodes from v0.8 graph mechanically.

Phase 2 should enrich the first 80–120 high-value nodes by adding readable summaries:

1. parent action / projection / charge ontology,
2. Maxwell zero-mode and mixed-sector firewall,
3. physical finite throat ontology,
4. moving-throat PDE stages,
5. 1PN / 2PN / 2.5PN / 3PN / 4PN chain,
6. quadrupole normalization open gate,
7. parent-wall action audit,
8. open gates and status firewalls,
9. atom / lepton / anomaly reduced branches,
10. 5PN continuation nodes.

Phase 3 should generate the curated Canvas maps.

---

## Completion criteria for v1.0 Obsidian atlas

Mark the Obsidian atlas as v1.0 when all are true:

- all v0.8 graph nodes have Markdown notes,
- all v0.8 graph edges are represented in YAML and body wikilinks,
- all major open gates have summary notes,
- all status firewalls have invalid/corrected inference notes,
- all major physical ontology objects have readable summaries,
- the core Dataview dashboards render without errors,
- Bases files exist for claims/open gates/status firewalls,
- curated Canvas maps exist and contain only valid file cards,
- validation passes,
- graph export regenerates from Markdown notes,
- README explains the read-only-for-human / AI-maintained workflow.
