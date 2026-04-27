# Codex Prompt — Build Fluid Universe Obsidian Atlas v0.9

You are working in the Fluid Universe git repository. Implement the Obsidian-facing atlas layer using the existing Fluid Universe Derivation Atlas artifacts.

## Goal

Create a read-only-for-human, AI-maintained Obsidian vault layer inside the repo. The vault should make the derivation atlas easy to browse with Markdown summaries, MathJax equations, Dataview dashboards, Bases views, and curated Canvas maps.

The user will not manually edit the atlas. Codex/AI agents must be able to regenerate, validate, and synchronize it.

## Source artifacts

Start from the latest available atlas files, preferring:

- `fluid_universe_derivation_atlas_graph_v08.json`
- `fluid_universe_derivation_atlas_graph_v08.yaml`
- `fluid_universe_derivation_atlas_v08_status_dashboard.md`
- `fluid_universe_derivation_atlas_status_firewall_register_v07.md`
- `fluid_universe_derivation_atlas_physical_register_v04.md`
- `fluid_universe_derivation_atlas_claim_register_v04.md`
- `fluid_universe_derivation_atlas_open_gates_v04.md`
- `fluid_universe_derivation_atlas_paper_backlink_register_v06.md`
- `fluid_universe_derivation_atlas_paper_insertion_manifest_v08.json`
- local Codex backlink application report, if present.

Also use the project source/summary files as anchors, especially:

- `4d_summary.md`
- `4d_em_fields_summary.md`
- `4d_plasma_summary.md`
- `4d_1pn_full_summary.md`
- `4d_2pn_summary.md`
- `4d_2_5pn_summary.md`
- `4d_3pn_summary.md`
- `4d_4pn_summary.md`
- `4d_5pn_summary.md`
- `moving_throat_pde_program_compact.md`
- `moving_throat_output_full.md`
- `pde_audit_full.md`
- `atom_work.md`
- `lepton_work.md`
- `lepton_mass_notes.md`
- `g2_full_output.md`
- `5pn_notes_full.md`

## Required outputs

Create:

```text
atlas/README.md
atlas/AGENTS.md
atlas/nodes/**/<NODE_ID>.md
atlas/views/*.md
atlas/canvas/*.canvas
atlas/bases/*.base
atlas/templates/*.md
atlas/scripts/*.py
atlas/exports/atlas_graph.json
atlas/exports/atlas_graph.yaml
atlas/exports/validation_report.md
```

## Architecture

1. Convert every v0.8 graph node into a Markdown note.
2. Preserve every node ID exactly.
3. Store machine-queryable metadata in YAML frontmatter.
4. Store readable summaries, equations, physical interpretation, status explanation, and wikilinks in Markdown body.
5. Generate graph JSON/YAML exports from the Markdown notes after conversion.
6. Generate curated Canvas maps from node notes. Canvas files should use file cards, not text-only cards, except for labels.
7. Create Dataview dashboard notes and Bases `.base` files.
8. Add validation scripts.

## Required node body sections

Every node should include:

```markdown
# Title

> **Status:** ...
> **Layer:** ...
> **Type:** ...

## Summary
## Physical meaning
## Mathematical role
## Atlas links
## Source anchors
## AI maintenance notes
```

Equation nodes must include `## Equation` and `## Variable dictionary`.
Open gates must include `## What remains open` and `## What would close it`.
Status firewalls must include `## Invalid inference` and `## Corrected inference`.

## Important summary policy

The user wants embedded summaries so the vault is readable without constantly opening source files. Therefore, do not make notes into bare metadata stubs. Add useful body summaries.

Recommended body depth:

- physical ontology nodes: 2–5 paragraphs,
- claim nodes: 1–3 paragraphs,
- equation nodes: equation + variable dictionary + derivation role,
- open gates: what remains open + what would close it + affected claims,
- status firewalls: invalid shortcut + corrected reading + affected branches,
- source nodes: paper purpose + claim taxonomy + main carry-forward outputs.

Keep equations in Markdown body with `$$ ... $$`, not YAML.

## Status firewalls to preserve

Do not allow the conversion to weaken these firewalls:

1. Electric charge sign is `eta_Q`, not circulation.
2. `q_star` / `q_eff` are EM charge objects; `kappa_rho=1` is a gravity-side scalar mass-dressing coefficient.
3. Zero-mode Maxwell suppresses mixed channels only in the controlled far-field brane limit; it does not erase `A_w`, `J^w`, `F_{mu w}`, `E_w`, or `C_a` from the microscopic ontology.
4. Projection is not reduction.
5. Effective wall closure is not strict parent-level moving-throat theorem unless `S_eta` or `S_Sigma` is promoted.
6. Local instantaneous 4PN closure is not unconditional full 4PN tail completion.
7. The angular source-map identity does not solve radial/axial quadrupole normalization.
8. Monomial/similarity orbit structure does not by itself close the full 5PN program.
9. Atom/lepton/anomaly branches remain reduced or conditional where the source files say they are.

## First implementation phases

### Phase 1 — scaffold

Create the folder structure, README, AGENTS, templates, and `.gitignore` entries.

### Phase 2 — mechanical conversion

Convert all v0.8 graph nodes to Markdown notes. Use the graph fields as the first-pass summaries when no better source summary is available.

### Phase 3 — enrichment pass

Enrich the high-value backbone nodes first:

- parent action / projection / charge ontology,
- Maxwell reduction and mixed-sector firewall,
- physical finite throat ontology,
- moving-throat PDE stages,
- parent-wall action audit,
- 1PN–4PN chain,
- 2.5PN/4PN quadrupole normalization gate,
- major open gates,
- major status firewalls,
- atom/lepton/anomaly reduced branches,
- 5PN continuation.

### Phase 4 — dashboards and canvas

Create Dataview dashboards, Bases files, and focused Canvas maps.

### Phase 5 — validation and report

Run validation and write:

```text
atlas/exports/validation_report.md
```

Also write a human-readable implementation report:

```text
atlas/exports/obsidian_atlas_application_report.md
```

The report should include:

- number of notes created,
- number of edges represented,
- number of Canvas files created,
- Dataview/Bases files created,
- validation pass/fail summary,
- unresolved or ambiguous nodes,
- source files that could not be found,
- status changes made, if any.

## Do not do

- Do not invent new physics claims.
- Do not promote claim statuses unless the source files explicitly support it.
- Do not make Canvas the only source of graph truth.
- Do not put long prose or equations in YAML frontmatter.
- Do not erase mixed-sector objects from the ontology because of zero-mode reduction.
- Do not collapse projection and reduction.
- Do not rewrite stable atlas IDs for readability.
- Do not treat summary-file section positions as exact full-paper locations unless the full paper is present.

## Final expected state

After this pass, opening the repository as an Obsidian vault should show:

1. searchable atlas node notes,
2. readable summaries with equations,
3. dashboards for claims/open gates/status firewalls,
4. curated Canvas maps for major derivation paths,
5. regenerated graph exports for AI tooling,
6. a validation report proving the vault is internally consistent.
