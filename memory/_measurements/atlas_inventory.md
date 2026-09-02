# Phase 0 Atlas Migration Inventory

Date: 2026-08-24
Superseded for current corpus selection: 2026-09-01 by
[`pde_ledger_exclusion_validation.md`](pde_ledger_exclusion_validation.md).
The historical audit below remains as recorded.
Repository snapshot: `30e96ee22245d4a7d0e873ad228cf5ce33de76f0`
Machine-readable inventory: [`../_meta/atlas-migration.yaml`](../_meta/atlas-migration.yaml)

This was a read-only audit of `atlas/` and `graph/`. No legacy file and no
physics status was changed. The inventory is a deletion gate: it identifies
what should be re-extracted into the new memory before the two legacy trees are
removed.

## Authority and age check

`atlas/AGENTS.md` says that generated Obsidian nodes must not be edited and that
the maintained graph YAML is the graph-level source. `atlas/README.md` likewise
identifies `graph/fluid_universe_derivation_atlas_graph.yaml` as the durable
Atlas content source, while Canvas, Bases, nodes, exports, and the reader site
are generated views.

The maintained graph and `atlas/topics.yaml` were last committed together on
2026-05-15 at `0ca49cd73f64452d45a8b3b47780fd5f7d1aa103`. The current tracked source
lineage is newer: for example, the v3 charter was changed on 2026-08-02 and the
v3 defect register on 2026-08-17. The v3 charter also records a requirements-
first method change, a self-contained v3 lineage, changed treatment of the
lepton tower, and explicit carry-forward decisions.

Consequence: Atlas statuses such as `active_v07`, `open`, and
`open_actual_branch_data` are historical metadata. The migration inventory
never treats them as current. Any migrated page must re-read the eligible
`research/` and `software/` sources and apply the memory precedence rules.

## Structural measurement

Commands:

```bash
git ls-files atlas | wc -l
git ls-files graph | wc -l
jq -r '.nodes | length' graph/fluid_universe_derivation_atlas_graph.json
jq -r '.nodes | group_by(.layer)[] | "\(.[0].layer)\t\(length)"' \
  graph/fluid_universe_derivation_atlas_graph.json
rg -c '^- topic_slug:' atlas/exports/future_paper_backfill.yaml
rg -c '^  - slug:' atlas/topics.yaml
```

Literal summary:

- 457 tracked paths under `atlas/` and 21 under `graph/`.
- 382 graph nodes and 1,348 graph edges.
- 382 generated Obsidian node notes.
- 12 hand-authored reader topic definitions.
- 92 generated future-paper backfill rows.
- Graph layers: 34 claim/theorem, 28 derivation, 25 equation, 27 math,
  23 physical ontology, 6 physical register, 16 open gate, 59 status audit,
  24 query validation, 20 file anchor, 76 section anchor, 16 paper backlink,
  15 Atlas meta, 8 Atlas completion phase, and 5 extension workflow nodes.

The Atlas validator's recorded final report has zero errors and zero warnings.
That establishes internal graph/view consistency, not current scientific
authority.

## What is genuinely worth migrating

The useful material is concentrated in a small fraction of the legacy system:

- Fourteen of the seventeen firewall-like records are re-anchorable and worth
  retaining. This includes the sixteen-node status register plus the separate
  mouth/interior ontology firewall. The atom, lepton, and g-2 firewalls are
  excluded because their exact support is the root `notes/` future-work corpus;
  the lepton record also conflicts with the newer v3 lineage if copied
  unchanged.
- Five of sixteen open-gate records are useful as explicit research questions:
  actual branch export, executable branch solving, material closure, mixed
  recirculation reconciliation, and the nonlinear/closed-parent-action handoff.
  Five more are already better covered by tracked sources, four are duplicates
  of those protocols or gates, and two are excluded notes-only topics.
- Eight reader topic shells are useful navigation seeds: foundations,
  projection/reduction, charge/EM, the PN ladder, quadrupole normalization,
  moving-throat dynamics, branch response/export, and status/reading rules.
  Their prose is not safe to copy; it must be rewritten from current sources.
- The response-readout ontology is uniquely useful: `D0`, `C`, `P0`, `N2`, and
  `N4` are compressed outputs, not the physical throat. The target-blind branch
  and exporter are upstream physical/protocol objects.
- The eight positive retrieval questions and sixteen negative/firewall queries
  are good regression-test seeds. The initial benchmark disables four that
  depend on excluded atom/lepton/g-2 notes. Expected answers must be regenerated
  rather than copied.

Everything else is either recoverable directly from canonical sources or is
Atlas machinery: graph mirrors, node notes, canvases, bases, reader output,
backlink blocks, completion dashboards, extension workflow nodes, and
source/section-node representations.

## Re-anchoring findings

The most successful re-anchoring target is `research/pde_audit/`. The old graph
often cites `notes/pde_audit_full.md:V2-xx`, but the repository contains one
tracked derivation note per audit stage, paired scripts, literal output, and a
simulation protocol. Examples include:

- Parent/wall status and no-post-hoc wall coefficients:
  `research/pde_audit/notes/stage_v2_01_parent_wall_action_derivation.md`.
- Gauge-localization caveat:
  `research/pde_audit/notes/stage_v2_02_maxwell_gauge_localization_derivation.md`.
- Actual branch, Packet A/B, no-refit, and exporter requirements:
  `research/pde_audit/notes/stage_v2_25_actual_branch_protocol_derivation.md`
  and `research/pde_audit/simulation/ACTUAL_BRANCH_PROTOCOL_V1.md`.
- Physical mouth/interior/readout ontology:
  `research/pde_audit/notes/stage_v2_28_physical_picture_and_ontology.md`.
- Material closure:
  `research/pde_audit/notes/stage_v2_29_superfluid_material_closure_gap.md`.
- Charge/circulation and mixed-plumbing status:
  `research/pde_audit/notes/stage_v2_30_electromagnetic_ontology_and_status.md`.

The paper-backed Atlas nodes also re-anchor cleanly to TeX labels in the 4D,
EM, plasma, and PN papers. The componentized PDE ledger has its own maintained
dependency map at `research/pde_ledger/paper/frontmatter/05_dependency_map.tex`,
so the old graph's per-stage and per-equation nodes are unnecessary.

The weak re-anchors are intentionally visible in the YAML. Atom, lepton, g-2,
and several circulation conclusions depend on root `notes/` files and have no
eligible source for the exact Atlas claim. They are marked `obsolete` with
`reanchor_complete: false`; they cannot silently enter the new semantic memory.

## Disposition totals for enumerated unique records

| Record family | Migrate | Covered by source | Generated duplicate | Obsolete |
|---|---:|---:|---:|---:|
| Firewalls, including mouth/interior | 14 | 0 | 0 | 3 |
| Open gates | 5 | 5 | 4 | 2 |
| Reader topics | 8 | 1 | 0 | 3 |

All 34 claim/theorem nodes are accounted for in four groups in the YAML:
twelve canonical-paper summaries, fourteen PDE-program summaries, one
paper-backed fluxoid identity, and seven excluded notes-only claims. Equations,
math objects, derivations, and ordinary physical-ontology nodes use a simpler
rule: eligible source-backed records are regenerated in source capsules or
topic synthesis; notes-only records are not migrated. The few unique ontology
guards are enumerated separately.

## Validation

The inventory was parsed with Ruby's safe YAML loader and checked for the
expected unique counts:

```text
firewalls=17
open_gates=16
topics=12
query_groups=2
```

Every path in an `original_sources` list resolves to an existing tracked
`research/` or `software/` path after removing its stable label/heading suffix.
The inventory preserves Atlas/graph only until its migration and benchmark
gates pass; it does not recommend keeping the legacy generators afterward.
