# Phase 0 Corpus Boundary and Structural Inventory

Status: initial Phase 0 boundary
Snapshot inspected: `30e96ee22245d4a7d0e873ad228cf5ce33de76f0`
Measurement record: [`../_measurements/phase0_corpus_scan.md`](../_measurements/phase0_corpus_scan.md)

## Confirmed boundary

Normal memory synchronization considers only Git-tracked content under:

- `research/`
- `software/`

The root `notes/` and `docs/` trees are excluded. They contain future ideas or
secondary writeups rather than the material this memory should catalog.

`atlas/` and `graph/` are not normal corpus roots. They are temporary legacy
migration inputs, and every useful item taken from them must be re-anchored to
an original source under `research/` or `software/` before deletion.

## Git and ignore policy

The committed Git tree is the discovery boundary. Initial ingestion and normal
updates enumerate tree/blob objects from the selected commit under
`research/` and `software/`; they do not read the mutable index or recursively
walk the filesystem.

- Dirty tracked inputs are reported by `memory status` but are not processed
  into the committed memory baseline.
- Untracked files are not eligible until they are deliberately committed.
- Git-ignored files are never proposed as new sources.
- Deleted tracked sources remain visible through synchronization state so their
  capsules and dependencies can be retired correctly.

The root [`.gitignore`](../../.gitignore) is semantic policy, not merely cleanup:

- Run trees, `_scratch/`, caches, solver data, and large numeric formats are
  intentionally non-source artifacts.
- Small curated Markdown/YAML reports may be tracked and can be evidence.
- `.out` is deliberately **not** blanket-ignored because committed Wolfram and
  arbiter transcripts may be dual-engine evidence.
- The v3 ledger PDF is deliberately untracked while that paper is under
  construction.

Consequently, neither `output/` nor `.out` can be rejected merely by name.
Tracked evidence is included only as a dependency of a selected report, step,
paper, or script family; it does not automatically receive its own capsule.

## What the scan found

The candidate roots contain 7,839 tracked paths: 7,235 under `research/` and
604 under `software/`. The raw filesystem is much larger—especially
`software/`, where ignored `_scratch/` content accounts for 120,979 files.

The tracked research set contains:

- 392 TeX files.
- 2,774 Markdown files.
- 518 Python files.
- 423 Wolfram Language files.
- 20 TeX document entry points.

Most of the apparent research volume is process history, not primary memory
content: 3,752 tracked red-team paths and 863 tracked directive paths. These
remain available for provenance disputes but are excluded from normal semantic
extraction.

## Research structure

### Monolithic paper families

Sixteen paper entry points are self-contained TeX documents with no `\input`
tree. Each is initially one source unit and gets one source capsule.

| Family | Paper entry point |
|---|---|
| 1PN hybrid | `research/1pn_hybrid/paper/1pn_hybrid.tex` |
| 1PN optics | `research/1pn_optics/paper/1pn_optics.tex` |
| 1PN orbital dynamics | `research/1pn_orbital_dynamics/paper/1pn_orbital_dynamics.tex` |
| 1PN spin and N-body | `research/1pn_spin_and_nbody/paper/1pn_spin_and_nbody.tex` |
| Unified 4D parent | `research/4d/paper/4d.tex` |
| 4D 1PN bridge | `research/4d_1pn_bridge/paper/4d_1pn_bridge.tex` |
| 4D full 1PN | `research/4d_1pn_full/paper/4d_1pn_full.tex` |
| 4D 2.5PN | `research/4d_2_5pn/paper/4d_2_5pn.tex` |
| 4D 2PN | `research/4d_2pn/paper/4d_2pn.tex` |
| 4D 3PN | `research/4d_3pn/paper/4d_3pn.tex` |
| 4D 4PN | `research/4d_4pn/paper/4d_4pn.tex` |
| 4D EM fields | `research/4d_em_fields/paper/4d_em_fields.tex` |
| 4D plasma | `research/4d_plasma/paper/4d_plasma.tex` |
| Brane-bulk ontology | `research/brane_bulk_ontology/paper/brane_bulk_ontology.tex` |
| EM fields | `research/em_fields/paper/em_fields.tex` |
| Moving-throat PDE | `research/pde/paper/pde.tex` |

Many of these families also have a small `scripts/` and/or `mathematica/`
verification archive. The memory will catalog the referee/master harness or
declared suite for each paper, not every theorem-step helper. Existing README
and manifest files take precedence when identifying those entry points.

### Componentized PDE lineages

The PDE ledgers cannot be treated as one capsule per file.

#### `research/pde_ledger/`

This is a large archive ledger with two paper entry points, shared Parts I–VIII,
and 253 stage cards. Its paper README declares the stage cards to be canonical
per-stage provenance/audit anchors.

Initial treatment: supporting historical lineage and Atlas-migration source.
Do not generate 253 stage capsules. Read a stage file on demand when a selected
current source or migrated firewall cites it.

#### `research/pde_ledger_v2/`

This is an intermediate rebuild lineage that remains directly cited by current
software reports and v3 provenance. Its paper README and current tree are not
sufficient to declare the whole lineage obsolete.

Initial treatment: supporting lineage. Create capsules only for records that
are dependencies of active v3 steps, current software results, or selected
Atlas migrations.

#### `research/pde_ledger_v3/`

This is the active, self-contained requirements-first rebuild. It has:

- Governance and status documents (`CHARTER.md`, `DEFECT_REGISTER.md`,
  `SUBSTRATE_REQUIREMENTS.md`, and `V3_STEP_PLAN.md`).
- A componentized TeX paper.
- Completed/in-progress step-record families under `steps/`.
- Paper-facing step cards under `paper/steps/`.
- Paired SymPy and Wolfram engines plus comparators and export scripts.
- Explicit `_measurements/` records.

Initial treatment:

- One governance/status capsule.
- One paper capsule whose hash covers the componentized paper bundle.
- One capsule and one script-family entry per completed step ID, grouping its
  record, paper card, engines, comparators, and selected evidence.
- Directives, review prompts, repair history, and `_scratch/` remain
  citation-only unless a dispute requires them.

This grouping is essential: a source unit may be a file or a bundle. Every
bundle member is hashed independently, and changing any member refreshes the
owning capsule without creating a second representation of every helper.

### PDE audit bundle

`research/pde_audit/` is executable verification infrastructure rather than a
paper family. Its simulation README defines a staged pipeline with structural
checks, target-blind frozen packets, post-hoc evaluation, obstruction
diagnostics, and an explicit physical-export guard.

Initial treatment: one audit-project capsule plus catalog entries for declared
top-level runners. Selected tracked reports are evidence dependencies. The
gitignored `simulation/output/` tree is never ingested.

## Software structure

### `software/em_charge_attribute/`

This is a collection of paired SymPy/Wolfram construction and audit pipelines,
comparators, run wrappers, directives, results, and curated reports. It does
not currently have a project README.

Initial treatment:

- Use `run_*.sh` files as the first entry-point inventory.
- Group paired engines, comparators, fixtures, schemas, and results behind the
  owning run entry point.
- Summarize top-level `*_result.md` and curated `reports/*.md` as result
  sources where they represent a current or historical verdict.
- Do not catalog prompts, every helper, or ignored generated artifact trees.

This project will likely benefit from a small generated script-domain page
rather than dozens of independent script pages.

### `software/force_visualizer/`

This is a compact visualization application with three clear entry points:

- `app.py` for the interactive application.
- `render_all.py` for headless GIF generation.
- `report.py` for numerical checks.

Initial treatment: one software-project capsule and one script catalog section.
Tracked GIFs are outputs, not semantic inputs. The README's research-paper
source map is retained, while its references to excluded root notes are not
promoted into the memory corpus.

### `software/stage1_solver/`

This is a full numerical solver and validation program with source modules,
tests, decisions, directives, reports, tools, frozen packets, and many ignored
run/scratch artifacts. `STAGE1_VERDICT.md` is a high-value current result and
the README documents both scope and executable harnesses.

Initial treatment:

- One software-project capsule from `README.md`.
- One result capsule for `STAGE1_VERDICT.md` and its cited frozen evidence.
- Script catalog entries for public harnesses and for coherent `pathA_*`
  tool/report families.
- Selected result reports as sources; directives and tests as dependencies or
  provenance, not standalone memory pages.
- Ignored `runs/`, `data/`, `figures/`, `output/`, and `_scratch/` never enter
  discovery.

## Initial corpus size and shape

The executable configuration now contains 32 initial source units:

- 16 monolithic paper capsules.
- 1 v3 governance capsule.
- 1 v3 paper capsule.
- Explicit S9, S10, S11, unified S11b, and in-progress S11c step-family
  capsules, with large outputs and export modules hashed identity-only.
- 1 PDE-audit project capsule.
- 3 software-project capsules.
- Selected result capsules for Stage 1 and four EM/charge pipelines.
- Script catalogs grouped by paper or software project.

This stays within the planned 20–50 source-unit pilot while covering papers,
componentized research, paired engines, explicit measurements, a supersession
lineage, and executable software.

## Execution order from here

1. Finish Phase 0 status/evidence rules and the Atlas migration inventory.
2. Define the exact source-unit schema, including bundle membership and
   dependency roles.
3. Build the deterministic scanner against this config and the committed Git
   tree/blob snapshot.
4. Run `memory status` before any semantic extraction; verify that it proposes
   only the source units above.
5. Pilot extraction on one monolithic paper, one v3 step family, and one
   software result family.
6. Review the pilot for citation accuracy, status restraint, and usefulness.
7. Populate the remaining initial corpus.
8. Migrate Atlas firewalls/open gates/source crosswalks into the resulting
   topics.
9. Pass the retrieval benchmark and legacy-deletion gate before removing
   `atlas/` and `graph/`.

No additional user decision is required to start these steps. Any ambiguity
about scientific precedence will be recorded as a candidate conflict rather
than guessed.
