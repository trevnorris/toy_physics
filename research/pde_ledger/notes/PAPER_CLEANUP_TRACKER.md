# PDE Ledger Paper Cleanup Tracker

## Purpose

This document is the working tracker for turning `research/pde_ledger/paper/`
from a generated monolithic derivation dump into a usable internal ledger and a
clean source base for future publishable papers.

The archive is not disposable scaffolding. It is the master reference for the
derivations accumulated during the moving-throat PDE exploration and is meant
to let future papers cite the completed derivation chain without re-deriving
the same material from scratch.

It is not a proof-status ledger for the mathematics itself.  It is a
document/product tracker for:

- readability,
- structure,
- provenance,
- build layout,
- and editorial cleanup.

Use this file as the master checklist for paper cleanup decisions and status.

## Authorial Charter

Clarified on `2026-04-20`:

1. `paper/pde_ledger.tex` is the authoritative master document for the PDE
   derivation program.
2. Its primary role is to preserve the derivation chain and make it citable by
   future papers.
3. Readability work should improve navigation without damaging archive value or
   provenance.
4. The verification stack is part of the product, not an optional sidecar:
   notes, SymPy audits, Mathematica audits, and saved outputs all contribute to
   trust in the master document.
5. Cleanup should therefore optimize for:
   - archival stability,
   - human navigation,
   - provenance clarity,
   - and verification traceability.

## Repo Snapshot

Snapshot taken on `2026-04-20`.

Observed facts:

| Item | Path | Observed state |
|---|---|---|
| Main TeX entry | `research/pde_ledger/paper/pde_ledger.tex` | `book` document with `frontmatter`, `mainmatter`, and appendices |
| Archive PDF size | `research/pde_ledger/paper/pde_ledger.pdf` | 682 pages |
| Reader PDF size | `research/pde_ledger/paper/pde_ledger_reader.pdf` | 187 pages |
| Frontmatter source size | `paper/frontmatter/*.tex` | 4,427 lines |
| Main parts source size | `paper/parts/*.tex` | 11,638 lines |
| Appendix source size | `paper/appendices/*.tex` | 19,203 lines |
| Stage file source size | `paper/stages/*.tex` | 6,372 lines across 253 files |
| Notes corpus | `research/pde_ledger/notes/` | master notes plus stage notes plus review docs |
| Verification scripts | `research/pde_ledger/scripts/` | stage SymPy audits plus numerical stress scripts |
| Mathematica mirror | `research/pde_ledger/mathematica/` | existing second-CAS audit layer plus saved outputs |

## Immediate Interpretations

1. This is not currently a conventional paper.
2. It is a mixed artifact: authoritative derivation ledger, provenance archive,
   verification index, and publication wrapper all at once.
3. The `.md` files in `notes/` are source artifacts, not reader-facing prose.
4. The `scripts/` tree is part of the verification archive, not the narrative.
5. The `paper/stages/` directory appears to be generated template inventory,
   not the canonical stage narrative currently used by the build.
6. The canonical built stage narrative currently lives in
   `paper/appendices/stage_appendix_part*.tex`.

## Why Chapter 1 Looks Like It Starts "Late"

The current document uses:

- `\frontmatter`
- `\mainmatter`
- `\appendix`

So Chapter 1 starts on printed page `1`, but only after a very large frontmatter
block with Roman numeral pages.  In a PDF viewer this makes Chapter 1 appear
far into the file even though the chapter numbering resets correctly.

This is a product-design problem, not a LaTeX bug.

## Working Diagnosis

The main cleanup problem splits into three tracks:

| Track | Question |
|---|---|
| Product definition | What document are we actually trying to ship? |
| Editorial structure | What should a human read in the PDF, and what should be archived elsewhere? |
| Provenance and audit | How do notes, scripts, outputs, and stage claims map cleanly to the ledger? |

We should not continue making local prose edits until the product definition is
stable enough to guide them.

## Product Hierarchy

The master archive remains the top-level product. Other outputs should be
derived from it rather than replacing it.

| Product | Role |
|---|---|
| Master archive ledger | authoritative derivation reference for future papers |
| Reader ledger | human-usable navigation layer over the same derivation program |
| Future papers | narrow downstream papers that cite the master archive |

The cleanup program should preserve this hierarchy.

## Verification Baseline

The master document is only as trustworthy as the verification apparatus
attached to it, so cleanup planning has to account for the actual audit state.

Current observed verification layers:

| Layer | Path | Observed state |
|---|---|---|
| Derivation notes | `research/pde_ledger/notes/` | canonical source/archive inputs |
| SymPy audits | `research/pde_ledger/scripts/` | broad stagewise audit corpus |
| Mathematica audits | `research/pde_ledger/mathematica/` | large existing mirror, not a greenfield task |
| Numerical stress | `research/pde_ledger/scripts/numerical/` and `research/pde_ledger/mathematica/numerical/` | existing stress infrastructure |

Observed Mathematica state on `2026-04-20`:

- `research/pde_ledger/mathematica/` already exists,
- it contains a large set of stagewise `.wl` audits,
- it has runner scripts and saved outputs,
- and `research/pde_ledger/mathematica/output/_summary.txt` currently reports
  `TOTAL: 167  PASS: 167  FAIL: 0  SKIPPED: 0`.

This changes the practical roadmap:

- we do not need to invent the Mathematica layer from scratch,
- we do need to audit its scope, trust model, and coverage relative to the
  notes and SymPy layers,
- and paper-side provenance should reflect that dual-CAS state cleanly.
- the detailed stage-coverage inventory is intentionally deferred until the
  archive/reader paper cleanup baseline is in place.

## Target Outputs

We likely need three distinct outputs, even if they share sources.

| Output | Audience | Purpose | Status |
|---|---|---|---|
| Archive ledger | internal and future-author use | preserve the authoritative derivation and audit trail | canonical |
| Reader ledger | internal collaborators | readable theorem/provenance companion | needed |
| Future papers | external | narrow, publication-scale claims | needed |

The current `pde_ledger.tex` is the archive ledger, but it is also carrying too
much reader-onboarding and publication-shaping burden.

## Recommended Product Decision

Default recommendation:

1. Keep `pde_ledger.tex` as the full authoritative archive build.
2. Create one slimmer reader build.
3. Treat future papers as separate downstream products, not edited versions of
   the 700+ page archive.

This avoids destroying provenance while still giving us something humans can
actually work with.

## Workstreams

### 1. Define the Canonical Role of Each Tree

Goal:
decide which directories are source-of-truth, generated artifacts, or archival.

Required decisions:

- `notes/`: canonical derivation source, or raw AI output archive?
- `scripts/`: canonical verification layer, or one audit implementation among several?
- `paper/parts/`: canonical human-readable theorem narrative?
- `paper/appendices/stage_appendix_part*.tex`: canonical per-stage ledger?
- `paper/stages/`: keep, regenerate, repurpose, or retire?

Acceptance criteria:

- every top-level subtree has one clear role,
- no file class is maintained in two places without explicit generation rules,
- redundant generated stubs are either retired or explicitly marked as such.

### 2. Split Archive Build from Reader Build

Goal:
stop forcing one PDF to satisfy incompatible use cases.

Recommended direction:

- keep `pde_ledger.tex` as the archive build,
- create a second entry point for a reader-oriented build.

Reader build should likely include:

- title page,
- short purpose/scope note,
- concise reading guide,
- parts 1--8,
- a short provenance appendix,
- a short reproducibility appendix.

Reader build should likely exclude:

- full stage ledger tables,
- full stage appendices,
- fill workflow,
- source file index,
- raw inventory-style audit material that is only useful for archival tracing.

Acceptance criteria:

- one PDF is readable end-to-end by a human collaborator,
- one PDF preserves the full archival record,
- both compile from the same repo without ad hoc editing.

### 3. Compress the Frontmatter Aggressively

Goal:
reduce frontmatter from a mini-book to an onboarding layer.

Current problem:

- the frontmatter is long enough to dominate the start of the PDF,
- much of it repeats policy that could live in an editorial guide instead.

Recommended actions:

- collapse duplicated reading-policy material,
- move detailed citation/status doctrine to an appendix or editorial guide,
- keep only the minimum required to orient a reader.

Acceptance criteria:

- reader build reaches Chapter 1 quickly,
- frontmatter introduces the ledger without exhausting the reader,
- rules remain available somewhere, but not all in the first block of pages.

### 4. Make the Main Parts Read Like Human Prose

Goal:
turn the main body into a navigable theorem ledger rather than a raw model dump.

Required improvements:

- add one short abstract at the start of each part,
- add a "what this part proves / does not prove" block,
- add a dependency summary at the top of each part,
- reduce repeated status language inside theorem exposition,
- move repetitive stage history back to appendices.

Acceptance criteria:

- a reader can skim the part openings and understand the global flow,
- theorem blocks are not buried under repeated policy prose,
- the main parts can be read without treating the document as a database dump.

### 5. Reframe Notes and Script References

Goal:
stop exposing raw AI-oriented filenames as if they were reader-facing sources.

Current problem:

- the paper references `.md` and script artifacts whose names are only useful to
  the machine-generation pipeline,
- this makes the document feel like an internal export log rather than a curated ledger.

Recommended actions:

- keep raw filenames in provenance tables only,
- in the main text, refer to curated source classes and stage IDs,
- create one human-readable source/provenance map that translates:
  stage ID -> note file -> script file -> output file -> paper section.

Acceptance criteria:

- the main text references stable human concepts,
- raw artifact names remain available for audit,
- provenance can be followed without cluttering narrative sections.

### 6. Decide What to Do with `paper/stages/`

Goal:
remove ambiguity around the 253 generated stage `.tex` files.

Current observed state:

- sampled stage files are template stubs,
- the main build does not include them,
- the active stage narrative already exists in `appendices/stage_appendix_part*.tex`.

Likely options:

- archive them as generated templates and stop editing them,
- use them as a future machine-readable source layer,
- or regenerate stage appendices from them and make that pipeline explicit.

Default recommendation:

- do not manually maintain both `paper/stages/` and stage appendices,
- treat `paper/stages/` as archival/generated until a real generation pipeline exists.

Acceptance criteria:

- there is one clear canonical stage narrative layer,
- dead or redundant stage files are not mistaken for active source.

### 7. Clean Up Provenance and Verification Metadata

Goal:
make the audit apparatus useful rather than overwhelming.

Required improvements:

- unify source-file references,
- normalize stage status language,
- maintain one reproducibility map,
- maintain one dependency ledger,
- maintain one assumption ledger.

Related existing docs worth reusing:

- `notes/review/ASSUMPTION_LEDGER.md`
- `notes/review/DEPENDENCY_LEDGER.md`
- `notes/review/PROOF_HARDENING_PLAN.md`

Acceptance criteria:

- provenance is centralized,
- verification docs complement the paper rather than spilling into it,
- stage status is consistent across the repo.

### 8. Verification Governance

Goal:
tie the paper cleanup plan to the actual trust model of the derivation program.

Required improvements:

- define what counts as the canonical claim for a stage:
  note, TeX ledger, SymPy audit, or Mathematica audit,
- inventory mismatches between note coverage, SymPy coverage, and Mathematica coverage,
- identify which stages are load-bearing for future papers,
- prioritize secondary audit review for stages whose scripts may be correct
  syntactically but weak epistemically,
- make saved audit outputs citable from the archive without flooding the main narrative.

Acceptance criteria:

- every load-bearing stage has an explicit verification status,
- SymPy and Mathematica coverage can be compared from one place,
- weak or tautological audits are discoverable,
- the master ledger can cite verification evidence without embedding raw transcripts.

### 9. Build Hygiene and Navigation

Goal:
make the build predictable and the PDF easier to navigate.

Needed items:

- standardize on `latexmk` for normal builds,
- decide whether the archive PDF should include bookmarks only or also deep TOC,
- add a short build note describing archive vs reader builds,
- keep unresolved refs and duplicate labels at zero.

Acceptance criteria:

- normal build path is documented,
- PDF navigation is usable,
- log is clean enough that new real problems are visible.

## Priority Queue

### Phase 0: Decision Pass

Status: `not started`

Tasks:

- decide whether `pde_ledger.tex` remains the archive build,
- record explicitly that the archive build is the authoritative master reference,
- decide whether we want a second reader build now,
- decide whether `paper/stages/` is active or archival,
- decide whether raw note/script filenames belong only in provenance appendices.

### Phase 1: Structural Refactor

Status: `in progress`

Tasks:

- create reader build entry point,
- trim frontmatter for reader build,
- move heavy provenance material out of the reader main path,
- preserve archive build unchanged except for hygiene fixes.

### Phase 2: Narrative Cleanup

Status: `in progress`

Tasks:

- add per-part abstracts,
- add dependency and claim-scope summaries,
- remove repetitive policy text from theorem flow,
- improve chapter openers and transitions.

Current note:

- Part I now uses the pilot opener format: short archive-facing introduction,
  explicit reader map, and grouped stage packets.
- Part II now uses the same opener pattern, replacing the row-by-row stage dump
  with grouped stage packets and a short reader map.
- Parts III--IV now use the same opener pattern, replacing large stage tables
  with grouped packets and short reader maps.
- Parts V--VIII now use the same opener pattern, replacing the generated
  source/status longtables with grouped packets and short reader maps.
- Next readability passes should remove repeated policy text from theorem flow
  and tighten transitions without compressing mathematical content.

### Phase 3: Provenance Cleanup

Status: `in progress`

Tasks:

- create a human-readable source map,
- standardize stage/source/script naming references,
- centralize audit references.

Current note:

- Part I is now the pilot for canonical stage sourcing: the appendix assembles
  `paper/stages/stage_001.tex` through `stage_023.tex` instead of carrying the
  derivations inline.
- Part II now uses the same pattern: the appendix assembles
  `paper/stages/stage_024.tex` through `stage_036.tex` instead of carrying the
  derivations inline.
- Part III now uses the same pattern: the appendix assembles
  `paper/stages/stage_037.tex` through `stage_090.tex` instead of carrying the
  derivations inline.
- Part IV now uses the same pattern for its stage cards: the appendix assembles
  `paper/stages/stage_091.tex` through `stage_163.tex` while keeping the
  audit-path derivation sections inline above the canonical stage cards.
- Part V now uses the same pattern for its stage cards: the appendix assembles
  `paper/stages/stage_164.tex` through `stage_200.tex` while keeping the
  audit-path derivation sections inline above the canonical stage cards.
- Part VI now uses the same pattern for its stage cards: the appendix assembles
  `paper/stages/stage_201.tex` through `stage_218.tex` while keeping the
  theorem-path derivation sections inline above the canonical stage cards.
- Part VII now uses the same pattern for its stage cards: the appendix assembles
  `paper/stages/stage_219.tex` through `stage_242.tex` while keeping the
  theorem-path derivation sections inline above the canonical stage cards.
- Part VIII now sources `paper/stages/stage_243.tex` through
  `stage_253.tex`; unlike Parts I--VII, those stage files carry the full
  derivation sections rather than compact stage cards.
- The raw note provenance for migrated stages now lives in
  `STAGE_PROVENANCE_INDEX.md`, while the stage files themselves carry only the
  executable verification references needed for human review.
- Format normalization is now decided: keep the mixed model. Parts I--VII use
  compact stage cards, while Part VIII uses full derivation-stage files.
- Verification coverage is now tracked in `STAGE_VERIFICATION_COVERAGE.md`,
  which gives the audit baseline without pushing raw artifact names back into
  the PDF.

### Phase 4: Verification Reconciliation

Status: `in progress`

Tasks:

- inventory stage coverage across notes, SymPy, and Mathematica,
- classify which existing audits are load-bearing vs routine,
- audit constant provenance and prohibit unexplained literals in verification scripts,
- identify weak, tautological, or stale audit implementations,
- decide how the archive ledger should expose verification status.

Current note:

- The mixed stage-source model is now intentional rather than provisional.
- The stage coverage baseline is now recorded in
  `STAGE_VERIFICATION_COVERAGE.md`.
- The future-paper citation-support subset is now defined in
  `CITATION_SUPPORT_SET.md`.
- Stages `001--002` now have dedicated SymPy audits and an initial
  constant-provenance log in `CHECKPOINT_CONSTANT_PROVENANCE.md`.
- Stages `090` and `096` now also have dedicated SymPy checkpoint audits, so
  the support-set Mathematica-only outliers are gone.
- Stage `200` now has a Mathematica mirror and a first review note, so the
  late-stage support frontier has narrowed again.
- Stage `203` now has a Mathematica mirror and a first review note too.
- Stage `218` now also has a Mathematica mirror and a first review note.
- Stage `221` now also has a Mathematica mirror and a first review note.
- Stage `239` now also has a Mathematica mirror and a first review note.
- Stage `242` now also has a Mathematica mirror and a first review note.
- Stage `243` now also has a Mathematica mirror and a first review note.
- Stage `248` now also has a Mathematica mirror and a first review note.
- Stage `253` now also has a Mathematica mirror and a first review note, so the
  late-stage support frontier is now closed.
- Stages `001--002` now also have Mathematica parity, so the checkpoint
  support set has full symbolic parity and review coverage.
- The checkpoint trust baseline now lives in `CHECKPOINT_TRUST_AUDIT.md`:
  `25` `strong`, `0` `moderate`, `0` `weak`.
- The next verification work is no longer parity; the checkpoint support set is
  symbolically closed and `248` / `253` now also have dedicated numerical
  stress.
- The paper-side verification surface is now exposed compactly through the
  frontmatter and reader verification summary, while raw per-stage script
  manifests are intentionally omitted from printed stage cards.

### Phase 5: Visual and Typographic Cleanup

Status: `not started`

Tasks:

- reduce overfull boxes where materially harmful,
- clean hyperref title warnings,
- review heading density,
- review table density and page-breaking strategy.

## Task Board

| ID | Task | Category | Status | Notes |
|---|---|---|---|---|
| P0-00 | Record archive charter in paper-side docs | product | open | authoritative master reference |
| P0-01 | Define archive vs reader build split | product | done | archive is authoritative; reader build added |
| P0-02 | Decide canonical status of `paper/stages/` | structure | done | all Parts I--VIII now source canonical stage files |
| P0-03 | Decide whether raw artifact names stay out of main narrative | editorial | done | migrated Stages 001--253 keep raw note names out of the PDF and in the repo-local provenance index |
| P1-01 | Add second TeX entry point for reader build | build | done | `pde_ledger_reader.tex` added |
| P1-02 | Compress frontmatter for reader build | editorial | done | reader build uses short orientation frontmatter |
| P1-03 | Keep archival appendices only in archive build | build | done | reader build keeps compact summary appendices only |
| P2-00 | Add per-part claim-scope summaries | readability | done | all Parts I--VIII now open with proves / does-not-prove / later-parts-need-it summaries |
| P2-01 | Tighten per-part opener prose into explicit abstract blocks | readability | done | all Parts I--VIII now use the short opener/readermap format |
| P2-02 | Add per-part dependency blocks | readability | done | reader-map format now covers Parts I--VIII |
| P2-03 | Remove repetitive policy text from theorem flow | readability | open | move repeated status doctrine out of chapter interiors |
| P2-04 | Remove per-stage forced page breaks from archive appendices | layout | done | Parts I--VII now neutralize stage-card `\clearpage` calls locally so derivations flow continuously without rewriting every canonical stage file |
| P3-01 | Create source/provenance map | provenance | done | `STAGE_PROVENANCE_INDEX.md` now covers Stages 001--253 |
| P3-02 | Normalize stage/source references | provenance | open | avoid raw filenames in prose |
| P4-01 | Inventory note/SymPy/Mathematica coverage | verification | done | baseline inventory now lives in `STAGE_VERIFICATION_COVERAGE.md` |
| P4-02 | Define verification status schema for load-bearing stages | verification | done | checkpoint trust baseline now lives in `CHECKPOINT_TRUST_AUDIT.md` |
| P4-03 | Audit constant provenance for load-bearing stages | verification | open | every constant must be derived, source-anchored, or probe-only |
| P4-04 | Define future-paper citation support subset | verification | done | control subset now lives in `CITATION_SUPPORT_SET.md` |
| P4-05 | Promote foundational checkpoints to dedicated audits | verification | done | Stages `001--002` now have dedicated SymPy scripts and provenance entries |
| P4-06 | Reconcile checkpoint Mathematica-only outliers | verification | done | Stages `090` and `096` now have dedicated SymPy mirrors and checkpoint provenance entries |
| P4-07 | Harden first late-stage support checkpoint | verification | done | Stage `200` now has a Mathematica mirror, provenance entry, and review note |
| P4-08 | Harden second late-stage support checkpoint | verification | done | Stage `203` now has a Mathematica mirror, provenance entry, and review note |
| P4-09 | Harden third late-stage support checkpoint | verification | done | Stage `218` now has a Mathematica mirror, provenance entry, and review note |
| P4-10 | Harden fourth late-stage support checkpoint | verification | done | Stage `221` now has a Mathematica mirror, provenance entry, and review note |
| P4-11 | Harden fifth late-stage support checkpoint | verification | done | Stage `239` now has a Mathematica mirror, provenance entry, and review note |
| P4-12 | Harden sixth late-stage support checkpoint | verification | done | Stage `242` now has a Mathematica mirror, provenance entry, and review note |
| P4-13 | Harden seventh late-stage support checkpoint | verification | done | Stage `243` now has a Mathematica mirror, provenance entry, and review note |
| P4-14 | Harden eighth late-stage support checkpoint | verification | done | Stage `248` now has a Mathematica mirror, provenance entry, and review note |
| P4-15 | Harden ninth late-stage support checkpoint | verification | done | Stage `253` now has a Mathematica mirror, provenance entry, and review note |
| P4-16 | Add Mathematica parity to foundational checkpoints | verification | done | Stages `001--002` now have Mathematica mirrors and refreshed review records |
| P4-17 | Baseline checkpoint trust audit | verification | done | stage-level checkpoint tiers now live in `CHECKPOINT_TRUST_AUDIT.md` |
| P4-18 | Repair Stage 185 trust defect | verification | done | Stage `185` now rebuilds the monomial laws from primitive-ratio compilers in both CAS layers and has a current PASS review |
| P4-19 | Resolve Stage 001/002 convention caveat | verification | done | canonical Stage `001--002` cards now state the densitized convention explicitly and the current review status is PASS |
| P4-20 | Repair Stage 003 checkpoint drift | verification | done | shared numerical-stress harness now resolves the repo-local sample path, the note/review bookkeeping is current, and Stage `003` is promoted to `strong` |
| P4-21 | Reconcile Stage 022 round-trip drift | verification | done | existing `N0/N2/N4` round-trip checks were confirmed in both CAS layers, the stale review caveat was removed, and Stage `022` is promoted to `strong` |
| P4-22 | Reconcile Stage 023 assembly drift | verification | done | existing one-port `Z_n/N_n` reconstruction checks were confirmed in both CAS layers, the stale assembled-input caveat was removed, and Stage `023` is promoted to `strong` |
| P4-23 | Reconcile Stage 024 overlap drift | verification | done | the explicit `H_r` rename note and existing unequal-lane witness checks were confirmed in both CAS layers, the stale review caveats were removed, and Stage `024` is promoted to `strong` |
| P4-24 | Reconcile Stage 036 boundary drift | verification | done | the live `xi >= 0` / `0 <= xi < 1` boundary assumptions were confirmed in both CAS layers, the stale onset-boundary caveat was removed, and Stage `036` is promoted to `strong` |
| P4-25 | Reclassify Stage 089 explicit verdict | verification | done | the explicit Family-1 verdict was confirmed in both CAS layers as a closed arithmetic theorem conditional on the upstream minimal module, and its carried-threshold provenance is now logged |
| P4-26 | Reclassify Stage 090 status boundary | verification | done | the reduced theorem-status boundary is now treated as a narrow but citation-grade checkpoint because its carried inputs are explicit, source-anchored, and replayed in both CAS layers |
| P4-27 | Reclassify Stage 163 off-family transport packet | verification | done | the theorem path was confirmed symbolic in both CAS layers, the Family-1 readbacks were isolated as explanatory outputs only, and Stage `163` is promoted to `strong` |
| P4-28 | Add numerical stress to Stages 248 and 253 | verification | done | shared JSON-driven Python + Mathematica stress harnesses now cover both late-stage checkpoints without adding committed output artifacts |
| P4-29 | Add paper-side verification routing summary | verification | done | the archive frontmatter and reader verification appendix now summarize the naming convention/baseline compactly, and printed stage cards omit raw verification filename manifests |
| P4-30 | Red-team audit pass on batch I.1 (stages 001-012) | verification | done | 2026-05-21: 12 stages reached `verified` with both engines independently checking claims; 9 new `.wl` mirrors created under `mathematica/`; tautological assertions replaced with real ones; `material_change: false` on every stage. See `redteam/batches/batch_I1.md`. |
| P4-31 | Red-team audit pass on batch I.2 (stages 013-023) | verification | done | 2026-05-21: 11 stages reached `verified` with both engines independently checking claims; 8 new native `.wl` mirrors created under `mathematica/` for stages 013-020 (using `SphericalHarmonicY`, `ThreeJSymbol`, `LinearSolve`, `Series`+`Coefficient`, `SphericalHankelH1`, `EulerEquations`); 3 mirrors (stages 021-023) rewritten off the SymPy algebraic path; tautological assertions across 8 stages (variable-independence traps, `A+(-A)=0`, solve-roundtrip self-substitution, hardcoded `m=0` Gaunt short-circuits) replaced with substantive checks; `material_change: false` on every stage. See `redteam/batches/batch_I2.md`. |
| P4-32 | Sweep for multi-line `lRed = ...` continuation defect class | verification | open | recurring across I.1 stage 003 and I.2 stage 021 (both caught only kinetic terms because the `lRed` RHS spans multiple lines without parentheses); audit all stage `.wl` files for silent continuation truncation before later batches surface the same defect again |
| P4-33 | Red-team audit pass on batch II.1 (stages 024-036) | verification | done | 2026-05-22: 13 stages reached `verified` with both engines independently checking load-bearing claims; 43 findings closed (`tautological_check` ~22, `mathematica_transliteration` 13, `insufficient_verification` ~4, `hardcoded_result` 2); `mathematica_transliteration` dominated the batch -- every stage's Mathematica mirror was a line-by-line port of the SymPy companion and had to be rewritten in native CAS idioms (`Integrate`, `Eigenvalues`, `LinearSolve`, `Solve`, `Factor`/`Apart`, `Limit`, `Reduce`); Stage 025 hardcoded-target `54/5` reanchored to Stage 023, Stage 031 9-term `radcrit` polynomial derived from `T0^2 R^2.subs(alpha, alpha_crit)`, Stage 031 SymPy abstract `Function` replaced with physical `s_-/lam_-`; zero codex iter-2 fixes needed (vs 1 in each of I.1 and I.2), partly because the mid-batch `fix_loop.sh` sanity-exec fix (commit `3534b80`) now actually runs refreshed transcripts and greps for `FAIL`/`Traceback`/`AssertionError`/`$Failed`, and `codex.md` now documents the `lRed` multi-line continuation defect and the `D[expr, f[t]] = 0` Mathematica quirk. `material_change: false` on every stage. See `redteam/batches/batch_II1.md`. |
| P3-03 | Propagate batch I.1 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-21" sentence to `paper/stages/stage_001.tex` through `stage_012.tex` |
| P3-04 | Propagate batch I.2 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-21" sentence to `paper/stages/stage_013.tex` through `stage_023.tex` |
| P3-05 | Propagate batch II.1 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-22" sentence to `paper/stages/stage_024.tex` through `stage_036.tex` |
| P4-34 | Red-team audit pass on batch III.1 (stages 037-048) | verification | done | 2026-05-22: 12 stages reached `verified` with both engines independently checking load-bearing claims; 27 findings closed across 10 dirty stages (`mathematica_transliteration` 10, `tautological_check` ~11, `insufficient_verification` 4, `hardcoded_result` 1, `symbol_assumption_error` 1); stages 042 and 048 verified clean on first read (042's closed-form-identity structure and 048's independent `Solve`/`Series` route already broke the line-by-line correspondence). Native CAS rewrites used: `Together`/`Apart` numerator-denominator extraction, `Solve[Det[...]]`/`NullSpace` eigenvector routes, `Series`+`Coefficient`, `Reduce[ForAll[...]]` sign claims, `PolynomialQuotientRemainder` factor checks, `Limit`-based endpoint asymptotics; stage 037 the closed-form Schur reconstruction was switched to "recover `xi`/`alpha` from two `sigmaWall` entries and cross-check the third"; stage 040 the eigenvector self-substitution was replaced with a perturbed-matrix residual against `M - alpha z z^T`. Zero codex iter-2 fixes needed (matches II.1). `material_change: false` on every stage. See `redteam/batches/batch_III1.md`. |
| P3-06 | Propagate batch III.1 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-22" sentence to `paper/stages/stage_037.tex` through `stage_048.tex` |
| P5-01 | Typography and overflow cleanup pass | formatting | open | after structure stabilizes |

## Questions That Must Be Answered Early

1. Is the primary goal an internal archive, a collaborator-readable companion,
   or something close to a preprint?
2. Do we want to preserve the full stage-by-stage material in the main PDF, or
   is it enough to preserve it in a separate archive build?
3. Are the `.md` notes canonical source material, or should the canonical human
   source move into the TeX files?
4. Is there any real consumer for the `paper/stages/` template files today?
5. Do we want one unified provenance map in the paper, or do we want the paper
   to link out to a repo-local audit document instead?
6. How should the archive expose verification status for a stage when the note,
   SymPy audit, and Mathematica audit do not all have the same strength?

## Recommended First Concrete Editing Pass

Recommended next implementation sequence:

1. keep the chosen mixed model: compact stage cards in Parts I--VII and full
   derivation-stage files in Part VIII,
2. continue provenance cleanup so stage files carry only the executable
   verification attachments intended for human review,
3. use `CITATION_SUPPORT_SET.md` as the hardening target for future-paper
   support,
4. audit constant provenance and no-magic-numbers compliance in that subset,
5. then audit SymPy and Mathematica trust strength against that subset.

This keeps the archive usable while moving it toward a single-source,
referee-checkable stage ledger.

## Change Log

| Date | Change |
|---|---|
| 2026-04-20 | Initial tracker created from repo survey and current build audit. |
| 2026-04-20 | Added separate reader build with compact frontmatter and summary appendices; archive build remains authoritative. |
| 2026-04-20 | Added reader-facing claim-scope summaries to all eight main parts and recompiled both archive and reader PDFs successfully. |
| 2026-04-20 | Reworked the Part I opener into a reusable reader-map and grouped stage-packet format, then recompiled both archive and reader PDFs successfully. |
| 2026-04-20 | Propagated the same opener pattern to Part II and recompiled both archive and reader PDFs successfully. |
| 2026-04-20 | Migrated Part V stage cards into canonical `paper/stages/stage_164.tex` through `stage_200.tex`, moved raw note provenance into `STAGE_PROVENANCE_INDEX.md`, and kept the audit-path derivation sections inline above the canonical cards. |
| 2026-04-20 | Migrated Part VI stage cards into canonical `paper/stages/stage_201.tex` through `stage_218.tex`, moved raw note provenance into `STAGE_PROVENANCE_INDEX.md`, and kept the theorem-path derivation sections inline above the canonical cards. |
| 2026-04-20 | Migrated Part VII stage cards into canonical `paper/stages/stage_219.tex` through `stage_242.tex`, moved raw note provenance into `STAGE_PROVENANCE_INDEX.md`, and kept the theorem-path derivation sections inline above the canonical cards. |
| 2026-04-20 | Migrated Part VIII into canonical `paper/stages/stage_243.tex` through `stage_253.tex`, moved raw note provenance into `STAGE_PROVENANCE_INDEX.md`, and kept the appendix as an assembler over full derivation-stage files plus downstream guardrails. |
| 2026-04-20 | Propagated the same opener pattern to Parts III--IV; archive and reader builds now stand at 707 and 193 pages respectively. |
| 2026-04-20 | Promoted Part I Stages 001--023 into canonical `paper/stages/` files, converted the Part I appendix into an aggregator, attached stage-level audit references, and recompiled both archive and reader PDFs successfully. |
| 2026-04-20 | Propagated the same opener pattern to Parts V--VIII, removed the generated source/status longtables from those chapter heads, and kept the mathematical derivations intact. |
| 2026-04-20 | Promoted Part II Stages 024--036 into canonical `paper/stages/` files, converted the Part II appendix into an aggregator, removed raw note filenames from migrated stage pages, and added the repo-local `STAGE_PROVENANCE_INDEX.md` for Stages 001--036. |
| 2026-04-20 | Promoted Part III Stages 037--090 into canonical `paper/stages/` files, converted the Part III appendix into an aggregator, removed raw note filenames from migrated stage pages, and extended `STAGE_PROVENANCE_INDEX.md` through Stage 090. |
| 2026-04-20 | Promoted Part IV Stages 091--163 into canonical `paper/stages/` files, converted the Part IV stage-card block into an aggregator, removed raw note filenames from migrated stage pages, and extended `STAGE_PROVENANCE_INDEX.md` through Stage 163. |
| 2026-04-20 | Chose the mixed canonical stage model intentionally: Parts I--VII remain compact stage cards, Part VIII remains full derivation-stage files, and the audit baseline now lives in `STAGE_VERIFICATION_COVERAGE.md`. |
| 2026-04-20 | Added the constant-provenance / no-magic-numbers rule to the paper-side verification queue so future audit work must justify every carry-forward constant. |
| 2026-04-20 | Added `CITATION_SUPPORT_SET.md` to define the future-paper support subset and focus the next audit wave on its checkpoint stages. |
| 2026-04-20 | Added dedicated SymPy audits for Stages `001--002`, updated their provenance/verification fields, and logged the first checkpoint constant-provenance findings. |
| 2026-04-20 | Added dedicated SymPy audits for Stages `090` and `096`, removed the last support-set Mathematica-only checkpoint gaps, and extended the checkpoint constant-provenance log accordingly. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `200`, updated the coverage/provenance baselines, and narrowed the late-stage support-set frontier to `203`, `218`, `221`, `239`, `242`, `243`, `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `203`, updated the coverage/provenance baselines again, and narrowed the late-stage support-set frontier to `218`, `221`, `239`, `242`, `243`, `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `218`, updated the coverage/provenance baselines again, and narrowed the late-stage support-set frontier to `221`, `239`, `242`, `243`, `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `221`, updated the coverage/provenance baselines again, and narrowed the late-stage support-set frontier to `239`, `242`, `243`, `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `239`, updated the coverage/provenance baselines again, and narrowed the late-stage support-set frontier to `242`, `243`, `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `242`, updated the coverage/provenance baselines again, and narrowed the late-stage support-set frontier to `243`, `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `243`, updated the coverage/provenance baselines again, and narrowed the late-stage support-set frontier to `248`, `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `248`, kept the Session-II numerics isolated as benchmark-only specialization inputs, and narrowed the late-stage support-set frontier to `253`. |
| 2026-04-20 | Added a Mathematica mirror and first review note for Stage `253`, kept the physical-calibration numerics isolated as benchmark-only readbacks, and closed the late-stage support-set frontier. |
| 2026-04-20 | Added shared JSON-driven Python + Mathematica numerical-stress harnesses for Stages `248` and `253`, updated the checkpoint review/provenance baselines, and moved the late-stage support set past its last explicit numerical-hardening gap. |
| 2026-04-20 | Added Mathematica mirrors for Stages `001--002`, refreshed the foundational review records, and closed symbolic parity across the checkpoint support set. |
| 2026-04-20 | Reconciled Stage `023` review drift by confirming the existing one-port `Z_n/N_n` reconstruction checks in both CAS layers, logging its grouped-bundle constants, and promoting Stage `023` from `moderate` to `strong`. |
| 2026-04-20 | Reconciled Stage `024` review drift by confirming the explicit `H_r` rename note and the existing unequal-lane witness checks in both CAS layers, logging its exact overlap constants, and promoting Stage `024` from `moderate` to `strong`. |
| 2026-04-20 | Reconciled Stage `036` review drift by confirming the live `xi >= 0` / `0 <= xi < 1` boundary assumptions in both CAS layers, logging its support-frontier constants, and promoting Stage `036` from `moderate` to `strong`. |
| 2026-04-20 | Reclassified Stage `089` as `strong` after confirming in both CAS layers that the explicit Family-1 verdict is a closed arithmetic theorem once the minimal-isotropic module is accepted, and logged its carried-threshold provenance. |
| 2026-04-20 | Reclassified Stage `090` as `strong` after making the narrow status-boundary rule explicit, confirming the dual-CAS replay of the carried threshold verdict, and extending its provenance record to the Mathematica mirror. |
| 2026-04-20 | Reclassified Stage `163` as `strong` after confirming in both CAS layers that the off-family normal-coordinate theorem packet is symbolic end-to-end and that the Family-1 coefficient readbacks are explanatory outputs rather than proof-critical inputs. |
| 2026-04-20 | Added the compact paper-side verification routing summary to the archive frontmatter and reader appendix, suppressed raw verification filename manifests in printed stage cards, and rebuilt both PDFs successfully. |
| 2026-04-20 | Removed the archive's per-stage forced page-break policy by locally neutralizing stage-card `\clearpage` calls in Parts I--VII, reducing the archive PDF from 736 pages to 682 pages while keeping the reader build at 187 pages. |
| 2026-05-21 | Red-team audit completed for batch I.1 (stages 001-012). 12 stages reached `verified` with both engines (SymPy + Mathematica) independently checking load-bearing claims. 9 new native `.wl` mirrors created under `mathematica/`; 3 existing mirrors upgraded from transliteration to native derivation (using `EulerEquations`, `VariationalD`, `SphericalHarmonicY`, `Coefficient`/`Series`, `ThreeJSymbol`). Auditor wrote 25 findings across the batch (`tautological_check`, `mathematica_transliteration`, `missing_verification_script`, `insufficient_verification`, `hardcoded_result`, `stale_output`). All findings closed substantively with `material_change: false` on every stage. Codex caught 3 directive bugs via the Block-and-delta loop. See `redteam/batches/batch_I1.md`. |
| 2026-05-21 | Red-team audit completed for batch I.2 (stages 013-023, "Part I.2 -- Maxwell bridge, parent throat action, reduced one-port"). 11 stages reached `verified` with both engines independently checking load-bearing claims. 8 new native `.wl` mirrors created under `mathematica/` for stages 013-020 using Mathematica-native idioms (`SphericalHarmonicY`, `ThreeJSymbol`, `LinearSolve`, `Series`+`Coefficient`, `SphericalHankelH1`, `EulerEquations` from `VariationalMethods`); 3 mirrors (stages 021-023) rewritten off the SymPy algebraic path. Auditor wrote 29 findings (`tautological_check`, `mathematica_transliteration`, `missing_verification_script`, `hardcoded_result`, `insufficient_verification`); tautologies across 8 stages (variable-independence traps, `A+(-A)=0`, solve-roundtrip self-substitution, hardcoded `m=0` Gaunt ratio short-circuits) replaced with substantive checks. Stage 021 iter 1 regressed because Mathematica's `D[expr, qFun[t]]` zeros terms containing `qFun[t]`; iter 2 switched to `EulerEquations` and incidentally caught a recurrence of the I.1 stage 003 multi-line `lRed = ...` continuation defect. `material_change: false` on every stage. See `redteam/batches/batch_I2.md`. |
| 2026-05-22 | Red-team toolchain hardening: `fix_loop.sh` post-codex sanity exec now actually runs (previous awk pattern silently no-op'd against the multi-line YAML emitted by `$RT paths`), writes refreshed transcripts to canonical `scripts/output/` and `mathematica/output/` paths, and greps each for `FAIL`/`Traceback`/`AssertionError`/`$Failed` to catch the case where `math -script` exits 0 despite a printed FAIL (silent `expectZero` helpers). `codex.md` gained a "Common Mathematica pitfalls" section documenting the multi-line `lRed` continuation defect and the `D[expr, f[t]] = 0` quirk codex hit in I.2 stage 021. Commit `3534b80`. |
| 2026-05-22 | Red-team audit completed for batch II.1 (stages 024-036, "Part II.1 -- Overlap isotropy through continuum kernel"). 13 stages reached `verified` with both engines independently checking load-bearing claims. 43 findings closed (`tautological_check` ~22, `mathematica_transliteration` 13, `insufficient_verification` ~4, `hardcoded_result` 2). `mathematica_transliteration` was the dominant theme -- every stage's Mathematica mirror was a SymPy line-by-line port and had to be rewritten in native CAS idioms (`Integrate`, `Eigenvalues`, `LinearSolve`, `Solve`, `Factor`/`Apart`, `Limit`, `Reduce`, `Series`+`Coefficient`). Tautologies removed include definition-by-construction traps, `sp.solve` round-trip self-substitution, `FullSimplify[E]-E` zero checks, and hardcoded-on-both-sides assertions; the two `hardcoded_result` findings reanchored Stage 025's `54/5` overlap target to Stage 023 and Stage 031's 9-term `radcrit` polynomial to its upstream source; Stage 031's SymPy abstract `Function` was replaced with the physical `s_-`/`lam_-` symbols. Zero codex iter-2 fixes were needed (vs 1 in each of I.1 and I.2). MATHEMATICA_MIRROR_POLICY policy prose updated to make transliteration screening a named first-pass step. `material_change: false` on every stage. See `redteam/batches/batch_II1.md`. |
| 2026-05-22 | Red-team audit completed for batch III.1 (stages 037-048, "Part III.1 -- Continuum kernel, generalized branch, rank-2"). 12 stages reached `verified`; 10 needed codex edits (27 findings: 10 `mathematica_transliteration`, ~11 `tautological_check`, 4 `insufficient_verification`, 1 `hardcoded_result`, 1 `symbol_assumption_error`) and 2 verified clean on first read (042 and 048). Native CAS rewrites used across the batch: `Together`/`Apart` numerator-denominator extraction, `Solve[Det[...]]`/`NullSpace` eigenvector routes, `Series`+`Coefficient` derivations, `Reduce[ForAll[...]]` sign claims, `PolynomialQuotientRemainder` factor checks, `Limit`-based endpoint asymptotics; Stage 037's closed-form Schur reconstruction was switched to "recover `xi`/`alpha` from two `sigmaWall` entries and cross-check the third"; Stage 040's eigenvector self-substitution was replaced with a perturbed-matrix residual against `M - alpha z z^T`; Stage 047's `rho_0 - chi_0` / `sigma_0 - chi_0` tautologies (where `lamW/lamW` and `lamphi/lamphi` factors cancelled by construction) were closed and `mSupp`/`sEnhance` rewritten through independent routes. Zero codex iter-2 fixes (matches II.1, vs 1 each in I.1 and I.2). MATHEMATICA_MIRROR_POLICY policy prose updated to record that transliteration was the dominant pattern in III.1 too (10 of 12 stages). `material_change: false` on every stage. See `redteam/batches/batch_III1.md`. |
