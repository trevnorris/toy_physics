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
| P4-32 | Sweep for multi-line `lRed = ...` continuation defect class | verification | done | 2026-05-22: full sweep across 196 `.wl` files (182 top-level audit scripts + 14 numerical) flagged 4 multi-line-RHS blocks at low risk (all inside `Module[...]`/`Function[...]` brackets where newlines are whitespace, so the defect signature cannot fire there); no unprotected top-level instances of the defect survive. The two known prior cases — I.1 stage 003 (lines 54-67) and I.2 stage 021 (lines 69-74) — remain confirmed patched. Report at `redteam/lred_sweep_report.md`. The four flagged blocks (stage 001 `lapEig`, stage 185_187_orbit_stress `dmu`/`dmuExpected`) are defense-in-depth candidates only; not patched in this pass because they sit inside outer brackets that already prevent the defect. |
| P4-33 | Red-team audit pass on batch II.1 (stages 024-036) | verification | done | 2026-05-22: 13 stages reached `verified` with both engines independently checking load-bearing claims; 43 findings closed (`tautological_check` ~22, `mathematica_transliteration` 13, `insufficient_verification` ~4, `hardcoded_result` 2); `mathematica_transliteration` dominated the batch -- every stage's Mathematica mirror was a line-by-line port of the SymPy companion and had to be rewritten in native CAS idioms (`Integrate`, `Eigenvalues`, `LinearSolve`, `Solve`, `Factor`/`Apart`, `Limit`, `Reduce`); Stage 025 hardcoded-target `54/5` reanchored to Stage 023, Stage 031 9-term `radcrit` polynomial derived from `T0^2 R^2.subs(alpha, alpha_crit)`, Stage 031 SymPy abstract `Function` replaced with physical `s_-/lam_-`; zero codex iter-2 fixes needed (vs 1 in each of I.1 and I.2), partly because the mid-batch `fix_loop.sh` sanity-exec fix (commit `3534b80`) now actually runs refreshed transcripts and greps for `FAIL`/`Traceback`/`AssertionError`/`$Failed`, and `codex.md` now documents the `lRed` multi-line continuation defect and the `D[expr, f[t]] = 0` Mathematica quirk. `material_change: false` on every stage. See `redteam/batches/batch_II1.md`. |
| P3-03 | Propagate batch I.1 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-21" sentence to `paper/stages/stage_001.tex` through `stage_012.tex` |
| P3-04 | Propagate batch I.2 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-21" sentence to `paper/stages/stage_013.tex` through `stage_023.tex` |
| P3-05 | Propagate batch II.1 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-22" sentence to `paper/stages/stage_024.tex` through `stage_036.tex` |
| P4-34 | Red-team audit pass on batch III.1 (stages 037-048) | verification | done | 2026-05-22: 12 stages reached `verified` with both engines independently checking load-bearing claims; 27 findings closed across 10 dirty stages (`mathematica_transliteration` 10, `tautological_check` ~11, `insufficient_verification` 4, `hardcoded_result` 1, `symbol_assumption_error` 1); stages 042 and 048 verified clean on first read (042's closed-form-identity structure and 048's independent `Solve`/`Series` route already broke the line-by-line correspondence). Native CAS rewrites used: `Together`/`Apart` numerator-denominator extraction, `Solve[Det[...]]`/`NullSpace` eigenvector routes, `Series`+`Coefficient`, `Reduce[ForAll[...]]` sign claims, `PolynomialQuotientRemainder` factor checks, `Limit`-based endpoint asymptotics; stage 037 the closed-form Schur reconstruction was switched to "recover `xi`/`alpha` from two `sigmaWall` entries and cross-check the third"; stage 040 the eigenvector self-substitution was replaced with a perturbed-matrix residual against `M - alpha z z^T`. Zero codex iter-2 fixes needed (matches II.1). `material_change: false` on every stage. See `redteam/batches/batch_III1.md`. |
| P3-06 | Propagate batch III.1 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-22" sentence to `paper/stages/stage_037.tex` through `stage_048.tex` |
| P4-35 | Red-team audit pass on batch III.2 (stages 049-060) | verification | done | 2026-05-22: 12 stages reached `verified` with both engines independently checking load-bearing claims; 27 findings closed across 11 dirty stages (`tautological_check` 13, `mathematica_transliteration` 6, `insufficient_verification` 5, `hardcoded_result` 3); stage 056 verified clean on first read (Series.removeO vs `Limit[pe^2(Omega-pi/2)]` already broke transliteration). Native CAS rewrites used: `Integrate[chiN, {s, 0, l}]` with integer assumption (049), `Solve` on Reals with explicit positivity (051 `xi_(2x)`), `Factor[Together[...]]` canonicalization (051 `Pi_tr`), `Reduce[ForAll[...]]` positivity (060 Onsager check), `Series`/`Coefficient` for small-Pe (053, 058), explicit BVP solve (060 `K_X=0` support). Two batch-wide toolchain patches landed: (a) `expectZero` helper patched in 10 .wl scripts to strip `ConditionalExpression[0, ...]` wrappers that `Solve` introduces under aggressive `$Assumptions`, and (b) stage 051 `Limit` infinity check switched from `pi1 =!= Infinity` to `1/pi1 == 0` to handle Mathematica's non-deterministic `Limit` output. Zero codex iter-2 fixes (matches II.1 and III.1). One stage flagged `material_change: true` (060: failing `sp.solve` replaced with explicit `Csol = a/(exp(a*L)-1)` + Jacobian-aware rescaling, plus `K_X=0` BVP confirmation); downstream Xi_micro consumers in batches III.3+ are still `pending` so no immediate cascade, but second pass should spot-check. See `redteam/batches/batch_III2.md`. |
| P3-07 | Propagate batch III.2 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-22" sentence to `paper/stages/stage_049.tex` through `paper/stages/stage_060.tex` |
| P4-36 | Red-team audit pass on batch III.3 (stages 061-072) | verification | done | 2026-05-22: 12 stages reached `verified` with both engines independently checking load-bearing claims; 27 findings closed across 10 dirty stages (`tautological_check` 14, `mathematica_transliteration` 9, `insufficient_verification` 3, `hardcoded_result` 1); stages 061 and 066 verified clean on first read. Native CAS rewrites used: `sp.solve` parent-action Gaussian elimination (062), `Reduce[... && gphiSq>0, ..., Reals]` (063), explicit profile instantiation for `chi_phi(y)`/`H(y)` (064), parity-driven shell integrals + `gphi` 1/ell scaling (065), `Integrate[sech^2]`/`Integrate[Gaussian^2]` transverse norms (067), `Solve`-derived `Wfail_res`/`Wfail_match` from resonance-corrected premises (068), parameterized `W_match` generator + monotonicity check plus `Cres2Prim`/`Pres = 1/Cres2`/`PresGap` route (069, CHECKPOINT), inlined primitive substitution removing mirrored intermediates (070), `K_m`-pinned `eta`-reconstruction (071), full-vs-leading-order ratio-limit checks with each engine's native `Limit` machinery (072). Zero codex iter-2 fixes (matches II.1, III.1, III.2). One stage (068) flagged `material_change: true`: `Wfail_res`/`Wfail_match` now derived via `Solve` from explicit resonance-corrected premises rather than postulated; the symbolic content of the derived expressions matches the prior postulated forms; downstream III.4+ stages are still `pending` so no immediate cascade. One directive-preordained Blocked (072 F2 `mathematica_transliteration`) resolved by orchestrator as won't-fix-here mitigated by F1's per-engine native-limit ratio checks. Transliteration share back up to 9/12 (vs 6/12 in III.2) reflecting a mirror-heavy cluster on asymptotic leading-order forms in stages 062-068. See `redteam/batches/batch_III3.md`. |
| P3-08 | Propagate batch III.3 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-22" sentence to `paper/stages/stage_061.tex` through `paper/stages/stage_072.tex` |
| P4-37 | Red-team audit pass on batch III.4 (stages 073-084) | verification | done | 2026-05-25: 12 stages reached `verified` with both engines independently checking load-bearing claims; 40 findings closed across all 12 stages (first batch with no clean-on-first-read stages): `tautological_check` 14, `hardcoded_result` 12, `mathematica_transliteration` 7, `insufficient_verification` 7. `hardcoded_result` rose sharply because the Family-1 numerology cluster 075-084 accumulates many literal constants (`alpha_r=10`, `1/20`, `4.06863235...`, `0.927552032...`, `136900`, `1369`). Native CAS rewrites used: parenthesized `Rule` precedence + symbolic `Lambda_ell - L/ell` identity (073), physical substitution chain `m_psi c_s L/hbar -> subs(c_s) -> subs(L/ell)` (074), free-symbol `Delta_0`/`Delta_inf` identity (075), `Integrate[P/rho^2]` + `Solve`-routed `mu_star`/`c_sw` + n=3 fail-check (076), symbolic `1 - alpha_r S(xi_*)^2 = 0` identity (077), symbolic `Sinh`/`Cosh` closed-form coefficients + branch-verdict orderings (078), `Limit[aF1*omega^2, pe -> 0/Infinity]` + `D[omega, pe]` slope check (079), independent `zetaTarget*` numeric paths + four orderings (080), `Solve[zeta == zetaExpr, piTr]` with retrofitted `ConditionalExpression` strip (081), `Solve`-derived `zetaReq` + two derivative assertions (082), `delta0Residual`/`deltaInfResidual`/`omegaResidual` cross-engine identities + monotonicity sign-check (083), cross-route Xi consistency + `FindRoot`/`Limit[zetaPhys, Pe -> Infinity]` (084). One orchestrator-applied mid-batch hot-fix on stage 081 (standard 3-line `ConditionalExpression[e_, _] :> e` strip retrofitted to the `expectZero` helper and to `piOfZeta`/`qq` assignments — same pattern documented in `codex.md` after III.2, but the directive for 081 did not preemptively include it). Two acceptable codex deviations (stage 076: `P = K*rho^n_poly` over the directive's literal — necessary because the directive's form was internally inconsistent; stage 078: removed a spurious `100` factor in the directive — codex caught a math error). One human-resolved docstring discrepancy (stage 076 F2: docstring `(25/4)` vs assertion `25` — cross-referenced paper at `appendices/stage_ledger.tex:221`, `appendix_part03.tex:130`, `stages/stage_077.tex:24`, `stages/stage_082.tex:46` all consistent on `25`; docstring was a stale typo; fix simplified to "delete the `/4` from the docstring" before fix_loop launch — no paper change, no downstream cascade). Zero codex iter-2 fixes. Zero `material_change: true` flags this batch — every derivation-route rewrite left printed symbolic and numeric content byte-identical. Transliteration share 7/12, slightly below III.3's 9/12, because the Family-1 numerology cluster distributed findings more across `hardcoded_result`. See `redteam/batches/batch_III4.md`. |
| P3-09 | Propagate batch III.4 verification status into paper stage cards | provenance | open | optional: add "red-team verified 2026-05-25" sentence to `paper/stages/stage_073.tex` through `paper/stages/stage_084.tex` |
| P4-38 | Red-team v2 paper-grounded re-audit pass on batch I.1 (stages 001-012) | verification | done | 2026-05-25: First batch processed under the v2 auditor, which reads `paper/stages/stage_NNN.tex`, per-stage notes under `notes/stages/`, and part-level appendix BEFORE opening scripts. Added 10th finding category `paper_misalignment` (subtypes: `target_mismatch`, `value_mismatch`, `script_missing_paper_claim`, `paper_missing_script_claim`, `notes_contradicts_script`) routing to user resolution rather than codex. All 12 stages reached `verified`. **5 substantive defects v1 missed**: 7 paper_misalignment items across 5 stages (001 F1+F2 source-coupling/gauge-fix sign disagreements; 006 F1 `Gauge_μ` placeholder in scripts not in paper; 007 F1 paper exports `xi_eff^proj` not verified in scripts; 010 F1+F2 `δu_n` vs `δP_n` distinction + 7 unanchored clusters; 011 F1 5 unanchored clusters with anchors in stages 022-024) plus 3 new script-side tautologies (004 Faraday/Bianchi symbol-substitution, 008 `Integrate[W*Z]` self-cancel, 012 M1 carried-forward self-checks). All 7 paper_misalignment items resolved via a structured **Codex-as-math-authority workflow**: questions markdown file (`redteam/resolutions/batch_I1_paper_alignment.md`) + Codex prompt; Codex read paper/notes/scripts/downstream stages, filled `## Recommendation` blocks with `direction: a|b|c|skip` and citations. Codex answered all 7, skipped 0; user reviewed and concurred on all 7; Codex apply session edited paper + scripts and ran each engine to confirm exit 0. Resolutions chosen: Q1(a) flip script source sign to match paper; Q2(b) flip script gauge-fix sign to match mostly-plus convention per `paper/frontmatter/03_notation_firewall.tex:88`; Q3(c) drop gauge_μ placeholders, reintroduce only via concrete H(w); Q4(a) extend stage 007 scripts with H(w) profile and `xi_eff^proj` checks; Q5(c) paper grows δP_n display, scripts grow δu_n assertions; Q6(a) paper grows 7 cluster anchor paragraphs; Q7(c) trim stage 011 scripts, destinations 022-024 already publish. Two stages flagged `material_change: true` (001 sign flips affect asserted EOMs; 004 Bianchi now exercises real cyclic identity); both downstream effects already absorbed by other stages in the same v2 sweep, so no separate upstream-stale cascade. One new toolchain pitfall documented as #5 in `codex.md` (`Part[]`-on-pattern-parameter inside `Do[Module[...]]` silently drops half the body — Mathematica evaluates Module body during analysis before locals bind; fix is precomputed immediate-valued expressions before any Do/Module scope opens). Stage 007 was v2-re-audited after the H(w) apply since the new code restructured the surrounding sections (line numbers from prior directive's F2/F3 were invalidated); re-audit confirmed F2/F3 didn't recur and only stale_output remained (resolved by output refresh). All 5 paper-alignment apply stages (001/006/007/010/011) had stale `.txt` outputs after codex apply; manually refreshed sympy + mathematica serially per the single-seat Mathematica rule. See `redteam/batches/batch_I1_v2.md` and `redteam/resolutions/batch_I1_paper_alignment.md`. |
| P3-10 | Propagate batch I.1 v2 verification status into paper stage cards | provenance | open | optional: add "red-team v2 verified 2026-05-25" sentence to `paper/stages/stage_001.tex` through `paper/stages/stage_012.tex`; especially stages 001, 004 with `material_change: true` |
| P4-39 | Red-team v2 paper-grounded re-audit pass on batch I.2 (stages 013-023) | verification | done | 2026-05-25: Second batch processed under the v2 auditor. All 11 stages reached `verified` (3 clean: 016, 017, 022). **6 paper_misalignment items + 4 script-side findings** — pattern notably different from I.1: dominated by **(b) trim scripts** because cross-stage cross-check revealed duplication (the EM-projected scripts were file-for-file ports of `notes/em_projected/step_NN_*` master notes, while the later compact paper cards distributed content across multiple stages). Trims removed: stage 013's `δP_2/δP_4`/sieve (destinations: stage 010 owns δP_n per I.1 v2 paper add; stage 014 owns sieve); stage 014's Xi_load/δP_n/Compat (destinations: 013, 010); stage 015's wall-only/Y20/grouped (~half each engine; destination 017 owns lane signature + b=3a + wall-only obstruction); stage 018's one-pole/gate/Xi_1 (destinations 019, 020); stage 020's Y20 block (destinations 010, 017). Sixth item (Q6 stage 021) was inverse direction — script_missing_paper_claim — paper Output enumerated three exports but scripts only asserted 2.5; resolved by **adding** a composed `δD_2^(odd)(ω)` assertion in both engines. After verifier wave caught the initial composed assertion as tautological (used bare `N0` symbol on RHS instead of closed form), Codex remediation pass replaced the RHS with the Section III closed form `(Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² - R²)²`. Stage 023's Schur-complement F2 produced real Codex iteration: directive's `+R_mix` prescription conflicted with existing `+2*g_U*g_W*R_mix` rational numerator; Codex deferred, then remediation pass had Codex read paper's `eq:app-stage023-full-lagrangian` to derive the sign convention from physics — option α chosen (`+R U W` Lagrangian → `-R` off-diagonal in frequency-space spring matrix → existing numerator sign confirmed correct). 4 stages flagged `material_change: true` (013, 014, 015, 018 — all script-side scope reductions removing duplicates; downstream consumers should rebind to destination stages that genuinely own each piece). **User re-direction was load-bearing**: Codex initially recommended (c) acknowledgement for Q1/Q2; user pushed "each step builds on prior" principle; orchestrator cross-check confirmed destination scripts already verified the content → Q1/Q2 flipped to (b) trim before apply. Apply prompt added a **destination-verification guardrail** ("grep destination script to confirm equivalent assertion exists before deleting from source") — worked as designed; no orphan trims. With this batch closed, **the entire linear projected-EM core range (004-021) is now paper-aligned at v2 depth**, plus grouped bookends 022-023. The 4 material_change cascades are contained within the same batch's verified stages (destinations 010/014/017/019/020 are all in batches I.1 v2 or I.2 v2). One stage 014 follow-up displaced from stage 013 F3 (K_1/H_even coefficient guards for literals 1/9, 2/3, -1/27) — optional future fix_loop on stage 014. See `redteam/batches/batch_I2_v2.md`, `redteam/resolutions/batch_I2_paper_alignment.md` (Codex's 6 recommendations with user-revised Q1/Q2), `redteam/resolutions/codex_apply_batch_I2.md`, and `redteam/resolutions/codex_remediation_batch_I2.md` (math-decision remediation pass for 021 and 023). |
| P3-11 | Propagate batch I.2 v2 verification status into paper stage cards | provenance | open | optional: add "red-team v2 verified 2026-05-25" sentence to `paper/stages/stage_013.tex` through `paper/stages/stage_023.tex`; especially stages 013/014/015/018 with `material_change: true` |
| P4-40 | Red-team v2 paper-grounded re-audit pass on batch III.4 (stages 073-084) | verification | done | 2026-05-27: Seventh batch processed under the v2 auditor. All 12 stages reached `verified` (5 clean: 073, 077, 079, 080, 083). **4 substantive paper_misalignment items + 4 banner-relabel paper_misalignment items (audit-flagged) + 14 v2 script-side findings + 8-stage orchestrator-direct banner-relabel sweep**. Substantive items: Q1=(a) Stage 074 `alpha = sqrt(kappa)` value mismatch — paper `128/sqrt(5)`, notes `179/sqrt(5)` (×2 — in 074 and 075 notes), engines both compute `111/sqrt(5)` (since `111² = 12321` exactly); orchestrator-direct fix to `111` in paper:074 line 31, notes:074 line 117, notes:075 line 63, plus added `expect_zero("alpha_ref - 111/sqrt(5)", alpha_ref - 111/sqrt(5))` in both engines. Q2=(a) Stage 075 `Upsilon_w` conversion factor — paper says `117 Theta_w`, notes section 3 says `168 Theta_w`, script uses `100 Theta_w`; paper's own boxed `Theta_fail = 3.626e-4` and the `Xi_F1 = 1369·Upsilon_w = 136900·Theta_w` carry-forward in notes/stages 082/083/084 are mathematically consistent ONLY with 100; orchestrator-direct paper:075 lines 7,24 `117 → 100`, notes:075 lines 108,116 `168 → 100`, notes:075 lines 124-128 `/168 arithmetic → /100`, plus added `alpha_r^2 == 100` assertion lock. Q3=(a) Stage 082 F1+F2 extend scripts to exercise `zeta_phys = Omega_Pe^2(kappa + pi^2/4)/(kappa + y(eta)^2)` closed form and instantiate Family-1 `(eta, kappa) = (37, 12321/5)`; orchestrator-direct added Omega_Pe symbolic definition, `Omega_Pe → pi/2` as `Pe → oo` verified, `y_F1 ≈ 1.52948` computed (via `mpmath.findroot(..., solver='bisect')` because `sp.nsolve` jumped to far roots near `pi/2` — new pitfall #10 candidate), and `zeta_phys(Pe→oo, kappa_F1, y_F1) ≈ 2.4675292` matched stage 084's upstream `zeta_max^(F1)` to `1.77e-13`. Q3-F3 (insufficient_verification): replaced tautological `dR_quad/dzeta_phys + 1 == 0` and `dR_quad/dPi_tr - dzeta_req/dPi_tr == 0` with numerator/denominator factorization `numerator - C_mix*(1-eps_blk) == 0` and `denominator - (C_mix - eps_blk*(2 C_mix - Pi_tr))^2 == 0` (exercises strict-positivity content from notes §4). **Banner-relabel items (audit-flagged):** 078 F2, 081 F1, 084 stale banner, plus 075 F2 banner subcomponent — applied orchestrator-direct. **Orchestrator-direct banner-relabel sweep (post-verifier note):** verifier on stage 076 flagged its stale banner as non-blocking observation; investigation revealed pervasive stale banners across ALL III.4 stages from global renumber (commit `0d09ef6`). Applied 23 banner edits across 11 stages (073, 074, 075, 076, 077, 078, 079, 080, 082, 083 — both .py and .wl; 081 .py only since .wl was already correct; 084 .wl only) aligning all III.4 self-banners with post-renumber numbering. **Other script-side fixes:** 075 F1 replaced tautological "round-trip" and "free-symbol identity" checks with asymptotic-limit `alpha * Delta_inf → 1` (large alpha) and `Delta_0 → 1/2` (small alpha) in both engines (genuinely independent Limit operators); 076 F1 replaced trivial `(2x)^2 = 4 x^2` check with closed-form target comparison exercising the enthalpy-lock `1/4` factor; 076 F2 replaced TODO(provenance) comment with citation to notes section 4; 078 F1 replaced `thetaSuffSym = thetaFailSym × decimal` with explicit Stage-75 symbolic closed form `-(45 cosh(α) + 27√5 sinh(α))/(2500 − 2500 cosh(α))` at `α = 111√5/5`. **Codex stall:** the first codex-chat consultation (Q1) stalled with no session log written; orchestrator killed the stalled processes and pivoted to orchestrator-direct apply — audit + grep evidence on all 3 substantive questions was conclusive without further Codex reasoning. **User-redirection rate 0** (fifth consecutive batch — Codex's stalled-but-orchestrator-recovered recommendation set held up: every direction approved was (a)). **Zero `material_change: true` flags** — paper/notes prose updates align text to the math the scripts already computed correctly, and the stage 082 closed-form pin produces only verification (no new constants). **New pitfall #10 candidate:** SymPy `sp.nsolve` is unstable for roots near a derivative singularity (`y tan y = 37` overshoots to far roots near `pi/2`); use `mpmath.findroot(..., (1.5, 1.55), solver='bisect')` with bracketing. Promote to `codex.md` if recurs. With III.4 v2 closed, **entire range 001-084 is now paper-aligned at v2 depth**. See `redteam/batches/batch_III4_v2.md`, `redteam/resolutions/batch_III4_paper_alignment.md` (3 substantive Codex-evidenced recommendations with orchestrator-direct apply logs). |
| P4-41 | Red-team v2 paper-grounded re-audit pass on batch V.1 (stages 164-175) | verification | done | 2026-05-28: 12 stages reached `verified` (both engines exit 0); 22 findings closed, 0 blocked. Cluster A only (stale `STAGE N-17` script banner across all 12 stages — 11 fixed in per-finding directives, 4 residuals mass-fixed; no Cluster B / Cluster C). One non-cluster paper_misalignment (170 F1, weak-axisymmetric (1,1/2,-1) signature unverified) resolved by user direction (a) = add the missing check to both engines. `mathematica_transliteration` 9/12 — 7 independent routes added, 2 accepted as policy mirrors (169, 175). 1 material_change (170, additive Sec.5 coverage; zero downstream propagation). Three orchestrator catches (166 round-trip scalarization, 175 F1 minimal resolution, 171 bundle-route rework after a `needs_rework` verification). Thirteenth consecutive batch clear of stop-cold. See `redteam/batches/batch_V1_v2.md`. |
| P3-12 | Propagate batch III.4 v2 verification status into paper stage cards | provenance | open | optional: add "red-team v2 verified 2026-05-27" sentence to `paper/stages/stage_073.tex` through `paper/stages/stage_084.tex`; especially stages 074, 075, 082 with paper-side prose updates |
| P4-42 | Red-team IV.x/V.1 orchestrator-direct integrity remediation, batch 1 (stages 108/116/151/170) | verification | done | 2026-05-29: re-verified the four stages whose first-pass fixes had been applied orchestrator-direct (Codex bypassed); all now `verified` with both engines exit 0/0 (116 & 151 in prior commit `e1cdfec`; 170 & 108 this session). **116** F1 mathematica_transliteration → independent eigenvalue solve (root of `Cos[u]==0`, `u=kW·lW`, on `(0,π)`; not a hand-typed `Pi/(2 lW)`), F2 tautological renorm block → labeled prints. **151** F1 insufficient_verification → SymPy rewritten as an EXACT 5-point multi-point cross-check at rational `Pi_star ∈ {1/2,1,3/2,2,5/3}` (symbolic in `r1,r2,A_T,B_T,gprime`), Mathematica is the full all-`Pi_star` symbolic authority; anti-footgun comment forbidding symbolic-`Pi_star` SymPy (it hangs); custom exact `∫₀¹ xⁿe^{ax}dx` integrator (verified sound). **170** F1 tautological_check → Section-5 weak-axisymmetric lane checks rerouted through the DERIVED Section-2 maps `dkappa_from_du2`/`dgamma_from_dP0` (re-typed helpers deleted). **108** F2 mathematica_transliteration → independent non-`Series` route (`chiGenAlt` from raw L-coefficients) confirming `chiGen`; F3 tautological → Class A anchored to literal fingerprint `1+z²/9+4z⁴/81+iz⁵/27`; F4 tautological → Class C anchored to `-Σ0/27` + general round-trip demoted to print. **108 F1 paper_misalignment resolved as direction (a), NO script change** (Claude+Codex concur, see `redteam/codex_reviews/_consult_108_f1.md`). **Deferred paper/doc cleanups surfaced (for the manual paper pass, no script change):** (1) **108** card Checks #2/#3 are block-level over-scoping — cross-reference to the actual verifying stages: **110** (Robin `chi_Q^R=3/(3-rho_R)`), **111** (mixed-pole no-go `kappa_W=-1/9`→`sigma_W=0`), **112** (compensated `chi_Q^hyb`, preservation iff `gamma_W=1/9`). (2) **106** card checks (ii)/(iii) genuinely verified at stages **102/104** — add a card cross-reference. (3) **139** card Checks line is forward-reference boilerplate; the self-matched susceptibility closure is established at stage **140** — defer it there. (4) **148** `redteam/directives/stage_148.md` carries a stale `168π²` that should be `100π²` (script's `100π²` is algebraically correct via `rF1 = √(4107−100π²)/(10π)`; the directive surd is a Codex false-positive) — log as a directive-doc typo to fix when stage 148 is processed. (5) **stage-112 precision caveat** (for 112's own fix loop, NOT a 108 item): "preservation iff `gamma_W=1/9`" holds only for `sigma_W ≠ 0` (at `sigma_W=0`, `chi_B=1` for any `gamma_W`); add the `sigma_W ≠ 0` qualifier / case split when 112 is processed. See `redteam/verifications/stage_{108,116,151,170}.md` and `redteam/codex_reviews/_consult_108_f1.md`. |
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
| 2026-05-22 | Red-team audit completed for batch III.2 (stages 049-060, "Part III.2 -- Tracking, zeta thresholds, asymmetry, boost"). 12 stages reached `verified`; 11 needed codex edits (27 findings: 13 `tautological_check`, 6 `mathematica_transliteration`, 5 `insufficient_verification`, 3 `hardcoded_result`) and 1 verified clean on first read (056). Native CAS rewrites used: `Integrate[chiN, {s, 0, l}]` under integer assumption (049), `Solve` on Reals with explicit positivity constraint (051 `xi_(2x)`), `Factor[Together[...]]` canonicalization (051 `Pi_tr`), forward-map via Stage 047/030 (051), `Reduce[ForAll[...]]` Onsager positivity (060), explicit `K_X=0` BVP solve (060 support). Two batch-wide toolchain patches landed: (a) the `expectZero` helper was patched in 10 `.wl` scripts to strip `ConditionalExpression[0, ...]` wrappers (`Solve` introduces these under aggressive `$Assumptions` and the original `=== 0` check missed them); (b) stage 051's `Limit` infinity assertion switched from `pi1 =!= Infinity` to `1/pi1 == 0` because Mathematica's `Limit` non-deterministically returns either `Infinity` or `Infinity/<positive>` for the same pole. Zero codex iter-2 fixes (matches II.1 and III.1). One stage (060) flagged `material_change: true`: the failing `sp.solve` for `Cnorm` was replaced with the explicit `Csol = a/(exp(a*L)-1)` closed form plus Jacobian-aware rescaling and a `K_X=0` BVP confirmation; downstream Xi_micro consumers in batches III.3+ are still `pending` so no immediate cascade, but the second pass should spot-check. Transliteration share fell from 13/13 (II.1) and 10/12 (III.1) to 6/12 (III.2), reflecting more pure-tautology / hardcoded findings rather than a true policy shift. See `redteam/batches/batch_III2.md`. |
| 2026-05-22 | Red-team audit completed for batch III.3 (stages 061-072, "Part III.3 -- Microclosure, gain thresholds, equilibrium, walls"). 12 stages reached `verified`; 10 needed codex edits (27 findings: 14 `tautological_check`, 9 `mathematica_transliteration`, 3 `insufficient_verification`, 1 `hardcoded_result`) and 2 verified clean on first read (061 and 066). Native CAS rewrites used: `sp.solve` parent-action Gaussian elimination (062), `Reduce[... && gphiSq>0, ..., Reals]` (063), explicit profile instantiation for `chi_phi(y)`/`H(y)` (064), parity-driven shell integrals + `gphi` 1/ell scaling (065), `Integrate[sech^2]`/`Integrate[Gaussian^2]` transverse norms (067), `Solve`-derived `Wfail_res`/`Wfail_match` from explicit resonance-corrected premises (068), parameterized `W_match` generator + monotonicity check plus `Cres2Prim`/`Pres = 1/Cres2`/`PresGap` route (069, CHECKPOINT), inlined `1/Hw -> rhoW/(m*cSw^2)` removing mirrored intermediates `J_1`/`gphi`/`I_1` (070), `K_m`-pinned `eta`-reconstruction (071), full-vs-leading-order ratio-limit checks with the algebraically-equivalent surd `2/Sqrt[5] + 1/(5 + 2 Sqrt[5])` confirming cross-engine independence (072). Zero codex iter-2 fixes (matches II.1, III.1, III.2). One stage (068) flagged `material_change: true`: `Wfail_res`/`Wfail_match` now derived via `Solve` from explicit resonance-corrected premises rather than postulated; the symbolic content of the derived expressions matches the prior postulated forms; downstream III.4+ stages are still `pending` so no immediate cascade. One directive-preordained Blocked (072 F2 `mathematica_transliteration`) resolved by orchestrator as won't-fix-here mitigated by F1's per-engine native-limit ratio checks. Transliteration share back up to 9/12 (vs 6/12 in III.2), reflecting a mirror-heavy cluster on asymptotic leading-order forms in stages 062-068. See `redteam/batches/batch_III3.md`. |
| 2026-05-25 | **First batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch I.1 (stages 001-012, "Part I.1 -- Geometry lift, BdG coupling, projected Maxwell setup") completed 2026-05-25. New auditor reads `paper/stages/stage_NNN.tex`, per-stage notes under `notes/stages/`, and part-level appendix BEFORE opening scripts; added 10th finding category `paper_misalignment` routing to user resolution rather than codex. All 12 stages reached `verified` (4 clean: 002/003/005/009). **5 substantive defects v1 missed**: 7 paper_misalignment items across stages 001/006/007/010/011 + 3 new script-side tautologies (004 Faraday/Bianchi symbol-substitution, 008 `Integrate[W*Z]` self-cancel, 012 M1 carried-forward self-checks). Resolved via Codex-as-math-authority workflow (questions markdown → Codex fills recommendations with file:line citations → user approves → Codex apply edits paper + scripts under explicit scope authorization). Codex answered 7/7 with cross-file evidence the auditor lacked (notation firewall for metric signature, stage 008's H(w) channel, stages 022-024 publishing stage 011's clusters). Two stages flagged `material_change: true` (001 source/gauge sign flips, 004 Bianchi restructure); downstream effects already absorbed by other I.1 stages in the same v2 sweep. One new toolchain pitfall documented as #5 in `codex.md` (`Part[]`-on-pattern-parameter inside `Do[Module[...]]` silently drops half the body). Stage 007 needed a v2 re-audit after the H(w) apply restructured surrounding code; result was only stale_output. All 5 paper-alignment apply stages had stale `.txt` outputs after codex apply; manually refreshed serially. See `redteam/batches/batch_I1_v2.md` and `redteam/resolutions/batch_I1_paper_alignment.md`. |
| 2026-05-25 | Red-team audit completed for batch III.4 (stages 073-084, "Part III.4 -- Family-1 geometry, thresholds, quadrupole"). 12 stages reached `verified`; **first batch with no clean-on-first-read stages** — all 12 needed codex edits (40 findings: 14 `tautological_check`, 12 `hardcoded_result`, 7 `mathematica_transliteration`, 7 `insufficient_verification`). `hardcoded_result` rose sharply from 1 in III.3 to 12 in III.4 because the Family-1 numerology cluster 075-084 packs many literal constants (`alpha_r=10`, `1/20`, `4.06863235...`, `0.927552032...`, `136900`, `1369`, etc.) and the auditors held each to either a derivation or a provenance comment. Native CAS rewrites used: parenthesized `Rule` precedence + symbolic `Lambda_ell - L/ell` identity (073, also fixed a real Mathematica precedence bug where `eta /. (len/ell) -> 37 - 37` parsed as `(len/ell) -> 0`), physical substitution chain `m_psi c_s L/hbar -> subs(c_s) -> subs(L/ell)` (074), free-symbol `Delta_0`/`Delta_inf` identity (075), `Integrate[P/rho^2]` + `Solve`-routed `mu_star`/`c_sw` + non-tautological n=3 fail-check (076), symbolic `1 - alpha_r S(xi_*)^2 = 0` cut-point identity (077), symbolic `Sinh`/`Cosh` closed-form coefficients + three branch-verdict orderings (078), `Limit[aF1*omega^2, pe -> 0/Infinity]` + `D[omega, pe]` slope check returning `(4-Pi)/(2*Pi)` (079), four independent `zetaTarget*` numeric paths + four orderings (080), `Solve[zeta == zetaExpr, piTr]` with retrofitted `ConditionalExpression` strip (081), `Solve`-derived `zetaReq` + two derivative assertions (082), `delta0Residual`/`deltaInfResidual`/`omegaResidual` cross-engine identities + monotonicity sign-check (083), cross-route Xi consistency + `FindRoot`/`Limit[zetaPhys, Pe -> Infinity]` returning `2.467529229456...` matching `zetaMaxF1` to ~14 digits (084). One orchestrator-applied mid-batch hot-fix on stage 081 (standard 3-line `ConditionalExpression[e_, _] :> e` strip retrofitted to the `expectZero` helper and to `piOfZeta`/`qq` assignments — same pattern documented in `codex.md` after III.2, but the directive for 081 did not preemptively include it). Two acceptable codex deviations, both verified necessary (stage 076: `P = K*rho^n_poly` over the directive's literal `P = K*rho^(1+1/n_poly)` — the directive's literal form is internally inconsistent with the `h = m cs^2/4` identity; stage 078: removed the directive's spurious `100` factor on `thetaSuffSym` — codex caught a math error in the directive itself). One human-resolved docstring discrepancy (stage 076 F2: docstring `(25/4)` vs assertion `25` — paper consistently uses `25` at `appendices/stage_ledger.tex:221`, `appendix_part03.tex:130`, `stages/stage_077.tex:24`, `stages/stage_082.tex:46`; docstring was stale typo; fix simplified before fix_loop launch). Zero codex iter-2 fixes (matches II.1, III.1, III.2, III.3 — a remarkable 5-batch streak). Zero `material_change: true` flags this batch — every derivation-route rewrite left printed symbolic and numeric content byte-identical. Transliteration share 7/12 (vs 9/12 in III.3), partly absorbed by the higher `hardcoded_result` share. See `redteam/batches/batch_III4.md`. |
| 2026-05-25 | **Second batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch I.2 (stages 013-023, "Part I.2 -- Maxwell bridge, parent throat action, reduced one-port") completed 2026-05-25. All 11 stages reached `verified` (3 clean: 016, 017, 022). **6 paper_misalignment items + 4 script-side findings** (1 symbol_assumption_error, 3 tautological_check, 4 insufficient_verification, 1 mathematica_transliteration — counting partial F4 K_eta resolution + wall-only auto-trim). Pattern notably different from I.1: **dominated by (b) trim scripts** because cross-stage cross-check revealed duplication — the EM-projected scripts were file-for-file ports of `notes/em_projected/step_NN_*` master notes, and the later compact paper cards distributed content across multiple stages. Trims: stage 013's δP_n/sieve (destinations 010, 014); stage 014's Xi_load/δP_n/Compat (destinations 013, 010); stage 015's wall-only/Y20/grouped (~half each engine; destination 017); stage 018's one-pole/gate/Xi_1 (destinations 019, 020); stage 020's Y20 (destinations 010, 017). Sixth item (Q6 stage 021) was inverse — paper claimed three exports but scripts only asserted 2.5; resolved by **adding** composed `δD_2^(odd)` assertion in both engines (after verifier wave caught the initial composed assertion as tautological, Codex remediation replaced bare `N0` RHS with the Section III closed form `(Ω_A² g_W + R g_A)² / (Ω_A² Ω_W² - R²)²`). Stage 023's Schur-derivation F2 produced real Codex iteration: directive's `+R_mix` prescription conflicted with existing `+2*g_U*g_W*R_mix` numerator; Codex deferred, then remediation pass had Codex read paper's `eq:app-stage023-full-lagrangian` and derive sign from physics — option α chosen (`+R U W` → `-R` off-diagonal → existing numerator sign correct). 4 stages flagged `material_change: true` (013, 014, 015, 018 — script-side scope reductions; downstream consumers should rebind to destination stages). **User re-direction load-bearing**: Codex initially recommended (c) acknowledgement for Q1/Q2; user pushed "each step builds on prior" principle; orchestrator cross-check confirmed destination scripts already verified the content → Q1/Q2 flipped to (b) trim. Apply prompt's **destination-verification guardrail** ("grep destination script to confirm equivalent assertion exists before deleting from source") worked as designed; no orphan trims. With this batch closed, **entire linear projected-EM core range (004-021) is now paper-aligned at v2 depth**, plus grouped bookends 022-023. The 4 material_change cascades are contained within batch-I.1-v2 + batch-I.2-v2 verified stages. One stage 014 follow-up displaced from stage 013 F3 (K_1/H_even coefficient guards for literals 1/9, 2/3, -1/27) — optional future fix_loop. Potential pitfall #6 candidate noted (Mathematica `Dt[..., Constants -> {list}]` may leave unevaluated `Dt[symbol, w, ...]` residuals if `symbol` has lurking `w` dependency not in Constants list — workaround: use ordinary `D` with explicit slot variables) — leave undocumented for now; promote to `codex.md` if recurs. See `redteam/batches/batch_I2_v2.md`, `redteam/resolutions/batch_I2_paper_alignment.md` (Codex's 6 recommendations with user-revised Q1/Q2), `redteam/resolutions/codex_apply_batch_I2.md` (apply session with destination-verification guardrail), and `redteam/resolutions/codex_remediation_batch_I2.md` (math-decision remediation pass for 021 and 023). |
| 2026-05-26 | **Third batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch II.1 (stages 024-036, "Part II.1 -- Overlap isotropy through continuum kernel") completed 2026-05-26. All 13 stages reached `verified` (3 clean: 027, 028, 034). **3 paper_misalignment items + ~15 script-side findings** (8 `insufficient_verification`, 6 `tautological_check`, 4 `mathematica_transliteration`). Pattern differs from I.2's duplication: II.1 was **script-thinness-dominated** — scripts verify the right object but assertions are thin, missing direct anchors. The 3 paper_misalignment items: Q1=(a) cosmetic relabel (stage 029 docstring "Stage 12" → 029, banner "STAGE 012" → 029); Q2=(b) trim (stage 029 `alpha_crit` removed; destination stage 031 owns the refined threshold per `paper/stages/stage_031.tex:43,65`); Q3=(a) paper-side polynomial coefficient fix (stage 035 `eq:app-stage035-F-derivative` bracket polynomial: `206 δ²ξ → 189`, `138 ξ³ → 121`; independent quotient-rule expansion confirmed scripts correct, paper had arithmetic errors). **User-redirection rate 0** this batch — Codex's first-pass recommendations all held up; cross-verification of Q2's destination claim by orchestrator independently confirmed via grep before user gate. **One Codex-as-math-authority remediation** (stage 024 F1/F2): directive prescribed `+rPair` off-diagonal conservative matrix; Codex correctly detected this inverts to `-2 g_U g_W R` mixed term contradicting paper's `Q_r = G_U² Ω_W² + 2 G_U G_W R + G_W² Ω_U²` (`paper/stages/stage_024.tex:108`); remediation pass had Codex read upstream `paper/parts/part01_parent_geometry.tex:956` `+R_l A_l W_l` Lagrangian, derived `-R` off-diagonal, applied with anchor checks. **Symbol-leakage Mathematica hang** (new pitfall #6 candidate): stage 024 Section IV `Table[tripleOverlap[...]]` over 5x5 entries with 6-fold sum + FullSimplify hung >18min CPU after F3/F4 additions introduced global-symbol contamination from earlier sections; resolved by adding `ClearAll[gU, gW, rPair, omegaU, omegaW, mPair, zFromMatrix, nFromMatrix, qRef, hRef, pRef, deltaRef, sRef, ...]` reset at top of Section IV plus memoization of `i4`/`i6` sphere integrals (`i6[args] := i6[args] = Integrate[...]`); runtime dropped from >1080s (killed) to 25.05s. **Zero `material_change: true` flags** this batch — all script-side edits were additions (anchors, eigenvector cross-checks, lane-collapse witnesses) or removals of tautological self-checks; no downstream-visible closed forms changed. With II.1 v2 closed, **entire range 001-036 is now paper-aligned at v2 depth**. Cosmetic follow-ups (non-blocking): stage 024 has an empty Section II.1 subbanner after F3 deletions; stage 026 Mathematica final banner still reads "Stage 9 Mathematica audit passed." (legacy renumbering); stage 032 `expectZero` label at `.wl:172` is malformed (SymPy+Mathematica labels concatenated). See `redteam/batches/batch_II1_v2.md`, `redteam/resolutions/batch_II1_paper_alignment.md` (3 Codex recommendations), `redteam/resolutions/codex_apply_batch_II1.md` (apply session), and `redteam/resolutions/codex_remediation_batch_II1.md` (stage 024 sign convention + Section IV performance fix). |
| 2026-05-26 | **Fifth batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch III.2 (stages 049-060, "Part III.2 -- Tracking, zeta thresholds, asymmetry, boost") completed 2026-05-26. All 12 stages reached `verified` (5 clean: 051 checkpoint, 052, 053, 056, 060). **2 paper_misalignment items + 14 script-side findings** (8 `insufficient_verification`, 3 `tautological_check`, 2 `mathematica_transliteration`, 1 `symbol_assumption_error`). The 2 user-gate items: Q1=(a) Stage 050 paper-card extension with fifth boxed equation `S_n^{twin}(x;ε) < S_n^{max}(ε) := 1 + (1-ε)/((2n+1)^2 - ε)` (label `eq:app-stage050-Sn-max`) — orchestrator independently confirmed no downstream consumers in 051-072 so option (b) "relocate" had no natural target; Q2=(a) Stage 057 script-side Pe-monotonicity numerical sweep added (`Pe ∈ {1/10, 1/2, 1, 2, 5, 10}` at `(kappa=1, y=π/4)`) — orchestrator destination check confirmed Stage 056's scripts only verify the covariance identity `dOmega_Pe/dPe = Cov(chi_0,s)/I_W`, NOT the sign, so option (b) "rely on Stage 056 carry-forward" was unsound at the script level. **User-redirection rate 0** this batch (third consecutive after II.1 and III.1 also 0). **One orchestrator hot-fix on stage 058**: Codex iter2's F2 prescription added a full `sp.dsolve` / `DSolve` symbolic BVP solve + boundary-condition `solve` + `simplify(phi_drop - Delta)` to verify the Green-kernel construction. The sympy version hung at 100% CPU for 7+ hours before being killed. Orchestrator replaced both engines' dsolve blocks with the equivalent kernel-integral identity: sympy uses a numerical sweep `Delta = integral(K * Sigma_Pe)` on 4 concrete `(α,η,Pe)` tuples; Mathematica relies on its pre-existing line 84 check `delta independent integral matches combination form` which compares the Green-function integral side against the `Ic`/`Is` combination closed form. Also added `Pe == α` singularity guards in sympy's monotonicity/IVT sweeps (Delta has removable 0/0 at `Pe = α`; `subs()` doesn't take the limit). **New pitfall #8 candidate** (heavy BVP `dsolve` via symbolic engines is not worth the cost): the natural equivalent is `Delta = ∫ K · Σ` which can be verified symbolically (when sympy can integrate the form) or numerically on a concrete parameter grid; the dsolve path adds nothing the integral path doesn't cover. Promote to `codex.md` if recurs. Stage 060 (v1 `material_change: true` carried forward) returned **clean (0 findings)** under v2 audit — the v1 gain definition `Xi_micro = Λ²L²/(Θ T_X)` is sound at v2 depth. **Zero `material_change: true` flags** this batch — Stage 050's paper card extension does not change any script export value, and Stage 057's added Pe-monotonicity sweep is an addition not a value change. With III.2 v2 closed, **entire range 001-060 is now paper-aligned at v2 depth**. Cosmetic follow-ups (non-blocking): stages 049, 050, 052, 053, 054, 055, 056, 057, 058, 059 all carry legacy banners ("STAGE 32" through "STAGE 042") from pre-renumbering; stage 050 mathematica F1 label has cosmetic `= 0` doubling; mathematica 058 emits benign `Power::infy`/`N::meprec`/`Limit::alimv` warnings at the `Pe = α` singular points but all assertions pass. See `redteam/batches/batch_III2_v2.md`, `redteam/resolutions/batch_III2_paper_alignment.md` (2 Codex recommendations with apply logs), `redteam/resolutions/codex_prompt_batch_III2.md`, `redteam/resolutions/codex_apply_batch_III2.md`. |
| 2026-05-26 | **Fourth batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch III.1 (stages 037-048, "Part III.1 -- Continuum kernel, generalized branch, rank-2") completed 2026-05-26. All 12 stages reached `verified` (4 clean: 038, 040, 041, 044). **3 paper_misalignment items + 1 insufficient_verification user-gate item** + ~10 script-side findings (5 `insufficient_verification`, 4 `tautological_check`, 2 `mathematica_transliteration`). The 4 user-gate items: Q1=(a) Mathematica stage 043 D_phi sign fix (`Det[{{kappa0,kappa1},{y0,y1}}]` row swap per stage 039 D_dir kappa-first convention; no downstream D_phi consumers in 044-048); Q2=(a) Stage 044 `F_cont` residual import into stage 045 (`scripts/moving_throat_pde_stage044_continuum_selected_rank2_sympy_audit.py:82-90,140-146` destination-verified); Q3=(b) stage 045 script + notes label relabel from "Stage 28/028" to "Stage 045" (scripts: docstring/banner/final-print; notes: line 232 heading); Q4=(a) stage 046 notes coefficient typos fix (P_R: 230→162×2, P_1: 248→180, 230→162, P_2: 237→220 — scripts unchanged because both engines independently verify the correct script-side values per `mathematica/output/moving_throat_pde_stage046_..._txt:14,20`). **User-redirection rate 0** this batch (second consecutive after II.1 also 0). One Codex iter2 remediation pass for stage 043 F2-Insertion2: directive prescribed `subs(sigma0=0, rho0=1)` symbolically but Mathematica primary symbols are primitive couplings `gB, gU, gS, gR, gW, kU`; iter2 fix used `gS → 0` for sigma_0=0 and `gW → gU·gR/kU` for rho_0=1, with expected numeric RHS `+1/4, -1/4, 0` per the closed-form `R_phi - R_U = deltaU·(rho_0 - sigma_0)/[(1+deltaU)(1+rho_0)(1+sigma_0)]`. **New pitfall #7 candidate** (primitive-vs-derived substitution): when directive prescribes `subs(<derived>, value)` but script uses primitive symbols underneath, must lift to the primitive symbols that realize the derived value. **One stage with `material_change: true`** — stage 045 F3 imported Stage-044's `F_cont` residual into 045's verification path; the export `F_tr` value is byte-identical to before the import, so this is **structural-only material_change**. Downstream cascade NOT triggered (running `$RT mark-stale-downstream 045` would have demoted 046-253 unnecessarily). With III.1 v2 closed, **entire range 001-048 is now paper-aligned at v2 depth**. Cosmetic follow-ups (non-blocking): stages 037/047 have legacy "STAGE 20/30" banners; stage 043 directive has two `## Applied: F2` blocks (one from Q1 apply session, one from fix_loop). See `redteam/batches/batch_III1_v2.md`, `redteam/resolutions/batch_III1_paper_alignment.md`, `redteam/resolutions/codex_apply_batch_III1.md`, `redteam/resolutions/codex_remediation_043_iter2.md`. |
| 2026-05-27 | **Seventh batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch III.4 (stages 073-084, "Part III.4 -- Family-1 geometry, thresholds, quadrupole") completed 2026-05-27. All 12 stages reached `verified` (5 clean: 073, 077, 079, 080, 083). **4 substantive paper_misalignment items + 4 audit-flagged banner relabels + 8-stage orchestrator-direct banner-relabel sweep** when investigation revealed pervasive stale banners across ALL III.4 stages from the global renumber. Substantive: Q1=(a) Stage 074 alpha = sqrt(kappa) = 128/sqrt(5) (paper) vs 179/sqrt(5) (notes, ×2 — 074 and 075) vs engines' 111/sqrt(5) — both engines independently compute `sqrt(12321/5) = 111/sqrt(5)` (111² = 12321 exactly), and the numerical value 49.6407091 cited in both notes files matches 111/sqrt(5) (not 179 or 128); orchestrator-direct fix to paper:074:31, notes:074:117, notes:075:63, plus added `alpha_ref - 111/sqrt(5) == 0` assertion in both engines. Q2=(a) Stage 075 Upsilon_w = 117 Theta_w (paper Inputs/body) vs 168 Theta_w (notes §3) vs 100 Theta_w (script's alpha_r=10 → alpha_r²=100); paper's own boxed Theta_fail = 3.626e-4 Pe_req and the carry-forward Xi_F1 = 1369·Upsilon_w = 136900·Theta_w in notes/stages 082/083/084 are mathematically consistent ONLY with 100; orchestrator-direct fix to paper:075:7,24, notes:075:108,116,124-128 (arithmetic too), plus added `alpha_r^2 == 100` lock. Q3=(a) Stage 082 F1+F2 extend scripts to exercise paper-claim `zeta_phys = Omega_Pe^2 (kappa + pi^2/4)/(kappa + y(eta)^2)` closed form and instantiate Family-1 (eta=37, kappa=12321/5); orchestrator-direct added symbolic Omega_Pe definition, verified Omega_Pe → pi/2 as Pe → oo, computed y_F1 ≈ 1.52948 (smallest positive root of y tan y = 37 in (0, pi/2)) via `mpmath.findroot(..., solver='bisect')` since `sp.nsolve` jumps to far roots near pi/2 (new pitfall #10 candidate), and matched zeta_phys(Pe→oo, kappa_F1, y_F1) ≈ 2.4675292 against stage 084's upstream zeta_max^(F1) constant to 1.77e-13. F3 replaced tautological derivative checks with numerator/denominator factorization exercising notes §4's strict-positivity content. **Banner relabels (4 audit-flagged + 8 orchestrator-direct sweep)**: 074, 078 (.py+.wl), 081 (.py), 084 (.wl) flagged by auditors; orchestrator-direct sweep also fixed 073, 075, 076, 077, 079, 080, 082, 083 banners (23 total label edits across 11 stages). **Codex stalled** mid-Q1 consultation (no session log written); orchestrator killed processes and pivoted to direct apply since audit + grep evidence on all 3 substantive questions was conclusive. **User-redirection rate 0** (fifth consecutive batch). **Zero `material_change: true` flags** — paper/notes prose updates align text to math the scripts already computed correctly; stage 082 closed-form pin produces only verification, no new constants. **New pitfall #10 candidate**: SymPy `nsolve` is unstable near singularities of the derivative (Newton overshoots `tan y` near pi/2); use bracketing solver (`mpmath.findroot(..., solver='bisect')`) instead. With III.4 v2 closed, **entire range 001-084 is now paper-aligned at v2 depth**. See `redteam/batches/batch_III4_v2.md`, `redteam/resolutions/batch_III4_paper_alignment.md`. |
| 2026-05-27 | **Ninth batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch IV.1 (stages 091-102, "Part IV.1 -- Grouped p2 geometry, decoupling, contamination") completed 2026-05-27. All 12 stages reached `verified` (1 clean: 093 — status-only Mathematica-only mirror per MANIFEST `is_status_only_candidate`). Checkpoint stage 096 (geometry-lane check verdict) passed at the higher-bar standard. **27 findings closed + 1 `blocked_legitimate` (100 F4 design-level Mathematica transliteration rewrite, deferred). 10 paper_misalignment items resolved per 3 user-gate questions, plus a 23-site orchestrator-direct banner-relabel sweep across all 12 IV.1 stages.** Substantive: **Q1=(a) Cluster A (091/095/097/099)** — paper card `\stagefield{Checks}` items "l=0/l=2 orthogonality + static-limit eps_2=eps_4=0 -> c_pole=1/4 + minimal-module hypothesis" were flagged by auditors as paper_misalignment in four stages. Resolution: SymPy docstring carry-forward annotations naming the upstream verifying stage(s) — orthogonality at 094 (15 angular integrals + Laplace eigenvalue); static-limit at 091/092/094/096; minimal-module at Part III chain 088-090. No new script-side assertions; no paper edit. Mirrors III.5 087 F1 consolidation pattern. **Q2=(c) Cluster B (100/101)** — stage 100 headline closure `mhat_0^2 chi_Q N_Q = 1` was being verified tautologically (A3/A8 forced by `K_n = N_Q * K_n_target` construction); stage 100/101 paper cards also listed DtN-fingerprint (Check 3, anchored at stage 097) and higher-odd-term placement (Check 2, anchored at stage 102) that scripts didn't exercise. Resolution: stage 100 substantive closure — IMPOSE `mhat_0^2 * Gamma_5 = Gamma_5_target` as the observable condition on script-derived Gamma_5(K_0, chi_Q, Omega), derive `closure_ratio - (mhat_0^2 chi_Q N_Q - 1) = 0`. Non-tautological — a wrong factor in Gamma_5 series derivation would break the closure ratio. F2/F3 close automatically; F5 (chiQ symbol_assumption_error) also fixed (declared real, not positive — chi_Q is pinned to 1 by upstream DtN comparison at stage 097, not script-side constrained here). F4 design-level Mathematica rewrite remains `blocked_legitimate`. Stage 101 docstring carry-forwards added naming stages 102/097 as upstream anchors. **Q3=(a) Cluster C (banner sweep, all 12)** — every paper card section title uses stale stage numbers (091→108, ..., 102→119) from a previous renumber, and every script banner uses stale "STAGE 74" through "STAGE 085". Resolution: script-side 23-site banner sweep applied (12 .py + 12 .wl minus 093 no-py = 23 sites); paper card section titles deferred to PAPER_CLEANUP_TRACKER for a future paper-side pass. Also docstring `Stage 75 → Stage 092` reference fixes in 095/096 SymPy. **Codex bypassed entirely** (seventh consecutive zero-redirection batch); orchestrator-direct math-authority workflow held cleanly because audit + grep + resolution evidence was conclusive on all three clusters. **One material_change: true at stage 100** (verification surface strengthened from tautological cross-check to substantive observable-condition closure; no derived numeric value changed — `mhat_0^2 chi_Q N_Q = 1` is the same Output as before; downstream stages > 100 not marked `upstream_stale`). All other stages zero material_change. **PASS-line counting confirmed expected substantive coverage** across all 12 stages (e.g., 094 Mathematica 34 PASS = 30 orthogonality + 1 Y00 norm + 3 static-limit; 097 Mathematica 9 PASS = 2 series-equiv + 1 Gamma5 closed form + 1 geometric target + 1 Gamma5_target + 4 R_i). No new pitfall candidates from IV.1 — all orchestrator-direct edits applied first-attempt clean. With IV.1 closed, **entire range 001-102 is now paper-aligned at v2 depth**. **Deferred paper-side cleanup** (held for future paper-cleanup pass): 12 paper card section titles in IV.1 (`paper/stages/stage_{091..102}.tex` line 1) all carry stale Stage 108-119 numbers from a previous renumber while `\label{stage:0NN}` and body text correctly use 091-102; rewriting these is a paper-only edit and was explicitly deferred per Cluster C direction (a). See `redteam/batches/batch_IV1_v2.md`, `redteam/resolutions/batch_IV1_paper_alignment.md`. |
| 2026-05-28 | **Fourteenth batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch IV.6 (stages 151-163, "Part IV.6 -- Correction, coevolution, traction, off-family") completed 2026-05-28. All 13 stages reached `verified` (3 clean: 153 status-only consolidation card per MANIFEST `is_status_only_candidate`, 159 and 162 substantive-clean). **No checkpoints** in IV.6 range. **19 findings closed across 10 dirty stages (151/152/154/155/156/157/158/160/161/163).** 19 resolved + 0 blocked. **One material_change** (stage 151 SymPy rewrite from symbolic-integration to mpmath numerical because `sp.integrate(Pi_star·exp(-Pi_star·x)·cos(πx/2)·R(x), (x,0,1))` with free `Pi_star` hung >30 min CPU; downstream impact zero — script verifies algebraic identities not numeric carry-forwards). **3 user-gate clusters resolved (all "Recommended")**: **Cluster A (mass renumbering, 20 edits)** — 18 banner edits across 9 stages × 2 engines at −17 offset (154/155/156/158/159/160/161/162/163; 151/152 use title-only banner convention with no offset; 157 banner was correct) + 2 notes H1 edits at −85 offset (159 `Stage 244→159`, 160 `Stage 245→160`); applied via `/tmp/iv6_mass_renumber.py`. **Cluster B (body-text forward-stage citation re-attribution, 53 edits across 13 notes files)** — three offsets surfaced (content-verified): −51 for 188-199 (IV.4 references, e.g., 188-189→137-138); −85 for 239-248 (IV.6 self-references, e.g., 244→159, 246→161, 248→163 — **new offset specific to IV.6**, disambiguated from −102 by content cross-check, e.g., stage 158's "Stage 241, exact co-evolving compensated Family-1 point" maps to stage 156 renormalized canonical branch via −85, not to stage 139 family1_actual_mouth_gains via −102); −102 for 221 (→119 parent compensation family, IV.3) and 249-250 (→147-148, IV.5). Applied via `/tmp/iv6_reattribute.py`. **Cluster C (158 paper-card Checks downgrade)** — items 2 (even-preservation constraints) and 3 (tangent motion δ⊥=0) describe the broader transport program but are verified DOWNSTREAM in stages 159, 162, 163. Edit applied to `paper/stages/stage_158.tex:23-24`: items 2 and 3 rewritten as forward-carry citations of `\ref{stage:159}` (even-preservation) and `\ref{stage:162}` / `\ref{stage:163}` (tangent motion). **First forward (downstream) carry-forward in v2 batches** — IV.4 stage 134 and IV.5 stage 144 both carried items UPSTREAM. **Mechanical findings (per-stage substance)**: 151 F1 → SymPy mpmath numerical integration at `Pi_star=1.50882951349316, r1=1.7, r2=-0.9` with ~40 dps tolerances on centering, moment shifts, and bias/traction retunings (replaces tautological definition-vs-self check); F2 subsumed by F1; F3 Mathematica derives `deltaSigma` via `Normal[Series[Exp[-Phi[x]]/Z, {epsilon, 0, 1}]]` independent of SymPy hand-form, cross-checks `deltaSigma + SigmaStar·(R-<R>) == 0`. 152 F1 → 4 anchored `expect_close` lines for deltaPi_act, deltaTm_act, lambda_eff^(Pi), lambda_eff^(T). 154 F1 → Mathematica `rShiftSeries = Normal[Series[rFun /. g -> gStar + dg, {dg, 0, 2}]]` and single-epsilon parameterization for `piLin`; orchestrator-rework after directive's multivariate `Series` retained cross-products and `dPi identity` failed (residual `-(dS·dSigma0·rStar) - dR·(dS·(dSigma0+sigma0)+dSigma0·sStar)`). 155 F1 → 4 numeric `assert` lines for g_fp/S_fp/R_fp/Pi_fp at 1e-12 tol; F2 banner via Cluster A. 156 F1 → 4 numeric asserts (Sigma0_can, S_can, Pi_can, T_hat_can) at 1e-10/1e-11 tols. 157 F1 → fresh-symbol `Solve[{dCsym - 9·sigmaStar·dKsym == 0, 5·dCsym - 72·sigmaStar·dKsym == 0}, {dCsym, dKsym}]` replaces `-16·sigmaStar·(dR /. dg -> gp·dr)` projector in both engines; F2 substantive rewrite via Solve (option B); F3 stale_output resolved by refresh. 158 F1 paper-card via Cluster C (scripts untouched); F2 delete `delta Ms law` definition + assert (tautological by construction); F3 add composed `delta Mq law` and `delta Pi law` symbolic identities in both engines. 160 F1 → Mathematica chain-rule total-differential `dKappa0/(1+rStar) - (kappa0Canon/(1+rStar)²)·deltaR` replaces Series+Coefficient recipe; full symbol rename for binding hygiene. 161 F1 → SymPy `epsg_exact = 9·gamma0/(1+rc) - 1`, take `sp.diff` for `depsg_direct`, substitute `gamma0_sym = (1+rc)/9` first to get `depsg_branch`, then substitute `dgamma0 -> (1+rc)·dln_gamma0/9` (orchestrator-catch: directive's prescribed `depsg_direct.subs(dgamma0, ...)` left free `gamma0` and yielded non-zero residual `drc·(rc+1-9·gamma0)/(rc+1)²` — fixed by using `depsg_branch` instead); F2 linearize exact `BW` via `BW.subs(eps_kappa: eps·deps_k, eps_gamma: eps·deps_g, rc: rc_star+eps·drc); dBW = sp.diff(BW_pert, eps).subs(eps, 0)`; F3 banner via Cluster A. 163 F1 → Mathematica `gPrimeImplicit = -D[fComp,r]/D[fComp,g]` (implicit-function ratio) + chain-rule `Series` on `Log` of parent expressions to derive `deltaR`/`deltaG`/`deltaPerp` independently of SymPy hand-form. **Codex bypassed entirely** (twelfth consecutive zero-redirection batch); orchestrator-direct math-authority workflow held cleanly across all three clusters. **Three orchestrator catches in the rework loop**: (1) stage 151 SymPy `sp.integrate` hang (>30 min, killed) → mpmath rewrite; (2) stage 154 Mathematica multivariate `Series` retains cross-products → switched to single-epsilon parameterization `piExprEps /. epsLin -> 1`; (3) stage 161 directive variable-substitution typo (used `depsg_direct` instead of `depsg_branch`) → fixed in both engines. **First v2 batch where the auditor-prescribed engine approach was infeasible at the engine level and the orchestrator had to redesign the verification.** With IV.6 closed, **entire range 001-163 is now paper-aligned at v2 depth**. **No new paper card display titles deferred** in IV.6 — previously deferred IV.1/IV.2 paper titles remain as P3-13. See `redteam/batches/batch_IV6_v2.md`, `redteam/resolutions/batch_IV6_paper_alignment.md`. |
| 2026-05-27 | **Thirteenth batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch IV.5 (stages 139-150, "Part IV.5 -- Susceptibility, branch, defect transport") completed 2026-05-27. All 12 stages reached `verified` (3 clean: 141, 145, 149 — all status-only consolidation cards per MANIFEST `is_status_only_candidate`). **No checkpoints** in IV.5 range. **31 findings closed across 9 dirty stages (139/140/142/143/144/146/147/148/150).** 30 resolved + 1 `blocked_legitimate` (144 F4 transliteration policy — user-accepted: symbolic forms of `r`, `g_∓^F1`, `g_Π`, `s_Q`, `r_Q` are upstream imports, second-engine value is the independent numerical root-finder). **Zero material_change** in IV.5 — every fix was renumbering, notes re-attribution, paper-card downgrade, or script-substance addition; no derived numerical constants moved. **3 user-gate clusters resolved (all "Recommended")**, mirroring IV.4: **Cluster A (mass renumbering, 21 edits)** — 11 `.wl` banner edits (off by −17: 139→122, 140→123, 142→125 + LEDGER, 143→126 + LEDGER, 144→127 + LEDGER, 146→129, 147→130, 148→131) + 6 `.py` banner sites (142, 143, 144 with LEDGER variants) + 4 notes/ H1 lines (146→248, 147→249, 148→250, 149→251; off by +102); applied via `/tmp/iv5_mass_renumber.py`. **Cluster B (body-text forward-stage citation re-attribution, 22 edits)** — 11 of 12 notes contain inline citations to pre-renumber stages 188-251; two offsets surfaced: −51 for 188-199 range (matches IV.4 body offset; 188-189 → 137-138 mouth-gain formulas; 188-191 → 137-140 mouth-gain block; 197-199 → 146-148 positive-deformation block) and −102 for 220-251 range (matches IV.5 notes-H1 offset; e.g., 220 → 118 shell coupling, 228 → 126 broadening fraction, 232 → 130 mouth-bias map, 237 → 135 outlet-consistent closure, 240→138, 241→139, 242→140, 244→142, 247→145, 248→146, 249→147, 250→148, 251→149); applied via `/tmp/iv5_reattribute.py`. **Cluster C (144 paper-card Checks downgrade)** — items (i) outlet-consistency (subject of stage 135) and (ii) self-matched susceptibility closure (subject of stage 140) were not exercised by either 144 script. Edit applied to `paper/stages/stage_144.tex:21-25`: items (i) and (ii) rewritten as carry-forward citations of `\ref{stage:135}` and `\ref{stage:140}` respectively; item (iii) (numerical fixed-point recording) unchanged. Mirrors IV.4's stage 134 downgrade pattern. **Mechanical findings (per-stage substance)**: 139 F1 → 9 anchored asserts each engine (r_F1, R_q^nat, M_s/M_q natural+compensated, outlet consistency, R_q^comp closed form); 139 F2 → `Solve[(gc-rF)^2 == (1+rF^2)/4 && gc<rF, gc, Reals]` independent gMinus derivation in Mathematica; 140 F2 → 3 numerics asserts in both engines (That_nat, That_comp, fractional enhancement); 142 F1 → R_q(Pi_*)−1/4 numeric anchor (replacing R_q(g_-)−1/4 tautology); 142 F2 → `expect_close` helper + 5 canonical-point anchors; 142 F3 → series-vs-closed-form + r_F1 algebraic identity in Mathematica; 143 F1 → 8 paper-deliverable asserts + helpers; 143 F2 → replaced hardcoded `sInf=1` with `Limit[]` over dynamical `rQ`/`sigma0`/`that`; 143 F3 → `Reduce[num>0, piM, Reals]` independent positivity; 144 F2 → 7 numerical-target asserts (Pi_*, Sigma_0, That, Pi_match, That_match, upper g_+^F1>1, lower bracket); 146 F1 → integral-form affine-law (eps-sample numeric fallback both engines: SymPy eps∈{1/10,1/2} 1e-15 tol, Mathematica same with 1e-6 tol commensurate with Integrate-with-numeric-pStar precision-9 residuals); 146 F2 → symbolic moment direct-formula checks with `.has(sp.Integral)` Pi-sample numeric fallback at Pi∈{7/10, 11/10, 17/10, 23/10}; 147 F1+F2+F3 → 3 paper-quoted A_T/B_T anchors + chain-rule consistency (`AT_30 = sp.N(AT, 30)` for full-precision compare) + centered-kernel structure check (`Chop[..., 10^-25]` wrap for precision-28 near-zero) + moment-stability resubstitution drift checks in both engines; 148 F1 → replaced hardcoded `xi_star = 0.183918405511538` with symbolic closed form `(-37√3 - 5π² + 2√(4107 - 100π²))/(5(8 - π²))` — **note**: this required correcting a directive-prescribed typo (auditor copied stage 148 notes' `4107 - 168π²` which numerically gave ~1.547; orchestrator caught via Mathematica's `FullSimplify`-yielded-symbolic-Sqrt-residual signature, cross-referenced stage 126 upstream notes which had the correct `4107 - 100π²` ≈ 0.184, and corrected both engines plus the stage 148 notes typo); 148 F3 → restructured Mathematica `dT` derivation via intermediate `dSigmaOfDeltas`/`dTOfDeltas` helpers; 150 F1 → replaced tautological `S_q = simplify(diff(T_q,x).subs(x,0))` with hand-derived `Aq*k - Cq*Pi`. **Codex bypassed entirely** (eleventh consecutive zero-redirection batch); orchestrator-direct math-authority workflow held cleanly across all three clusters. **Pitfall #13 re-confirmed** at stage 139 — Mathematica comment `(* Pi_* and S_q(Pi_*) imported from Stage 134. *)` parsed `Pi_*)` as comment-terminator; ASCII-safe rewrite required. Four orchestrator catches in the rework loop: (1) stage 148 directive-typo; (2) stage 139 pitfall #13 recurrence; (3) stage 142 SymPy tolerance `1e-20` → `1e-15` for nsolve's actual `~1.95e-18` precision; (4) stage 147 SymPy `sp.N(AT)` default-15-digit truncation → explicit `sp.N(AT, 30)`. Plus stage 146 needed two SymPy numeric-sample fallbacks (Pi-sample for F2 `Integrate` unevaluation, eps-sample for F1 `simplify` non-reduction). With IV.5 closed, **entire range 001-150 is now paper-aligned at v2 depth**. **No new paper card display titles deferred** in IV.5 — previously deferred IV.1/IV.2 paper titles remain as P3-13. See `redteam/batches/batch_IV5_v2.md`, `redteam/resolutions/batch_IV5_paper_alignment.md`. |
| 2026-05-27 | **Twelfth batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch IV.4 (stages 127-138, "Part IV.4 -- Penetration, mouth boundary, fixedpoint") completed 2026-05-27. All 12 stages reached `verified` (2 clean: 128, 129 — 128 is status-only consolidation card per MANIFEST `is_status_only_candidate`; 129 had only a cosmetic banner note not raised as a formal finding). **No checkpoints** in IV.4 range. **22 findings closed across 10 dirty stages (127/130/131/132/133/134/135/136/137/138).** 21 resolved + 1 `blocked_legitimate` (137 F3 matrix-Schur — directive itself acknowledged the prescribed block was tautological without sign-convention bookkeeping; F1's anchors cover the substantive claim). **Zero material_change** in IV.4 — every fix was renumbering, notes re-attribution, paper-card downgrade, or script-substance addition; no derived numerical constants moved. **3 user-gate clusters resolved (all "Recommended")**: **Cluster A (mass renumbering)** — 9 `.wl` banner edits (`STAGE N` where N = actual − 17: 127→110, 129→112, 130→113, 131→114, 133→116, 134→117, 135→118, 137→120, 138→121) + 4 `.py` banner sites (127, 133×2, 134, 135) + 10 notes/ H1 lines (127-136, where H1 said `Stage M` with M = actual + 102: 229→127, …, 238→136); 19 mechanical text edits via `/tmp/iv4_mass_renumber.py`. Files 137 and 138 H1 were already correct. **Cluster B (status-only carry-forward re-attribution)** — stages 132 and 136 notes attributed load-bearing constants to *downstream* stages (132 cited 180-182; 136 cited 184-186). The constants actually originate within IV.4 itself (Π_* at 131, g_Π at 130, σ_Π at 129, coupled fixed-point law at 133, F1 first explicit branch at 134, Σ_m^* at 135). Resolution: notes 132 line 6 `Stages 180–182` → `Stages 129–131`; notes 136 line 6 `Stages 184–186` → `Stages 133–135`. **Cluster C (134 paper-card Checks downgrade)** — items 1 (outlet consistency, subject of stage 135) and 2 (susceptibility closure, runs through 137-138) were not exercised by either 134 script. Edit applied to `paper/stages/stage_134.tex:21-25`: items 1 and 2 rewritten as carry-forward citations of `\ref{stage:135}` and `\ref{stage:137}` respectively; item 3 (numerical fixed-point recording) unchanged. Resolves 134's F3 paper_misalignment without script duplication. **Mechanical findings (per-stage substance)**: 127 F2 mathematica_transliteration → independent `Integrate[...]` derivations for `g_slab`/`g_exp` plus closed-form-match `pass[]` lines; 130 F1+F2 → 6-point `dg/dPi > 0` monotonicity sweep at Π = 1/10, 1/2, 1, 1.5088, 3, 10 + `gPi == gPi_boxed` closed-form-vs-paper-form assertions in both engines; 131 F1+F2+F4 → 4 anchored SymPy asserts (Π_*, slope at Π_*, threshold identity, lower-branch discrimination) + Mathematica banner fix + tautological `expectApprox["Pi_* compensation point", ...]` replaced with the same 4 anchored checks (rewritten with ASCII labels — `piStar`, `slope at piStar` — because Mathematica's parser chokes on the `_*)` substring in `Pi_*)` near a comment terminator; pitfall #13 candidate) + closed-form `g_minus = (2√(4107-100π²) - 37√3)/(20π)` derivation with literal-vs-closed-form PASS check (Anchor 3's `Simplify === 0` test required `Chop[..., 10^-30]` because FindRoot's numeric `piStar` produces precision-79 near-zero residue rather than exact symbolic zero); 133 F1 → hand-ansatz block replaced with `DSolveValue[{ODE, BCs}, uFun[x], x]` so Mathematica derives `u(x)` from PDE+BCs independently; 134 F1+F2+F4 → SymPy 4 new substantive asserts (S_shell=1, S_q at Pi=1/2/1/2 vs externally-verified mpmath targets, S_q(Pi_*) value, gain line coefficients) + Mathematica 3 `expectClose` lines against same targets ((**orchestrator directive-target catch**: auditor agent fabricated literal targets for `S_q(1/2)`, `S_q(1)`, `S_q(2)` that disagreed with mpmath-computed exact values by orders of magnitude; orchestrator recomputed at 50 digits and substituted 0.608336415687717…, 0.633127670034487…, 0.681366857005321…)); 135 F1+F2 → SymPy 5 new substantive asserts (outlet substitution identity, `0 < S_q(Π_*) < 1`, `Sigma_m^* ≈ 0.451485`, `M_s^* ≈ 1.80594`, `M_q^* ≈ -0.451485`, mixed-lane correction `M_q^* * S_q(Π_*) ≈ -0.297112`) replacing the tautological closure residual `raise`; 137 F1+F2 → 3 new anchored asserts in both engines (paper-card closed-form `M_s`/`M_q` via independent route, Schur static limit `delta_Lambda_core(z→0) → rho_c - sigma_c` with SymPy `sp.limit` and Mathematica `Normal[Series[..., {z, 0, 0}]]` — distinct algorithmic routes, outlet consistency `Pi = M_s` at `S_q = 0`); 138 F1 → banner mass-fix only. **One follow-up squashed in**: notes 134 line 140 trailing-9 typo (`605429` → `605428` to match boxed forms at lines 86/92 and the scripts). **Codex bypassed entirely** (tenth consecutive zero-redirection batch); orchestrator-direct math-authority workflow held cleanly across all three clusters. **Pitfall #13 candidate added** (Mathematica parser fails on comment substring `g'(Pi_*)` adjacent to `*)` — workaround is ASCII labels like `piStar`). With IV.4 closed, **entire range 001-138 is now paper-aligned at v2 depth**. **No new paper card display titles deferred** in IV.4 — previously deferred IV.1/IV.2 paper titles remain as P3-13. See `redteam/batches/batch_IV4_v2.md`, `redteam/resolutions/batch_IV4_paper_alignment.md`. |
| 2026-05-27 | **Eleventh batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch IV.3 (stages 115-126, "Part IV.3 -- Core balance, DtN mixed, outlet, positive source") completed 2026-05-27. All 12 stages reached `verified` (2 clean: 120, 124 — status-only consolidation cards per MANIFEST `is_status_only_candidate`). **No checkpoints** in IV.3 range. **27 findings closed across 10 dirty stages.** **7 paper_misalignment items resolved per 4 user-gate clusters**, plus a 10-site orchestrator-direct banner-relabel sweep. **Cluster A (notes-side numerical typos)**: 4 notes files fixed in place (`notes/stages/moving_throat_pde_stage{121,122,126}_*.md` `168π² → 100π²` in 5 boxed surd occurrences; `notes/stages/moving_throat_pde_stage123_*.md` `228 → 160` denominator in 2 occurrences). Scripts confirmed correct in all four cases by 4 independent auditors who each re-derived the correct values; the notes' own quoted decimals match the scripts, not the typo'd surds. **Cluster B (118 λ sign flip — substantive script-side material_change)**: Internal script inconsistency where section IV's bilinear derivation gave `λ = −q* v_w0 𝓘_sq` (matching notes), but section V's closure asserted `λ_uniform = +(8√2/3)·…` (plus). Resolution: flipped script section V sign in both engines; F2 added 3 new asserts (K_q, g_s, λ from bilinear) with consistent signs. Downstream Schur reductions use `K_s K_q + λ²` and `(K_s g_q − λ g_s)²` — both sign-invariant under squaring — so no numerical propagation downstream; `upstream_stale` NOT flagged. **Cluster C (125 integral inequality)**: Paper card's `Output` is the integral inequality `0 ≤ 𝔤[σ] ≤ 1`; script only proved pointwise kernel bound. Resolution: parametric family `σ_a(z) = (a+1)(z/L)^a/L` added in both engines with `g_a = ∫σ_a·cos(πz/(2L)) dz` symbolic integration + endpoint asserts (a=0 → `g = 2/π`; a→∞ → 0). SymPy required a numerical proxy at a=100 because the hypergeometric closed form blocks `sp.limit`; Mathematica's `Limit` works directly. **Cluster D (117 consolidation card)**: 3 blocked items (transliteration + 2 tautological κ_c=1/3, γ_c=1/9) resolved via direction (a) cite-upstream-and-downgrade: added comment block citing stages 115 (Schur) and 116 (D/N eigenvalue, patched this batch); replaced misleading `expect_zero` tautological wrappers with `print("carrying forward (Stage 116/119): ...")` lines; F4 wired `classification_rows` booleans from sections 1-5 residuals (`nontrivial_compensated` anchored to `delta_core - delta_core_expected` series residual). **Other mechanical findings (17)**: 115 F1 mathematica_transliteration (independent parent-overlap reparametrization block via `frakR/frakG` with multiplicative-factor equivalence + Solve uses fresh `gVar`); 116 F3 hardcoded_result (`q'' + k²q = 0` BVP solved from first principles in both engines via `q(x) = sin(k x)` with k_W = π/(2L_W)); 116 F1 tautological_check (kappa_0 round-trip through geometric `4 L_W_required²/(π²a²)`); 116 F2 tautological_check (γ_0 extraction from D_bare's z⁵ coefficient — **orchestrator-corrected sign in directive**: `+I·coeff` not `-I·coeff`); 116 F4 insufficient_verification (SymPy tube-length assert added); 119 F1 tautological_check (rc → rhat² substitution link); 119 F2 insufficient_verification (T_m ± branch closed-form matches; Mathematica `stripCE` to remove ConditionalExpression heads); 121 F2 insufficient_verification (4 new closed-form asserts) + F3 script_missing_paper_claim (Ω_W identification); 122 F2 insufficient_verification (6 new asserts: compensation quadratic, defect closed form, natural off-compensation, traction ratios); 123 F2 banner only; 125 F2 mathematica_transliteration (.wl uses `Solve[balance==0, gSym]` instead of hand-written closed forms); 126 F2 paper_missing_script_claim (positivity asserts: SymPy `sp.calculus.util.minimum`; Mathematica boundary-value `expectZero` checks because `Minimize[]` returned unevaluated under monotone-decreasing assumption); 126 F3 insufficient_verification (interval check hardened to raise); 126 F4 stale_output banner. **Codex bypassed entirely** (ninth consecutive zero-redirection batch); orchestrator-direct math-authority workflow held cleanly across all four clusters. **Pitfall #12 candidate added** (Mathematica `Solve[expr == 0, frakG]` fails with "frakG is not a valid variable" when `frakG` is bound to its definition — always introduce a fresh symbol for Solve's target, then substitute back). Two directive-correction catches by orchestrator at edit time: (a) stage 116's `-I·coeff(z,5)` → `+I·coeff(z,5)` sign; (b) stage 115's multiplicative factor `(kS·kQ)/gS²` → `(kS·kQ + lam²)/(gS²·kQ)`. With IV.3 closed, **entire range 001-126 is now paper-aligned at v2 depth**. **No paper card display titles deferred** in IV.3 — the previously deferred IV.1/IV.2 paper titles remain as P3-13 / future paper-cleanup pass. See `redteam/batches/batch_IV3_v2.md`, `redteam/resolutions/batch_IV3_paper_alignment.md`. |
| 2026-05-27 | **Tenth batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch IV.2 (stages 103-114, "Part IV.2 -- Outgoing DtN, deformation, robustness, robin") completed 2026-05-27. All 12 stages reached `verified` (5 clean: 103 and 113 — status-only consolidation cards per MANIFEST `is_status_only_candidate`; 104, 110, 114 — clean substantively). Checkpoint stage 105 (chi_Q fix from outgoing DtN) passed at the higher-bar standard via genuinely independent Mathematica path (Apart-based partial-fraction round-trip + SeriesCoefficient operator form + Reduce-over-reals for chi_Q + polynomial inversion of Lambda·Y = -3 for deformed branch). **16 findings closed across 7 dirty stages (105/106/107/108/109/111/112). 5 paper_misalignment items resolved per 3 user-gate clusters, plus a 24-site orchestrator-direct banner-relabel sweep across all 10 IV.2 scripted stages.** Substantive: **Q1=(a) Cluster A (106 F1, 109 F3)** — script-side docstring carry-forward annotations naming neighbour stages that own the deferred Checks: 106's higher-odd item carry-forwards to stage 102, l=2 DtN fingerprint item carry-forwards to stage 104; 109's Robin/mixed-pole/even-coeff items reference downstream stages 110/111/112. No new asserts, no paper edits. Mirrors IV.1 Cluster A pattern. **Q2=(a) Cluster B (108 F1)** — substantive β-parameterized preservation submanifold added to both engines: builds `Lambda_gen(z) = S*Lambda_out(beta*z) + Sigma0 + Sigma2 z^2 + Sigma4 z^4 + i Sigma5 z^5`, re-solves (Σ₂,Σ₄) under canonical-even matching (now β-dependent), asserts the locus `Sigma_5 = S(1-β^5)/9 - Sigma_0/27`. Existing β=1 reduction (Class C) is subsumed as sanity. **One material_change: true at stage 108** (verification surface widened from β=1 reduction to general β; no derived value changed — Class C `Sigma_5 = -Sigma_0/27` still holds and is now a sanity-reduction check; downstream stages > 108 not marked `upstream_stale`). **Q3=(a) Cluster C (banner sweep)** — script-side 24-site banner sweep applied across 10 IV.2 scripted stages (104.py+104.wl, 105.py:3+28+105.wl, 106.py:3+25+106.wl, 107.py+107.wl, 108.py+108.wl, 109.py+109.wl, 110.py+110.wl, 111.py+111.wl, 112.py:3+54+112.wl, 114.py+114.wl) replacing stale "STAGE 87-097/Stage 87-97" labels from the prior renumber with current "STAGE 104-114". Paper card display titles `\section[Stage~108-119]{...}` deferred to PAPER_CLEANUP_TRACKER P3-13. **Other script-side findings (11)**: 105 F1 mathematica_transliteration rewrite (full .wl re-author with structurally distinct path); 106 F2 tautological_check (switched `K4 - 4 K2^2/K0` from script-construction tautology to `K0_target K4_target - 4 K2_target^2` testing the four target literals' mutual consistency); 106 F3 mathematica_transliteration (full retarded one-pole + omega^7 series-expansion rewrite); 106 F4 insufficient_verification (Δ_Q := chi_Q - 1 first-order sensitivity exercising the linearization slope -2G/(5c^5)); 107 F1 insufficient_verification (SymPy added Sigma2/Sigma4 exact-formula asserts for engine parity with Math twin); 108 F2 insufficient_verification (Math parse-bug fix: `chiArg /. beta -> 1 - 1` parsed as `beta -> 0`, never tested χ_arg(β=1) - replaced with `(chiArg /. beta -> 1) - 1`); 109 F1 tautological_check (added a5 closed-form `-5b/9 - a0/27` assert before substitution); 109 F2 mathematica_transliteration (linearization rewrite via numerator/denominator-separate series); 111 F1 mathematica_transliteration (independent chi_Q^mix re-derivation via geometric-series pole route + routes-agree cross-check); 112 F2 mathematica_transliteration (independent Stage-92 linearized cross-check extracting (a_0, a_5) from solB's deformation and solving the preservation condition for gamma_W = 1/9). **Codex bypassed entirely** (eighth consecutive zero-redirection batch); orchestrator-direct math-authority workflow held cleanly across all three clusters. **Pitfall #11 PASS-line discipline re-confirmed** — the 108 F2 Mathematica parse bug (`chiArg /. beta -> 1 - 1` is `chiArg /. (beta -> 0)`) had been passing the buggy assertion `0^5 = 0` for every prior run; only the auditor's structural read caught it. No new pitfall candidates from IV.2. With IV.2 closed, **entire range 001-114 is now paper-aligned at v2 depth**. **Deferred paper-side cleanup**: 12 paper card section titles in IV.2 (`paper/stages/stage_{103..114}.tex` line 1) use display Stage 108-119 from a prior renumber; see P3-13 above. See `redteam/batches/batch_IV2_v2.md`, `redteam/resolutions/batch_IV2_paper_alignment.md`. |
| 2026-05-27 | **Eighth batch processed under v2 paper-grounded auditor**: red-team v2 first-pass audit of batch III.5 (stages 085-090, "Part III.5 -- Quadrupole cancellation, loading ratio, verdict") completed 2026-05-27. All 6 stages reached `verified` (1 clean: 085). **2 substantive paper_misalignment items + 12-script orchestrator-direct banner-relabel sweep** (every III.5 stage carried stale "STAGE 68/069/70/070/71/071/72/072/73/073" banners from the pre-renumber numbering). Substantive: Q1=(a) Stage 087 status/checkpoint consolidation — paper's `\stagefield{Purpose}` reads "Stage 087 records that..." and `\stagefield{Inputs}` names stages 085-086; the cancellation chain that collapses dependence on `s_-, lambda_-, beta_0, mhat_-, Pi_tr, C_mix, Pe_req` to `rho_alpha` was performed and verified upstream in stages 081-086. Both engines' docstrings now name the upstream verifying files; F2 tautological "anchors" (comparing `rho_X - 1` to literals equal to `rho_X - 1`) replaced with `expect_close` cross-checks of rho_X against upstream stage-082 paper-quoted values; F3 transliteration marked won't-fix per F1=(a) — consolidation stages don't need a second independent derivation. Q3=(a) Stage 089 CHECKPOINT — paper boxes `Pe_req = 0` (eq. app-stage089-Pe-zero) as `\stagefield{Output}` but the link `Omega(Pe=0) = 1 ⇒ zeta_F1(0) = A_F1` was not script-verified (Omega is 0/0 at Pe=0, requires symbolic limit). Added `sp.limit(Omega, Pe, 0) - 1`, `sp.limit(zeta_F1, Pe, 0) - A_F1` (SymPy l'Hopital returns residual 2.17e-101, confirming the limit is genuinely computed and not a definitional substitution), explicit `Pe_req = sp.Integer(0)` + assertion at chain closure; both engines independently. F2 (Mathematica transliteration) rederives `peSuffChi`/`peFailChi` via `FindRoot[zetaF1[pe] == zetaTarget]` from notes-quoted `rho_target - 1` — independent path from SymPy which keeps the literals with provenance comment per pitfall #10 (sp.nsolve near `tan y` singularity instability documented during III.4). F3 (both engines) replaces tautological `rho_X - (1 + zeta_X)` with `expect_close`/`expectApprox` vs upstream stage-082 values. F4 (sympy hardcoded) resolved by provenance comment naming `scripts/output/moving_throat_pde_stage082_*_sympy_audit.txt` as source. **Two orchestrator hot-fixes on stage 088**: (a) SymPy — `Y_rho.subs(omega**2/Omega_Q**2, u)` failed silently because `sp.simplify` reshapes the denominator into `(Omega_Q**2 - omega**2)` form and the combined ratio is no longer a syntactic subexpression; fix substitute `omega**2 -> u * Omega_Q**2` then `sp.simplify` (more robust to canonicalization). (b) Mathematica — initial comment `(* ... stage085_*). Substitute rho_min. *)` was prematurely closed by the embedded `*)` substring at `stage085_*)`, causing `Syntax::sntx`; Mathematica's parser recovery reached `Exit[0]` with rc=0 but silently skipped the F1 assertion `Pi_tr_from_rho - (4/3) C_mix` and the regime trichotomy. **Verifier caught it from the missing PASS line count** (7 vs expected 9), confirming the verifier-prompt warning that "passing exec log is necessary but not sufficient." Fix: reword to "stage 085 files" — strip `*)` substring. **New pitfall #11 candidate**: Mathematica comments cannot contain `*)` substrings; verifiers must count expected PASS lines, not just check rc. **User-redirection rate 0** (sixth consecutive batch — Codex was bypassed entirely in III.5 per the III.4 availability lesson; orchestrator-direct math-authority workflow held up because audit + grep evidence was conclusive on both substantive questions). **Zero `material_change: true` flags** — all derived numerics unchanged (rho_alpha=4/3, zeta_req=1/3, Pi_tr=(4/3)C_mix, Pe_req=0); edits restructure paths to those values, add limit checks, add cross-engine independence. With III.5 closed, **entire range 001-090 is now paper-aligned at v2 depth**. See `redteam/batches/batch_III5_v2.md`, `redteam/resolutions/batch_III5_paper_alignment.md`. |
| 2026-05-26 | **Sixth batch processed under v2 paper-grounded auditor**: red-team v2 re-audit of batch III.3 (stages 061-072, "Part III.3 -- Microclosure, gain thresholds, equilibrium, walls") completed 2026-05-26. All 12 stages reached `verified` (5 clean: 061, 063, 065, 069 checkpoint, 071). **4 paper_misalignment items + 9 script-side findings** (4 `mathematica_transliteration`, 4 `insufficient_verification`, 1 `tautological_check`). The 4 paper_misalignment items: Q1=(a) Stage 062 scripts extended with second equality of boxed `G_micro = (ρ_* g_φ² N_φφ/(m c_{s,*}² K_X))·C_{σφ}²` plus Cauchy parameterization `Osp = cos(θ)·√(N_ss·N_pp)` proving `C_{σφ}² = cos²(θ) ∈ [0,1]` — orchestrator destination-verified that stage 063's audit directly consumes `C_sp_sq` at line 36 and uses Cauchy saturation at line 117; Q2=(a) Stage 062 σφ coupling sign flipped from `+Λ_φ σ φ` to `−Λ_φ σ φ` to match notes' `δV_conf = −g_φ χ_φ φ` derivation — orchestrator destination-verified that stage 064 already uses the minus convention at scripts/064 lines 63 and 146; Q3=(direct, no user gate) Stage 067 banner relabel `STAGE 50/050 → STAGE 067` (auditor written `Resolve` block was inconsistent with auditor's own report that explicitly stated "does not require user resolution"; orchestrator-applied directly with 3 string replacements); Q4=(direct via codex-invoke) Stage 072 banner relabel `STAGE 55/055 → STAGE 072` (5 string replacements). **User-redirection rate 0** this batch (fourth consecutive after II.1, III.1, III.2 also 0). **One orchestrator hot-fix on stage 064 Mathematica**: Codex iter1's F1 (general-H two-point check) prescribed `Integrate[hFun*chiSigmaFun^2, ...]` vs `gPhi^2*Integrate[chiPhi^2/hFun, ...]`. Algebraically the integrands are identical (`gPhi^2 chiPhi^2/hFun`), but Mathematica's `FullSimplify` cannot pull the constant `gPhi^2` outside `Integrate[...]` when the integrand contains unspecified symbolic functions. The exec-mathematica failed with surface form `-gPhi^2*Integrate[c^2/h] + Integrate[gPhi^2*c^2/h]`. Hot-fix: verify integrand equality first via `FullSimplify[hFun*chiSigmaFun^2 - gPhi^2*chiPhi^2/hFun] == 0` (which reduces cleanly), then define `thetaGeneral = lambdaGeneral = gPhi^2*i1Integral` directly. **New pitfall #9 candidate** (Mathematica `Integrate[]` does not factor constant multipliers when integrand has unspecified symbolic functions): when comparing integral identities, verify integrands FIRST, integral equality follows by pointwise equality. Promote to `codex.md` if recurs. **Pitfall #8 (heavy BVP `dsolve`) was preemptively promoted to `codex.md` before III.3 launched as defense**; no BVP verifications were prescribed in III.3 so #8 did not recur. Stage 068 (v1 `material_change: true` carried forward) returned **clean** under v2 audit — Solve-derived `Wfail_res`/`Wfail_match` preserved on resonance-corrected premises; Mathematica section 2 upgraded from `Solve` to `Reduce` without value change. **Zero `material_change: true` flags** this batch — Stage 062's added second-equality and Cauchy checks are additions; σφ sign flip changes only the displayed sign of `σ_*(φ)` (gain magnitude invariant); the 067/072 banner relabels are pure strings. With III.3 v2 closed, **entire range 001-072 is now paper-aligned at v2 depth**. Cosmetic follow-ups (non-blocking): stages 061-066, 068, 070, 071 carry pre-renumbering legacy banners ("STAGE 44" through "STAGE 53"); saved `.txt` outputs under `scripts/output/` and `mathematica/output/` for 066, 067, 068, 072 are stale relative to script mtimes (orchestrator's `redteam/exec_logs/` are the authoritative post-fix transcripts); stage 070 Mathematica `Print` block has a stale `(analytic 8/15 = 0.5333)` comment for `I_g` — actual integral is `14/15 ≈ 0.9333` (comment-only; both assertions use `IgNum` consistently so no false pass). See `redteam/batches/batch_III3_v2.md`, `redteam/resolutions/batch_III3_paper_alignment.md` (2 Codex recommendations with apply logs), `redteam/resolutions/codex_prompt_batch_III3.md`, `redteam/resolutions/codex_apply_batch_III3.md`. |
