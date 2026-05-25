# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the PDE ledger archive.

Use it together with `STAGE_PROVENANCE_INDEX.md`:

- `STAGE_PROVENANCE_INDEX.md` holds the exact raw note and audit artifact paths.
- `STAGE_VERIFICATION_COVERAGE.md` summarizes the current verification surface,
  exposes the main gaps, and gives us a stable baseline for audit planning.

Snapshot date: `2026-05-25` (batch III.4 close)

## Scope

- Canonical stage range: `001--253`
- Canonical stage source model:
  - Parts `I--VII`: compact stage cards
  - Part `VIII`: full derivation-stage files

## Coverage Totals

| Metric | Count |
|---|---:|
| Total stages in archive | 253 |
| SymPy audits present | 240 |
| Mathematica audits present | 165 |
| Numerical stress artifacts present | 15 |
| Stage-specific review notes present | 178 |
| Stages with no SymPy audit | 13 |
| Stages with no executable audit | 11 |
| Stages with Mathematica but no SymPy | 2 |

Mathematica counts above are presence counts only. They are **not**
independence counts. See `MATHEMATICA_MIRROR_POLICY.md` for the current rule
that distinguishes secondary replay coverage from genuinely independent
Mathematica mirrors.

The SymPy runner summary reports `241` passing files because it includes the
repo-level `moving_throat_pde_master_sympy_audit.py` in addition to the `240`
stage-level SymPy audits counted above.

### Red-team verification status

Distinct from "audit file present" (counted above) is "red-team verified" —
the stage has passed the `audit → directive → codex → verifier` pipeline run
out of `redteam/`, with both engines independently checking load-bearing
claims and a clean-context verifier agent confirming the directive's intent
was honored. See `redteam/BATCHES.md` for the live batch table.

As of `2026-05-25`: 84 of 253 stages red-team verified.

| Batch | Range | Stages | Verified | Date |
|---|---|---:|---:|---|
| I.1 | `001--012` | 12 | 12 | 2026-05-21 |
| I.2 | `013--023` | 11 | 11 | 2026-05-21 |
| II.1 | `024--036` | 13 | 13 | 2026-05-22 |
| III.1 | `037--048` | 12 | 12 | 2026-05-22 |
| III.2 | `049--060` | 12 | 12 | 2026-05-22 |
| III.3 | `061--072` | 12 | 12 | 2026-05-22 |
| III.4 | `073--084` | 12 | 12 | 2026-05-25 |
| III.5 onward | `085--253` | 169 | 0 | pending |

Cumulative findings closed: ~219 (`tautological_check` dominant overall,
`mathematica_transliteration` second; `hardcoded_result` rose sharply in
III.4 to 12 because the Family-1 numerology cluster 075-084 packs many
literal constants). Two stages still carry `material_change: true`: III.2
stage 060 (a `Csol = a/(exp(a*L)-1)` restructure of a previously-failing
`sp.solve`) and III.3 stage 068 (`Solve`-derived `Wfail_res`/`Wfail_match`
matching prior postulated forms symbolically); III.4 added no new
`material_change` flags (every derivation-route rewrite left printed
symbolic and numeric content byte-identical). Downstream consumers of 060
and 068 are still `pending` so no immediate cascade. See per-batch
summaries in `redteam/batches/batch_<ID>.md`.

### Linear projected-EM update

Snapshot addendum: `2026-05-11`

Stages `004--021` are now canonical linear stages for the projection-first
Maxwell integration, parent-action packet, and retained reduced one-port normal
form.  They are no longer counted as Stage `004` substages.  The old compact
Stage `004` reduced Maxwell/mixed calculation is retained as Stage `021`.

Stages `004--020` have file-for-file SymPy migrations from the derivation-only
`notes/em_projected` scripts through
`step_18_parent_throat_action_weak_axisym_packet`.  The `step_19_*`
branch-export packet and `step_20+` computational/runtime diagnostics remain
excluded from paper coverage.  File-for-file Mathematica mirrors for
Stages `004--020` have not been independently derived yet; Stage `021` retains
the legacy reduced Mathematica audit.

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | `001--023` | 23 | 23 | 6 | 2 | 6 |
| II | `024--036` | 13 | 13 | 13 | 0 | 13 |
| III | `037--090` | 54 | 53 | 54 | 2 | 54 |
| IV | `091--163` | 73 | 61 | 59 | 8 | 69 |
| V | `164--200` | 37 | 37 | 25 | 1 | 25 |
| VI | `201--218` | 18 | 18 | 2 | 0 | 2 |
| VII | `219--242` | 24 | 24 | 3 | 0 | 3 |
| VIII | `243--253` | 11 | 11 | 3 | 2 | 3 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica | 163 | See `STAGE_PROVENANCE_INDEX.md` for the file-by-file list. |
| SymPy only | 77 | Includes projected-EM Stages `004--020` and later stages whose Mathematica mirrors are not yet present. |
| Mathematica only | 2 | `084`, `093` |
| No executable audit | 11 | `103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153` |

## Constant Provenance Rule

Coverage counts are not trust grades.

Likewise, `Mathematica present` is not the same thing as `independent second
CAS derivation`. Repo-wide counts track coverage breadth; independence is now a
separate policy classification in `MATHEMATICA_MIRROR_POLICY.md`.

For this archive, an audit should be treated as insufficient until every
nontrivial constant used in it is classified as one of:

- `derived in audit`
- `carried forward with source anchor`
- `probe-only numeric value labeled`

Any unexplained literal should be treated as a verification defect, not a style
issue.

## Immediate Gaps

### 1. No executable audit

These stages currently have neither SymPy nor Mathematica nor numerical stress:

`103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`, `153`

These are the first places where the archive has mathematical content without an
executable backstop.

### 2. Mathematica without SymPy

These stages have a Mathematica artifact but no SymPy mirror:

`084`, `093`

These need reconciliation before we can make strong claims about dual-CAS
coverage.

### 3. SymPy-only frontier

The current SymPy-only region is:

`121--123`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

Operationally, the main contiguous Mathematica mirror still ends at Stage `187`,
with isolated later checkpoints now hardened at Stages `200`, `203`, `218`,
`221`, `239`, `242`, `243`, `248`, and `253`.

### 4. Review-coverage gap

Stage-specific review notes are missing for:

`121--124`, `188--199`, `201--202`, `204--217`, `219--220`, `222--238`, `240--241`, `244--247`, `249--252`

That means the late-stage archive is not only missing Mathematica coverage, but
also lacks the earlier review-note pattern that exists for most of Stages
`001--187`.

### 5. Numerical stress remains sparse

Numerical stress coverage exists only for the following stage families:

`003/021`, `045--048`, `058/170`, `106`, `125`, `142--144`,
`146--148`, `150--152`, `154`, `155--156`, `157`, `185--187`,
`248`, `253`

This is a narrow slice relative to the symbolic verification surface and should
be treated as targeted spot-checking, not broad numerical validation.

## How To Use This Baseline

1. Use `STAGE_PROVENANCE_INDEX.md` when you need the exact note or script path
   for a stage.
2. Use this baseline when deciding which verification gaps are structural and
   which are only missing a second engine or review note.
3. Use the part-level counts to prioritize the next audit wave.
4. Use the coverage classes to define the future-paper citation support set.
5. Use `CITATION_SUPPORT_SET.md` when deciding which gaps are most important for
   downstream citation hardening.
6. Use `CHECKPOINT_CONSTANT_PROVENANCE.md` for the growing no-magic-numbers log
   on the checkpoint subset.
7. Use `CHECKPOINT_TRUST_AUDIT.md` for the current stage-level trust baseline
   on that checkpoint subset.

## Recommended Next Verification Sequence

1. Use `CHECKPOINT_TRUST_AUDIT.md` as the current checkpoint trust baseline.
2. Reconcile the remaining repo-wide Mathematica-only outliers:
   `084`, `093`.
3. Use `MATHEMATICA_MIRROR_POLICY.md` when deciding whether a mirror gap is an
   execution-coverage gap or an independence gap.
4. Then widen the audit wave to the remaining repo-wide gaps.
5. Backfill executable audits for the remaining no-executable stages:
   `103`, `113`, `120`, `124`, `128`, `132`, `136`, `141`, `145`, `149`,
   `153`.
6. Use the new `248` / `253` harnesses as the template if numerical-stress
   coverage is widened beyond the current spot-check set.
