# Stage Verification Coverage Baseline

This document is the stage-coverage control sheet for the PDE ledger archive.

Use it together with `STAGE_PROVENANCE_INDEX.md`:

- `STAGE_PROVENANCE_INDEX.md` holds the exact raw note and audit artifact paths.
- `STAGE_VERIFICATION_COVERAGE.md` summarizes the current verification surface,
  exposes the main gaps, and gives us a stable baseline for audit planning.

Snapshot date: `2026-04-21`

## Scope

- Canonical stage range: `001--236`
- Canonical stage source model:
  - Parts `I--VII`: compact stage cards
  - Part `VIII`: full derivation-stage files

## Coverage Totals

| Metric | Count |
|---|---:|
| Total stages in archive | 236 |
| SymPy audits present | 223 |
| Mathematica audits present | 165 |
| Numerical stress artifacts present | 15 |
| Stage-specific review notes present | 175 |
| Stages with no executable audit | 11 |
| Stages with Mathematica but no SymPy | 2 |

Mathematica counts above are presence counts only. They are **not**
independence counts. See `MATHEMATICA_MIRROR_POLICY.md` for the current rule
that distinguishes secondary replay coverage from genuinely independent
Mathematica mirrors.

## Coverage By Part

| Part | Stage Range | Total | SymPy | Mathematica | Numerical | Review |
|---|---|---:|---:|---:|---:|---:|
| I | `001--006` | 6 | 6 | 6 | 2 | 6 |
| II | `007--019` | 13 | 13 | 13 | 0 | 13 |
| III | `020--073` | 54 | 53 | 54 | 2 | 54 |
| IV | `074--146` | 73 | 61 | 59 | 8 | 69 |
| V | `147--183` | 37 | 37 | 25 | 1 | 25 |
| VI | `184--201` | 18 | 18 | 2 | 0 | 2 |
| VII | `202--225` | 24 | 24 | 3 | 0 | 3 |
| VIII | `226--236` | 11 | 11 | 3 | 2 | 3 |

## Coverage Classes

| Coverage class | Count | Stage ranges |
|---|---:|---|
| SymPy + Mathematica + numerical | 15 | `003--004`, `028`, `041`, `089`, `108`, `125`, `129`, `133`, `137--138`, `140`, `168`, `231`, `236` |
| SymPy + Mathematica only | 148 | `001--002`, `005--027`, `029--040`, `042--066`, `068--073`, `074--075`, `077--079`, `080--085`, `087--088`, `090--100`, `101--102`, `109--110`, `112--114`, `116--118`, `120--123`, `126--127`, `130--131`, `134--135`, `139`, `141--167`, `169--170`, `183`, `186`, `201`, `204`, `222`, `225--226` |
| SymPy only | 60 | `104--106`, `171--182`, `184--185`, `187--200`, `202--203`, `205--221`, `223--224`, `227--230`, `232--235` |
| Mathematica only | 2 | `067`, `076` |
| No executable audit | 11 | `086`, `096`, `103`, `107`, `111`, `115`, `119`, `124`, `128`, `132`, `136` |

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

`086`, `096`, `103`, `107`, `111`, `115`, `119`, `124`, `128`, `132`, `136`

These are the first places where the archive has mathematical content without an
executable backstop.

### 2. Mathematica without SymPy

These stages have a Mathematica artifact but no SymPy mirror:

`067`, `076`

These need reconciliation before we can make strong claims about dual-CAS
coverage.

### 3. SymPy-only frontier

The current SymPy-only region is:

`104--106`, `171--182`, `184--185`, `187--200`, `202--203`, `205--221`, `223--224`, `227--230`, `232--235`

Operationally, the main contiguous Mathematica mirror still ends at Stage `170`,
with isolated later checkpoints now hardened at Stages `183`, `186`, `201`,
`204`, `222`, `225`, `226`, `231`, and `236`.

### 4. Review-coverage gap

Stage-specific review notes are missing for:

`104--107`, `171--182`, `184--185`, `187--200`, `202--203`, `205--221`, `223--224`, `227--230`, `232--235`

That means the late-stage archive is not only missing Mathematica coverage, but
also lacks the earlier review-note pattern that exists for most of Stages
`001--170`.

### 5. Numerical stress remains sparse

Numerical stress coverage exists only for:

`003--004`, `028`, `041`, `089`, `108`, `125`, `129`, `133`, `137--138`,
`168`, `231`, `236`

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
   `067`, `076`.
3. Use `MATHEMATICA_MIRROR_POLICY.md` when deciding whether a mirror gap is an
   execution-coverage gap or an independence gap.
4. Then widen the audit wave to the remaining repo-wide gaps.
5. Backfill executable audits for the remaining no-executable stages:
   `086`, `096`, `103`, `107`, `111`, `115`, `119`, `124`, `128`, `132`,
   `136`.
6. Use the new `231` / `236` harnesses as the template if numerical-stress
   coverage is widened beyond the current spot-check set.
