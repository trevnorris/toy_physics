# Future-Paper Citation Support Set

This document defines the load-bearing subset of the PDE ledger archive that
future papers are most likely to cite.

It is not a replacement for `STAGE_PROVENANCE_INDEX.md` or
`STAGE_VERIFICATION_COVERAGE.md`.

- `STAGE_PROVENANCE_INDEX.md` keeps the exact note/script paths.
- `STAGE_VERIFICATION_COVERAGE.md` describes repo-wide coverage.
- `CITATION_SUPPORT_SET.md` says which anchors and checkpoint stages deserve the
  strongest audit attention because they are the most likely downstream support
  surface.

Snapshot date: `2026-04-20`

## Selection Rule

The support set is chosen using three rules:

1. Prefer the public result anchors already named in the paper's downstream
   citation crosswalk.
2. Within each anchor family, choose checkpoint stages that either:
   - define the reusable packet,
   - give the terminal local verdict,
   - or establish the status boundary that downstream papers must not overclaim.
3. Treat checkpoint stages as audit control points, not preferred public
   citations. Downstream papers should still cite the stable `MTDC-T*` anchors
   first and use stage numbers only for detailed audit routing.

## Current Support-Set Snapshot

The current checkpoint subset contains `25` stages:

- `25` have SymPy coverage,
- `25` have Mathematica coverage,
- `4` have numerical stress coverage,
- `25` have stage-specific review notes,
- `0` have no executable audit,
- `0` are Mathematica-only,
- `0` are SymPy-only and also missing stage review.

This makes the support set much narrower than the full `001--236` archive, but
still large enough to represent the derivation surface future papers are most
likely to reuse.

## Tier A: Recurrent Public Citation Surface

These are the anchors most likely to recur across framework, PN, retarded, and
quotient papers.

| Public handle(s) | Checkpoint stages | Why this is load-bearing | Current control gap |
|---|---|---|---|
| `MTDC-T1`, `MTDC-T2` | `001`, `002` | Parent ontology, geometry lift, and breathing reduction are the entry conditions for nearly every downstream derivation. | Canonical stage prose now matches the dual-CAS verification baseline; numerical stress is still absent. |
| `MTDC-T3` | `003` | First microscopic BdG support kernel and first pole-shift packet. | Dual-CAS and shared numerical stress now run cleanly; no open review gap. |
| `MTDC-T4` | `005` | Outgoing-normalization bridge and invariant grouped-response product. | Exact symbolic bridge and Stage-4 dictionary round-trip now close cleanly in both CAS layers; numerical layer still absent. |
| `MTDC-T5` | `006`, `007`, `019` | Full grouped bundle, isotropy/splitting law, and support-feasibility frontier. | All three checkpoints now have clean dual-CAS symbolic support; numerical stress is still absent. |
| `MTDC-T8` with emphasis on `T8.1`, `T8.2`, `T8.3`, `T8.5` | `079`, `088`, `095`, `146` | Geometry-lane verdict, exact fixing of `chi_Q`, compensation law, and first explicit off-family drift. | All checkpoints now have clean dual-CAS symbolic support; numerical stress is still absent. |
| `MTDC-T9` with emphasis on `T9.3`, `T9.6` | `168`, `183` | Microscopic monomial coordinates and the four-scalar final verdict packet `Delta_full`. | `168` now has strong symbolic + numerical coverage; `183` still lacks numerical stress. |

## Tier B: Branch, Search, and Placement Surface

These are the anchors future realization and same-charge papers are most likely
to cite directly.

| Public handle(s) | Checkpoint stages | Why this is load-bearing | Current control gap |
|---|---|---|---|
| `MTDC-T10` with emphasis on `T10.1`, `T10.4` | `186`, `201` | Scalar graph-slice theorem and final local mixed-ray closure theorem define the declared realization-search class. | Both checkpoints now have symbolic mirrors and review; numerical stress is still absent. |
| `MTDC-T11` with emphasis on `T11.1`, `T11.6`, `T11.7` | `204`, `222`, `225` | Resonance/linewidth gate, physical `(U,V)` chart, and final twin-support placement/orbit-lock compiler. | All three checkpoints now have symbolic mirrors and review; numerical stress is still absent. |
| `MTDC-T12` with emphasis on `T12.1`, `T12.3`, `T12.5` | `226`, `231`, `236` | Relaxed recovery map, dynamic event-chain compiler, and material-threshold/export packet. | `231` and `236` now also have dedicated dual-CAS numerical stress; `226` remains symbolic-only, which is acceptable unless future papers start using it quantitatively. |

## Tier C: Dependency Anchors That Still Need Explicit Control

These are less likely to be the public headline citation, but they remain
load-bearing dependencies inside the archive logic.

| Public handle(s) | Checkpoint stages | Why this is load-bearing | Current control gap |
|---|---|---|---|
| `MTDC-T6` with emphasis on `T6.4` | `034` | Lowest-twin sufficiency criterion is the early support-feasibility checkpoint that later branch arguments build on. | Symbolic mirrors and review exist; numerical layer absent. |
| `MTDC-T7` with emphasis on `T7.2`, `T7.7` | `052`, `072`, `073` | Reduced support/source verdict, explicit minimal-isotropic Family-1 verdict, and theorem-status boundary. | All checkpoints now have clean dual-CAS symbolic support; no numerical stress at the checkpoints. |

## What This Means For Audit Order

The first audit wave should not spread evenly across all stages.

It should proceed in this order:

1. Use `CHECKPOINT_TRUST_AUDIT.md` as the current checkpoint trust baseline.
2. Decide which checkpoint stages need numerical-stress hardening beyond the
   four cases that now have stress artifacts: `003`, `168`, `231`, and `236`.

## Operational Rule

Until a checkpoint stage has at least:

- a trustworthy symbolic audit,
- explicit constant provenance,
- and a human review classification,

it should be treated as useful for internal orientation but not yet as
fully-hardened citation support.
