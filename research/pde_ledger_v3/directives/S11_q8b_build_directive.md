# S11 §Q8b build — the strata. ⭐ A BUILD, not a fix: the section is unimplemented at every object.

## Authority and boundary

Edit `research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py` **only**. ⛔ Change no other
file. `CLAUDE.md` binds. `research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md` is the sole physics
authority; **§Q8b is `:653-778`** and it, not this directive, defines every object — where this
directive and the spec could be read differently, the spec wins and the divergence is a §10 report.

⭐ **Baseline is `9fb45365`** (fix round 3). The feeders are live: admissibility tokens decide, every
locus solve is attempted, and `emit_locus` returns the live locus object
(`{"equations", "solution", "real_admissible", "base_quantity", "variables"}`, `:899-905`). The
`stratum=` scope already threads through tag naming (`:237`), and a non-local `MAIN` tag in stratum
scope is automatically an export candidate (`:290`). `NOT_DEFINED_ON_COMPONENT` is declared (`:34`) and
never used. ⚠ Measured, and it is the dormant defect this build removes: **a candidate entry driven to
`ADMISSIBLE` through the live `emit_q8` still produces an empty `STRATUM_ORDERING` and zero downstream
stratum tags.** Measurements: `directives/_measurements/S11_q8b_build_directive.md`.

## What must become TRUE

1. **The candidate pool is the one §Q8b`:653-657` names — all three families.** Measured at baseline:
   `emit_q8` (`:1429-1446`) pools `RANK_DROP` and `STACKED_DROP` (K/COEFF/JOINT per root) and discards
   the pooled list; the `ROOT_COINCIDENCE` family (both `_K` and `_COEFF`, including the pair-named
   `ROOT_COINCIDENCE_R<p>_R<q>_*` loci) is emitted but never pooled.
2. **A stratum is what `:657-660` says it is** — a candidate component whose `_REAL_ADMISSIBLE` entry is
   `ADMISSIBLE` under §3 and for which the engine emits the exact point. ⭐ A candidate with `UNDECIDED`
   admissibility **remains an explicit locus-protocol coverage finding** — ⛔ neither silently discarded
   nor promoted. Emit the admitted components **as your CAS returns them, without canonicalising**.
3. **The per-stratum identity objects of `:664-678` exist**: `STRATUM_ORDERING`, and per component
   `_SOURCES` (tokens from exactly `RANK_DROP · STACKED_DROP · ROOT_COINCIDENCE`, that order, each at
   most once), `_SOURCE_LOCUS_TAGS`, `_DEFINING_EQUATIONS`, `_FREE_PARAMETERS`, `_POINT` (an assignment
   to **exactly** the solve variables the source loci name, satisfying every §3 assumption). ⭐ The point
   comes from the defining equations, so per §5 corollary 3 (`:676-678`) ⛔ **no residual between them is
   emitted**.
4. **The complete Q3 and Q4 objects are computed ON each component** (`:680-684`): the engine's own `M`
   and derived objects restricted to `_DEFINING_EQUATIONS` with `_FREE_PARAMETERS` left free, under the
   `STRATUM<s>` / `STRATUM<s>_ROOT<r>` scopes. ⛔⛔ **None of these payloads may be obtained by
   substituting `_POINT`** — that is the typed-object defect wearing a point evaluation, and it has an
   observable signature: **a component-wide payload, status, certificate or change locus may not depend
   on which point was chosen.** ⭐ The engine emits its own chart (`:686-719`): what it actually
   retained, unchanged; ⛔ no elimination is pinned, and the cross-engine comparison is the
   orchestrator's operation, ⛔ not this engine's obligation.
5. **Every component-wide integer count carries the `:721-757` record** — the three companion families
   (`_STATUS` exactly `CONSTANT · VARIES · UNDECIDED`, `_CONSTANCY_CERTIFICATE`, `_CHANGE_LOCUS_*`) and
   the single payload shape `{STATUS_TOKEN, VALUE}` with `NOT_DEFINED_ON_COMPONENT` where no
   component-wide count was obtained. ⛔ Never a bare integer **as a component-wide answer**; ⛔ never a
   payload shape that varies by status; ⛔ no component count copied or inferred from a point evaluation.
   ⚠ **Scope, per `:721` and `:759-762`: this record applies to component-scoped counts.
   `POINT_EVIDENCE_` counts retain the ordinary point-local Q3/Q4 payload form** — no `_STATUS`,
   certificate or change-locus families exist for them. ⭐ Every named tag emits **unconditionally**
   (corollary 4). `STRATUM<s>_COMPONENT_Q3_Q4_COVERAGE` carries exactly one of its two tokens
   (`:752-757`). ⛔ Per `:748-750`, a `_CHANGE_LOCUS_*` sub-locus is reported **only as an object** —
   never promoted to a stratum, never recursed on, never canonicalised.
6. **The point-local evidence set exists under the `:759-771` naming rule** (`POINT_EVIDENCE_` inserted
   before each unscoped name, after `ROOT<r>_` for root-scoped) — evidence at the point, ⛔ never the
   component's answer, and emitted even beside `INCOMPLETE_COVERAGE`.
7. **A package with no stratum emits `STRATUM_ORDERING` as the empty list** (`:773-775`) — the tag is
   never omitted, and the candidates' typed branch and real-status objects are what distinguish proved
   exclusion from undecided coverage from never-computed.

## Boundaries

- ⛔ No memory cap, no timeout, no handler that swallows a failure to make a run finish.
- ⛔ Do not change `PACKAGE_DIMS`, the `D` range, or any package.
- ⭐ **The §§4–8 objects and the round-3 behaviors are the feeders — they must not change**: attempted
  solves, decided tokens, equation-derived degenerate tests, the publish placement. A Q8b need that
  seems to require changing one is a §10 report, ⛔ not an edit.
- ⛔ **No expected value and no acceptance criterion referencing one** (rule 5): not a stratum count, not
  an admissibility outcome, not a rank, not a timing. ⛔ Do not treat any prior output as a reference — a
  changed value is a §10 finding, never something to tune away.
- ⛔ Do not add a registry, completeness map or exit policy. The run's own observed bookkeeping is not a
  registry.
- ⭐⭐ **Runtime (user, 2026-08-13): a cell may exceed 600 s provided its output is streaming.** The
  engine emits and flushes per tag; a long stretch with **no output** is the failure
  (spec `:1042-1048`), elapsed time is not. ⛔ Do not narrow scope, drop a component, or weaken an object
  to fit a time number (CLAUDE.md rule 11).

## Acceptance

⚠ Checked against the items. ⛔ The acceptance cells may legitimately admit no live component — that is
why item 3 below drives the path deliberately: **at baseline it fails by measurement** (an `ADMISSIBLE`
entry does not promote), so no placeholder or disposition-only half-build can green this section.
1. `/home/trevnorris/.s11_build/repro_d5.py` runs `run_cell('MAIN', 5)` to completion; literal stdout
   and timing. ⭐ The run must stream; report what it measured, ⛔ no timing target exists.
2. **The pool and its dispositions, `MAIN` D2 and `XKIN_ANISO` D2**: report every candidate disposition
   **in 1–1 correspondence with the emitted `_REAL_ADMISSIBLE` entries of exactly the three `:653-657`
   families** — including the pair-named `ROOT_COINCIDENCE_R<p>_R<q>_*` tags — citing each literal tag:
   promoted (with its complete per-stratum tag set), excluded, or an undecided coverage finding.
   ⛔ An accounting against the builder's own candidate list, or an empty `STRATUM_ORDERING` beside an
   unaccounted entry, fails.
3. **The promotion-path control, on a `/tmp` copy — ⛔ never in-tree**: for **each** of the three source
   families, drive at least one candidate entry to `ADMISSIBLE` through the engine's **real** pool and
   promotion path (never by hand-writing the output tags) and demonstrate with literal stdout: the
   candidate enters `STRATUM_ORDERING`; the complete identity objects emit; the component-scoped Q3/Q4
   payloads carry the retained free parameters; one count shows its full item-5 record; the coverage
   token emits; the point-evidence set emits at the point. ⭐ **Then change only the chosen point and
   re-run: every component-wide payload, status, certificate and change locus must be unchanged, while
   the point-evidence set recomputes at the new point.** ⛔ A build whose component payloads move with
   the point fails; a build that cannot run this control has not wired the path.
4. **Payload-form criteria, wherever stratum tags appear in any acceptance stream**: a bare integer as a
   **component-wide answer** fails; a `POINT_EVIDENCE_` count in the ordinary point-local form is
   **correct** and must not be flagged; a payload shape that varies by status fails; a missing companion
   family for a component-wide count fails. Report the emitting line numbers.
5. **The publish carries the new volume**: on a reduced-`PACKAGE_DIMS` `/tmp` copy whose `MAIN` rows
   include the new stratum tags, demonstrate the round-3 publish path end-to-end — the mid-loop publish
   completes with the Q8b rows aboard, the publish-time record and post-sweep `RUN_PAIRS`/
   `SKIPPED_PAIRS` semantics hold, and a publish failure still attributes to `PUBLISH`. Literal stdout;
   long and streaming is fine.
6. ⛔ Do not run the full package loop on the in-tree engine (OOM hazard, now heavier). ⭐ Any engine-cell
   run happens with `/tmp/s11_watchdog.sh` armed by the builder (capture its pid; kill that pid when the
   run ends — ⛔ never `pkill` by name).

## Deliverable

The edited script, and a note saying per item what changed, plus every value that moved anywhere and
where you reported it.
