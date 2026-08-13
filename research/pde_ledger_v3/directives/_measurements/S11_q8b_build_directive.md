# Measurements behind `S11_q8b_build_directive.md`

Commands and literal output. CLAUDE.md rule 2. Regenerate; do not transcribe.
⛔ NO S11 RESULT VALUE APPEARS HERE — rule 5. Taken at HEAD `a8c8ddef` (engine `9fb45365`), 2026-08-13.

## 1 · The spec section

```
$ grep -n "Q8b" research/pde_ledger_v3/directives/S11_SHARED_PHYSICS.md
653:**Q8b — ⭐⭐ STRATA, defined.** All three `ROOT<r>_RANK_DROP_K_*`,
853:                                                  ⚠ point-local EVIDENCE under §Q8b's rule, ⛔ never the
```
§Q8b runs `:653-778` (through the corollary-4 empty case and orchestrator alignment note). Read in full
by the orchestrator before the directive was written.

## 2 · The engine's Q8b surface is a placeholder

```
$ grep -n "stratum_candidates\|STRATUM" research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
237:        parts.append(f"STRATUM{stratum}")
1432:    stratum_candidates = []
1442:            stratum_candidates.append(("RANK_DROP", locus))
1445:            stratum_candidates.append(("STACKED_DROP", locus))
1446:    emit_quantity(package, n, "STRATUM_ORDERING", sp.Tuple())

$ for t in _SOURCES _SOURCE_LOCUS_TAGS _DEFINING_EQUATIONS _FREE_PARAMETERS _POINT POINT_EVIDENCE \
      _CONSTANCY_CERTIFICATE _CHANGE_LOCUS COMPONENT_Q3_Q4_COVERAGE; do grep -c "$t" <engine>; done
0 0 0 0 0 0 0 0 0
```
Every §Q8b downstream object greps to **zero**. `STRATUM_ORDERING` is a literal empty `sp.Tuple()`
emitted unconditionally (`:1446`). `stratum_candidates` is written and never read.

## 3 · The pool collects two of the spec's three families

`emit_q8` (`:1429-1446`) appends `("RANK_DROP", locus)` and `("STACKED_DROP", locus)` for the K/COEFF/
JOINT variants per root. ⛔ The `ROOT_COINCIDENCE` family (§Q8b`:656` — both `_K` and `_COEFF`) is
emitted via `emit_locus` in the root section but its return value is discarded there — never pooled.

## 4 · The plumbing that already exists (supply, not to be rebuilt)

- `emit_locus` returns the live locus object (`:899-905`):
  `{"equations", "solution", "real_admissible", "base_quantity", "variables"}` — post-round-3 these
  carry attempted solutions and decided admissibility entries.
- The `stratum=` scope parameter threads through `emit_quantity`/`tag_name` (`:237`:
  `parts.append(f"STRATUM{stratum}")`).
- `NOT_DEFINED_ON_COMPONENT` is declared at `:34` and used nowhere.
- Round 3 (`9fb45365`) made the feeders real: admissibility tokens decide (ADMISSIBLE/EXCLUDED reachable)
  and every COEFF solve is attempted. Verified by two script legs + orchestrator probes
  (`_measurements/S11_engine_fix_round3_brief.md`).

## 5 · Runtime ruling (user, 2026-08-13)

> "we're computing something difficult. If it takes more than 600s that's fine. 600s was a guide, but if
> output is streaming so we can see progress then we can go over."

Matches spec `:1042-1048` (observable progress, not elapsed time; the engine emits and flushes per tag).
Baseline cell timings at `9fb45365`: MAIN D2 5.6 s · MAIN D5 ~371 s · XKIN_ANISO D2 ~171-188 s ·
XKIN_ANISO D3 >600 s streaming.

---

# Fold-time verifications (after both directive legs reported)

Leg reports: `~/.s11_build/q8b_directive_leg_grok.log` (probes `/tmp/s11q8b_leg_QN32/`),
`~/.s11_build/q8b_directive_leg_codex.log` (probes `/tmp/s11q8b_leg_8snV/`). Every folded finding was
re-verified by the orchestrator (rule 13).

## 6 · The dormant defect: an ADMISSIBLE entry does not promote (both legs; re-run)

```
$ timeout 300 python3 /tmp/s11q8b_leg_8snV/11_admitted_candidate_control.py
PY_S11_MAIN_D2_ROOT1_RANK_DROP_JOINT_REAL_ADMISSIBLE ['ADMISSIBLE']
PY_S11_MAIN_D2_STRATUM_ORDERING: Tuple()
ADMISSIBLE_ENTRY_COUNT 2
DOWNSTREAM_STRATUM_TAG_COUNT 0
```
And on the live acceptance cells the three families carry ZERO `ADMISSIBLE` entries (both legs,
independently; consistent with the orchestrator's round-3 postfix token counts) — so the first draft's
"per component promoted (if any)" acceptance was vacuous. ⇒ acceptance item 3 (promotion-path control)
now fails at baseline by construction.

## 7 · The bare-integer ban was OVERBROAD (Codex finding 2; verified against the spec read)

Spec `:721` scopes the count record to "**component-scoped** Q3/Q4 tag whose payload is an integer
count"; `:759-762` make `POINT_EVIDENCE_` tags the ordinary point-local recomputation. The first draft
banned bare integers in all stratum scope, which would have flagged a spec-CORRECT point-evidence
payload. Items 4/5 and acceptance 4 now scope the ban to component-wide answers.

## 8 · run_cell never publishes; stratum tags auto-enter the export (Codex finding 4; verified)

```
$ grep -n "if export and" research/pde_ledger_v3/scripts/S11_stray_longitudinal_sympy_audit.py
290:    if export and not local and package == PRIMARY_PACKAGE:
$ sed -n '/^def run_cell/,/^def /p' <engine> | grep -c "write_exports\|merged_export"
0
```
⇒ acceptance 5 (publish demo with the Q8b rows aboard) added; round 3's publish demos did not carry the
new volume.

## 9 · Rule-5 leak in the first draft (Codex finding 5)

"It will be slower than ~371 s; that is expected" was a forecast a builder could tune toward — and
already false (Codex measured the unchanged cell at 347.4 s). The line now demands literal timing with
no target. ⚠ The ~371 s figure remains above as a baseline observation ONLY; it is not in the directive.

## 10 · Notes, non-adoptions

- Grok finding 3 (change-locus non-promotion, spec `:748-750`): folded into item 5 verbatim-by-pointer,
  though the first draft's `:721-757` range already contained it — an explicit clause beats a mid-range
  implication.
- Codex's watchdog note: the directive tells the BUILDER to arm/kill the watchdog; leg prompts tell LEGS
  not to manage it. Different actors, both correct; no change.
- Codex counted MAIN D2 dispositions {EXCLUDED 12, UNDECIDED 8, NO_EXPOSED_BRANCH 1} vs Grok
  {EXCLUDED 12, UNDECIDED 10}: convention difference (opaque/no-branch entries), same substance —
  ADMISSIBLE 0 in both.
