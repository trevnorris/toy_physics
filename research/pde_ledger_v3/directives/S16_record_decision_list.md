# S16 record: the moving-region step record and its register pass (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-24. **Author of the record:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S16_record_setup_lookups.md` ("L-n").

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). The record
author starts only after that. The record is physics-bearing prose. Two non-author legs (fresh Claude + Grok) review
it, source-first, until clear on what may be claimed. No build directive is in the record author's or the record
legs' packet (`CLAUDE.md`, step-record row). This follows the S14a part 1 record
(`directives/S14a_1_record_decision_list.md`).

## Deliverable
- `steps/S16_moving_region.md`: the S16 step record.
- `steps/_measurements/S16_record_measurements.md`, produced by a committed generator script. Every statement about a
  printed object, a source row or a register entry cites its retrieval there.
- An S16 pass in `SUBSTRATE_REQUIREMENTS.md`.
- The status and provenance lines of item 10.

None of the three files exists yet, and the register has no S16 pass (L2).

## Sources (read-only)
- **The step's decision list:** `directives/S16_decision_list.md`, revision 4, cleared at `0464b80a` and unchanged
  since (L1, L4). It holds the filing and framing, the user decisions, the supplied equations, the region, what the
  build computes, the model point, what goes to the record and its register pass, the inventories, what stays OPEN and
  the intended claim (L4). Its review dispositions are beside it in `directives/_measurements/`
  (`S16_decision_list_review_p{1,2,3,4}_disposition.md`), with its setup lookups `S16_setup_lookups.md`.
- **The accepted constructions** (accepted at `1d990ee6`, L1): `scripts/S16_moving_region.py`,
  `mathematica/S16_moving_region.wl`, their transcripts `scripts/out/S16_moving_region.out` and
  `mathematica/out/S16_moving_region.out`, and the run record `directives/_measurements/S16_build_run_record.md`.
  The transcripts are git-annex content: in a fresh checkout, run `datalad get` first (L3).
- **The build's acceptance:** `directives/_measurements/S16_build_review_r{0,1}_disposition.md`. Round 1 lists what it
  carries to the record (L5).
- **The accepted comparator** (accepted at `f08db34a`, L1): `scripts/S16_cross_engine_comparator.py`, its decision
  list `directives/S16_comparator_decision_list.md`, its production output
  `scripts/out/S16_cross_engine_comparator.out` (git-annex) with
  `directives/_measurements/S16_comparator_production_run.md`, its run record and builder report, and its review
  dispositions `directives/_measurements/S16_comparator_review_r{0,1,2}_disposition.md`. Round 2 lists what it
  carries to the record (L6).
- **The inventory sources** named in the step's decision list, under "Inventories" (L4).
- **The plan:** `V3_STEP_PLAN.md`, the S16 rows the step's decision list cites (L4, L12).
- **The O2 spec rows** the step's decision list cites: `directives/O2_SHARED_PHYSICS.md:133`, `:138` and `:607` (L8).
- **The gravity prior-art section** (accepted at `9d069620`, L1): `../../docs/medium_requirements_and_prior_art.md`,
  the section headed "ATTRIBUTION PASS — the gravity sector and S16", with its sources file
  `../../docs/_measurements/gravity_prior_art_sources.md` and its review dispositions
  `../../docs/_measurements/gravity_prior_art_review_r{0,1,2}_disposition.md` (L7).
- **The register:** `SUBSTRATE_REQUIREMENTS.md`: its entry schema, its requirement/postulate/defect distinction, and
  its population method, including both routes (L9).
- **The S14a part 2 record's acceptance:** commit `79511272` and
  `directives/_measurements/S14a_2_record_r2_review_disposition.md`, which lists what it carries to the next record
  list (L1, L10).

## User decisions this list records
- **Framing** (2026-10-19) and **"Keep S12"** (2026-10-20), as the step's decision list carries them ("Filing and
  framing").
- The **decisions of 2026-10-21**, as the step's decision list states them ("User decisions (2026-10-21)").
- **How the comparator review finished** (2026-10-24): a citation-only repair with a scoped check, then acceptance
  with each citation's claim stated as containment, as the comparator's round-1 and round-2 dispositions record them.

## What must be true
1. **The object and its standing.** The record states:
   - what the two accepted constructions print for each of the step's decision list items 1–5 ("What the build
     computes"), part by part, and the domain those items state;
   - what the build computes nothing of, as that section says;
   - the region and the model point, as the step's decision list states them, each labelled choice reported as the
     labelled choice it is;
   - which inputs are supplied, which are adopted premises, which are flagged identifications, and which stay live or
     OPEN.

   An adopted premise or a flagged identification is never presented as derived. The record claims no more than the
   step's decision list's intended claim. It states what of the plan's S16 rows this step prints and what stays with
   later steps, as its sources state them; where they do not say, the record says so.
2. **Cross-engine evidence, as printed.** Every statement about the two constructions' agreement or difference rests
   on the filed comparator output or the transcripts. Each statement cites its retrieval command, and the commands
   with their literal output go in the measurements file.
   - Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named stored
     object.
   - A question that needs computation beyond retrieval is listed as open, not computed.
   - The record covers every row class the comparator prints, every printed difference, and every tag or component
     that was not joined, not compared, not formed, or compared below its parsed count.
   - A difference is classified as representational only where a cited source settles it. Otherwise it is listed as
     an open cross-engine difference, with its owner or "no owner named". No difference is reconciled or explained
     away (`CLAUDE.md` M1).
3. **The instruments' limits travel with their output.** The record reproduces the comparator's printed statement of
   what it compares and does not compare. It also reproduces the limits and notes named in the build's and the
   comparator's acceptance dispositions, including what the build's round-1 disposition and the comparator's round-2
   disposition carry to this record (L5, L6). A claim reaches only as far as the compared objects reach.
4. **The inventories.** The record inventories what the step's decision list names under "Inventories", kept apart as
   the list requires, by pointing at each source row.
5. **Premises and their standing.** This covers the supplied equations, the region's choices, each user decision, and
   each premise the model point names. For each one, the record states:
   - its source and date;
   - its standing as the step's decision list gives it, quoted from that list;
   - which printed objects depend on it.

   The record gives no standing the list does not give. One that enters no computed object is listed as such.
6. **What remains open.** Each item the step's decision list retains as OPEN is stated with its owner, as the sources
   state them. None is closed here.
7. **The requirements this step files.** The step's decision list names what goes to the record and its register pass
   ("For the record and its register pass") and the debts for S22's table ("Filing and framing") (L4). The record and
   the pass file each one as that section states: once, against the owner the sources assign, saying so where none is
   assigned. They keep the printed objects apart from those debts, as "Filing and framing" does. Each entry:
   - cites what it rests on: a user decision, a source row or a printed object;
   - names no means of meeting it: no law, material, premise or candidate;
   - has a step id or "no owner named" as its target, and a status from the register's schema.

   If the register's rest-on test admits further entries whose requirement is a printed object, the pass adds them
   the same way. If it admits none, the record says why.
8. **Prior art.** The record cites the gravity prior-art section's rows as the section labels them, as that
   section's round-2 disposition carries (L7). A row is an oracle, never a premise (`CLAUDE.md` M3).
   - For each row the record cites, it quotes the row's conditions and S16's side by side.
   - Where a comparison between a printed object and a row's prior result needs only retrieval, the record makes it.
     Where it needs computation, or the conditions differ, the record says which, and lists the comparison as open.
   - S16 computes no passive factor. The record names the row the section names as the prior result for the open
     passive-factor requirement, with that row's conditions. It does not identify the two (the section's round-0
     disposition, R5, L7).
9. **The S16 register pass.** Populate the pass by the register's own schema, rest-on test and method, as written,
   both routes included. The second route covers the prior art the gravity section identifies (L7, L9).
   - Keep the distinction between a derived requirement and an adopted substrate premise.
   - Merge into existing entries wherever one covers the same requirement; don't duplicate.
   - Where this list and the register disagree, the register wins, and the disagreement is reported.
   - A reader of the register can tell that `R-S14a-01`'s id names its source step, not its deliverer, as the
     register's note does for `R-S1-01` (L9, L10).
10. **Status and provenance lines.**
    - The S14a part 2 record and register pass were accepted at `79511272` (L1). In the S14a part 2 record, in the
      register (its status and its S14a part 2 pass) and in `../../STATUS.md`, every statement that presents them as
      not yet accepted, or their review or commit as still to come, states that acceptance instead, heading or not.
      The statements found so far are in L10; that list is not claimed to be complete.
    - Nothing else in the S14a part 2 record or in the register's S14a part 2 pass changes.
    - `../../STATUS.md`'s open-questions table has a row for `R-S14a-01`, with its owners as the register's entry
      states them (L10, L11).
    - `../../STATUS.md` states where S16 stands, as the record does.
11. **The account of edits.** The record names every file this work changes outside the record, its measurements
    file and its generator, with each change. That account matches `git diff` against the commit the work starts
    from, and the measurements file retrieves that diff.

## Exclusions
- No new CAS, derivation, closure or profile choice. A comparison with prior art is made by retrieval only (item 8).
- No engine, comparator or transcript edit, and no rerun of any.
- Edit only:
  - the record, its measurements file and generator;
  - `SUBSTRATE_REQUIREMENTS.md`;
  - `steps/S14a_2_slow_limit.md`, for item 10 only;
  - `../../STATUS.md`.
- No commit and no push.

## STOP report (at most one page)
- The record's headline claims, each with its source.
- Each printed object of items 1–5, with its cross-engine status.
- Every printed cross-engine difference and unresolved coverage, and its classification.
- The requirements filed, each with its target.
- Register entry counts before and after, new entries by target, and existing entries that gained a source.
- The prior-art comparisons made, and those listed as open.
- Disagreements between sources, with both sides quoted.
- Open questions that need computation beyond retrieval.
- The account of edits (item 11).
