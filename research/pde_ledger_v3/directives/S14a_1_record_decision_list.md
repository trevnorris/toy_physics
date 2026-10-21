# S14a part 1 record: the projection step record and its register pass (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the record:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S14a_1_record_setup_lookups.md` ("L-n").

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). The record
author starts only after that. The record is physics-bearing prose. Two non-author legs (fresh Claude + Grok) review
it, source-first, until clear on what may be claimed. No build directive is in the record author's or the record
legs' packet (`CLAUDE.md`, step-record row). This follows the S12 record (`directives/S12_record_decision_list.md`).
Pass 1 was folded once (`directives/_measurements/S14a_1_record_decision_list_review_p1_disposition.md`; the
reviewed text is at `dcd6fd44`).

## Deliverable
- `steps/S14a_1_projection.md`: the S14a part 1 step record.
- `steps/_measurements/S14a_1_record_measurements.md`, produced by a committed generator script.
- An S14a part 1 pass in `SUBSTRATE_REQUIREMENTS.md`.
- The status and provenance lines of item 9.

None of the three files exists yet (L3).

## Sources (read-only)
- **The step's decision list:** `directives/S14a_1_decision_list.md`, revision 3, cleared at `e9eed614` (L1). It holds
  the filing and framing, the supplied equations, the user decision on what `J` counts, what the build computes, the
  model point, the inventories, what stays OPEN and the intended claim (L4). Its review dispositions are beside it in
  `directives/_measurements/` (`S14a_1_decision_list_review_p{1,2,3}_disposition.md`).
- **The accepted constructions** (accepted at `94fddf2f`, L1): `scripts/S14a_1_projection.py`,
  `mathematica/S14a_1_projection.wl`, their transcripts `scripts/out/S14a_1_projection.out` and
  `mathematica/out/S14a_1_projection.out`, and the run record `directives/_measurements/S14a_1_build_runs.md`. The
  transcripts are git-annex content: run `datalad get` first (L2).
- **The build's acceptance:** `directives/_measurements/S14a_1_build_review_r{0,1,2}_disposition.md`. Round 2 gives
  the authorship and order of the constructions and the known item N1 (L6).
- **The accepted comparator** (accepted at `de81b9fc`, L1): `scripts/S14a_1_cross_engine_comparator.py`, its decision
  list `directives/S14a_1_comparator_decision_list.md`, its production output
  `scripts/out/S14a_1_cross_engine_comparator.out` (git-annex) with
  `directives/_measurements/S14a_1_comparator_production_run.md`, its run record and builder report, and its review
  dispositions `directives/_measurements/S14a_1_comparator_build_review_r{0,1}_disposition.md`. Round 1's disposition
  lists what it carries to this record (L7).
- **The inventory sources** named in the step's decision list, under "Inventories" (L2).
- **The plan:** `V3_STEP_PLAN.md`, the rows on S13, S14a and S15 (L5).
- **The register:** `SUBSTRATE_REQUIREMENTS.md`: its entry schema, its requirement/postulate/defect distinction, and
  its population method, including both routes (L9).
- **The S12 record's acceptance:** commit `ac909aaf` and `directives/_measurements/S12_record_r3_review_disposition.md`
  (L1, L8).

## User decisions this list records
- **Framing** (2026-10-19) and **"Keep S12"** (2026-10-20), as the step's decision list carries them ("Filing and
  framing").
- **S13 is folded into S14a** (2026-10-21), as the step's decision list states it ("Filing and framing").
- **What `J` counts** (2026-10-21): both channels, each printed separately, as the step's decision list states it.

## What must be true
1. **The object and its standing.** The record states:
   - what the two accepted constructions print for the step's decision list items 1–4: the projected balances and
     their parts, and, in the specialization, each channel, `J`, the order lost, the outer-sphere flux, storage, the
     balance that relates them, the face channel minus the order lost, the outer-radius dependence and the
     conditioned velocity objects, and the dimensions;
   - the domain and the model point, as the step's decision list states them;
   - which inputs are supplied, which are adopted premises, which are flagged identifications, and which stay live or
     OPEN.

   An adopted premise or a flagged identification is never presented as derived. The spherical case is reported as
   the labelled specialization it is. The record claims no more than the step's decision list's intended claim. It
   states which of the plan's four S14a objects (`V3_STEP_PLAN.md:636–637`) this part prints and which stay with
   later parts, as its sources state them; where they do not say, the record says so. It states the standing of S14,
   S15 and S20 as the plan gives it (`:639–640`).
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
   what it compares and does not compare, and the limits and notes named in the build's and the comparator's
   acceptance dispositions, including what the comparator's round-1 disposition carries to this record (L6, L7). A
   claim reaches only as far as the compared objects reach.
4. **The inventories.** The record inventories what the step's decision list names under "Inventories", kept apart as
   the list requires, by pointing at each source row. Where the list asks the record to compare (v2's labels, the
   O2 spec's earlier projected material balance), the record does so by retrieval, or lists the comparison as open
   where it needs computation.
5. **Premises and their standing.** This covers U1 and its encoding, the material factor at a face, the
   specialization, the user decision on what `J` counts, and each premise the model point names as not imposed. For
   each one, the record states:
   - its source and date;
   - its standing as the step's decision list gives it, quoted from that list;
   - which printed objects depend on it.

   The record gives no standing the list does not give. One that enters no computed object is listed as such.
6. **What remains open.** Each item the step's decision list retains as OPEN is stated with its owner, as the sources
   state them. None is closed here.
7. **What a later law must meet.** If the register's rest-on test admits it, the pass adds entries whose requirement
   is a printed object. Each entry:
   - cites the printed objects it rests on;
   - names no means of meeting the requirement: no law, material, premise or candidate;
   - has a step id or "no owner named" as its target, and a status from the register's schema.

   If the rest-on test admits none, the record says why.
8. **The S14a part 1 register pass.** Populate the pass by the register's own schema, rest-on test and method, as
   written, both routes included.
   - Keep the distinction between a derived requirement and an adopted substrate premise.
   - Merge into existing entries wherever one covers the same requirement; don't duplicate.
   - Where this list and the register disagree, the register wins, and the disagreement is reported.
9. **Status and provenance lines.**
   - The S12 record and register pass were accepted at `ac909aaf` (L1). In the S12 record, the register (its status
     and its S12 pass) and `../../STATUS.md`, every line that describes S12's record or register pass as written for,
     under or awaiting review states that acceptance instead, heading or not (L8).
   - The S12 known items J1 and J2 (the S12 record's round-3 disposition, L8) are repaired, in the record and in the
     S12 pass. The account of edits names the two status-only edits J1 lists. The "18 → 18" sentence for the round-3
     repair cites the section that holds its count.
   - Nothing else in the S12 record or in the register's S12 pass changes.
   - The run record's construction-order sentence (N1) states the authorship and order the build's round-2
     disposition gives (L6).
   - `../../STATUS.md` states where S14a part 1 stands, as the record does.
10. **The account of edits.** The record names every file this work changes outside the record, its measurements
    file and its generator, with each change. That account matches `git diff` against the commit the work starts
    from, and the measurements file retrieves that diff.

## Exclusions
- No new CAS, derivation, oracle check, closure or profile choice.
- No engine, comparator or transcript edit, and no rerun of any.
- Edit only:
  - the record, its measurements file and generator;
  - `SUBSTRATE_REQUIREMENTS.md`;
  - `steps/S12_order_source.md`, for item 9 only;
  - `directives/_measurements/S14a_1_build_runs.md`, line 5 only;
  - `../../STATUS.md`.
- No commit and no push.

## STOP report (at most one page)
- The record's headline claims, each with its source.
- Each printed object of items 1–4, with its cross-engine status.
- Every printed cross-engine difference and unresolved coverage, and its classification.
- Register entry counts before and after, new entries by target, and existing entries that gained a source.
- Disagreements between sources, with both sides quoted.
- Open questions that need computation beyond retrieval.
- The account of edits (item 10).
