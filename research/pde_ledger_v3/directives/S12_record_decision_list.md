# S12 record: the order-source step record and its register pass (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-20. **Author of the record:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`.

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). The record
author starts only after that. The record is physics-bearing prose. Two non-author legs (fresh Claude + Grok) review
it, source-first, until clear on what may be claimed. No build directive is in the record author's or the record
legs' packet (`CLAUDE.md`, step-record row). This follows the S9b record (`directives/S9b_record_directive.md`).
Pass 1 was folded once (`directives/_measurements/S12_record_decision_list_review_p1_disposition.md`).

## Deliverable
- `steps/S12_order_source.md`: the S12 step record.
- `steps/_measurements/S12_record_measurements.md`, produced by a committed generator script.
- An S12 pass in `SUBSTRATE_REQUIREMENTS.md`.
- The status lines that describe the accepted light-sector closeout as awaiting review, updated to its acceptance
  (item 9).

## Sources (read-only)
- **The step's decision list:** `directives/S12_decision_list.md`, revision 4, cleared at `b08673a2`. It holds the
  supplied equations, the user premises U1–U4, how U1 enters, what the build computes, the inventories, the model
  point, what stays OPEN, and the intended claim. Its review dispositions are beside it in
  `directives/_measurements/` (`S12_decision_list_review_p{1..5}_disposition.md`).
- **The accepted constructions** (accepted at `a3da43d9`): `scripts/S12_order_source.py`,
  `mathematica/S12_order_source.wl`, their transcripts `scripts/out/S12_order_source.out` and
  `mathematica/out/S12_order_source.out`, and the run record `directives/_measurements/S12_build_run_record.md`.
  The transcripts are git-annex content: run `datalad get` first.
- **The build's acceptance:** `directives/_measurements/S12_build_review_r{0,1,2}_disposition.md`. Round 2's notes
  name limits for the record and the comparison.
- **The accepted comparator** (accepted at `fe485b73`): `scripts/S12_cross_engine_comparator.py`, its decision list
  `directives/S12_comparator_decision_list.md`, its production output `scripts/out/S12_cross_engine_comparator.out`
  (git-annex) with `directives/_measurements/S12_comparator_production_run.md`, its run record and builder report,
  and its review dispositions `directives/_measurements/S12_comparator_build_review_r{0,1,2}_disposition.md`. Round
  2's disposition records the user's acceptance and the engine-row note.
- **The inventory sources** named in the step's decision list, under "Inventories".
- **The register:** `SUBSTRATE_REQUIREMENTS.md`: its entry schema, its requirement/postulate/defect distinction,
  and its population method, including both routes.
- **The light-sector closeout acceptance:** commit `e070c7f3` and
  `directives/_measurements/light_sector_closeout_review_r2_disposition.md`.

## User decisions this list records
- **Framing** (2026-10-19): S12 states what any conversion law must satisfy and leaves the forms open. It picks no
  drain law and does not define the brane. Each force states the flow it needs, S21 compares them, and a mismatch
  reopens S9b (the step's decision list, "Filing and framing").
- **Order carried along a face** (2026-10-20): "Leave it out for now". It is not in the object, and it stays an open
  possibility (the step's decision list, "Retained OPEN").
- **The comparator's engine row** (2026-10-20): "Accept, note the engine row" (comparator review round 2
  disposition).

## What must be true
1. **The object and its standing.** The record states:
   - what the two accepted constructions print for the step's decision list items 1–3, on each part of the domain:
     outside below, inside, outside above, and each face, including what the total-material balance requires at each
     face;
   - the domain and model point, as the step's decision list states them;
   - which inputs are supplied, which are adopted premises, which are flagged identifications, and which stay live or
     OPEN.

   An adopted premise or a flagged identification is never presented as derived. Item 2 is reported as the labelled
   conditional specialization it is. The record does not present the object as a drain law, a throat solve or a flow.
2. **Cross-engine evidence, as printed.** Every statement about the two constructions' agreement or difference rests
   on the filed comparator output or the transcripts. Each statement cites its retrieval command, and the commands
   with their literal output go in the measurements file.
   - Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named stored
     object.
   - A question that needs computation beyond retrieval is listed as open, not computed.
   - The record covers every row class the comparator prints, every printed difference, and every row that was not
     joined, not compared, not formed, or compared below its parsed count.
   - A difference is classified as representational only where a cited source settles it. Otherwise it is listed as
     an open cross-engine difference, with its owner or "no owner named". No difference is reconciled or explained
     away (`CLAUDE.md` M1).
3. **The instruments' limits travel with their output.** The record reproduces the comparator's printed statement of
   what it compares and does not compare, and the limits and notes named in the build's and the comparator's
   acceptance dispositions, the engine-row note included. A claim reaches only as far as the compared objects reach.
4. **The inventories.** The record inventories what the step's decision list names under "Inventories", kept apart as
   the plan requires, by pointing at each source row. It answers, or lists as open with an owner, the question the
   list raises there about v2's return controllers.
5. **Premises and their standing.** This covers each of U1–U4, how U1 enters (the support condition and U1's
   re-ordering clause, which it does not encode), the "Fields at a face" rule, and the fast-ordering assumption of
   item 2. For each one, the record states:
   - its source and date;
   - its standing as the step's decision list gives it, quoted from that list;
   - which printed objects depend on it.

   The record gives no standing the list does not give. One that enters no computed object is listed as such.
6. **What remains open.** Each item the step's decision list retains as OPEN is stated with its owner, as the sources
   state them. None is closed here.
7. **What a conversion law must meet.** If the register's rest-on test admits it, the S12 pass adds entries whose
   requirement is a printed object: what any conversion law, face law or order-flux form must satisfy on each part
   of the domain. Each entry:
   - cites the printed objects it rests on;
   - names no means of meeting the requirement: no law, material, premise or candidate;
   - has a step id or "no owner named" as its target, and a status from the register's schema.

   If the rest-on test does not admit one, the record says why.
8. **The S12 register pass.** Populate the S12 pass by the register's own schema, rest-on test and method, as
   written, both routes included.
   - Keep the distinction between a derived requirement and an adopted substrate premise.
   - Merge into existing entries, `R-S12-01`, `R-S12-02` and `R-S21-01` included; don't duplicate them.
   - Carry the framing decision above to the entry that owns the flow comparison.
   - Where this list and the register disagree, the register wins, and the disagreement is reported.
9. **Status lines.** The light-sector closeout was accepted at `e070c7f3`. In the register, the S9b flow-only oracle
   comparison record and `../../STATUS.md`, every line that describes that work as awaiting review, in whatever
   words, states its acceptance instead. This includes a line that makes later work wait on that review.
   `../../STATUS.md` also states where S12 stands, as the record does.

## Exclusions
- No new CAS, derivation, oracle check, closure or profile choice.
- No engine, comparator or transcript edit, and no rerun of any.
- Edit only the record, its measurements file and generator, `SUBSTRATE_REQUIREMENTS.md`, the status line of
  `steps/S9b_flow_oracle_comparison.md`, and `../../STATUS.md`.
- No commit and no push.

## STOP report (at most one page)
- The record's headline claims, each with its source.
- Each part of the domain's printed object, with its cross-engine status.
- Every printed cross-engine difference and unresolved coverage, and its classification.
- Register entry counts before and after, new entries by target, and existing entries that gained a source.
- Disagreements between sources, with both sides quoted.
- Open questions that need computation beyond retrieval.
