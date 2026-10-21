# S14a part 2 record: the slow limit and far flux, recorded as a requirement (pre-builder decision list)

**Author:** Claude (orchestrator), 2026-10-21. **Author of the record:** Codex. `AGENTS.md` and `CLAUDE.md` apply.
Paths are relative to `research/pde_ledger_v3/`. The lookups behind this list are in
`directives/_measurements/S14a_2_record_setup_lookups.md` ("L-n") and in the review lookups
`directives/_measurements/S14a_2_record_decision_list_review_p1_lookups.md` ("N-n").

**Gate.** This is a pre-builder decision list. It gets one Codex + Grok pass and is folded once (G2). The record
author starts only after that. The record is physics-bearing prose. Two non-author legs (fresh Claude + Grok) review
it, source-first, until clear on what may be claimed. This follows the S14a part 1 record
(`directives/S14a_1_record_decision_list.md`). The G2 pass is complete and this revision is its one fold
(`directives/_measurements/S14a_2_record_decision_list_review_p1_disposition.md`; reviewed text at `ee6791c1`).

## Deliverable
- `steps/S14a_2_slow_limit.md`: the S14a part 2 step record.
- `steps/_measurements/S14a_2_record_measurements.md`, produced by a committed generator script.
- An S14a part 2 pass in `SUBSTRATE_REQUIREMENTS.md`.
- The status and provenance lines of item 6.

None of the three files exists yet (L1).

## Sources (read-only)
- **The step's decision list:** `directives/S14a_2_decision_list.md`. Revision 1 is at `d72bb551` and revision 2 at
  `945a1cd2` (L2). Its review dispositions and lookups are beside it in `directives/_measurements/`:
  `S14a_2_decision_list_review_p{1,2}_{disposition,lookups}.md`, and its setup lookups `S14a_2_setup_lookups.md`.
- **The review witnesses:** the review legs' scripts and their literal stdout, filed byte for byte in
  `directives/_measurements/S14a_2_review_witnesses/`, with their sha256 in
  `directives/_measurements/S14a_2_review_witnesses.sha256` (L3). Each file name starts with its pass and leg.
- **S14a part 1,** accepted: the step's decision list `directives/S14a_1_decision_list.md` (cleared at `e9eed614`),
  the record `steps/S14a_1_projection.md` (accepted at `4a4dfaf5`), and its record review dispositions
  `directives/_measurements/S14a_1_record_r{0,1}_review_disposition.md` (L4).
- **The plan:** `V3_STEP_PLAN.md`, the rows on S12, S14a, S14 and S15 (L5).
- **The register:** `SUBSTRATE_REQUIREMENTS.md`: its entry schema, its requirement/postulate/defect distinction, its
  population method including both routes, and `R-S12-02` (L6).
- **The S12 record** `steps/S12_order_source.md` (accepted at `ac909aaf`), for the owners of the forms (L6).
- `../../STATUS.md`.

## User decisions this list records
- **Framing** (2026-10-19) and **"Keep S12"** (2026-10-20), as the step's decision list carries them.
- **This part covers both remaining objects** (2026-10-21), and the two decisions **"Steady throughput"** and
  **"Exchange live"** (2026-10-21), as the step's decision list quotes them.
- **"Record as requirement"** (2026-10-21), after the step's decision list's second review pass. The user selected
  this option; the wording is the orchestrator's option text: "No build. Record that the slab's own dynamics (the
  open S12 forms) must pick the slow limit and the far behaviour, as a register requirement on the brane.
  S14/S15/S20 stay CONDITIONAL, now on that named requirement."

## What must be true
1. **What was attempted, and its standing.** The record states what the step's decision list named, in each
   revision, and what each review pass found and how it was disposed of. No build ran, by the user's decision. The
   record presents no computed object of S14a part 2 and no result of a build.
2. **The evidence, as filed.** Every statement about why the slow limit or the far behaviour was not determined rests
   on the review witnesses and the dispositions, as filed. Each statement cites its retrieval command, and the
   commands with their literal output go in the measurements file.
   - Only retrieval is allowed: existence, verbatim retrieval, literal-match counts, and the shape of a named stored
     object.
   - Each witness is stated with its pass, its leg and its standing: a review leg's instrument, not an accepted
     build.
   - Each witness is described by what its own source and its disposition say it constructs, checks and compares,
     and by no more (N4). Neither a general theorem nor the orchestrator's reading in the pass-2 disposition is
     presented as a computed result. That reading is quoted as a reading, with its source.
3. **The plan's objects.** The record states which of the plan's four S14a objects (`V3_STEP_PLAN.md:636–637`) are
   printed, and by which part, as its sources state them. It states that the remaining two are not derived, and the
   standing of S14, S15 and S20 as the plan gives it (`:639`). Nothing is claimed for them beyond what the user's
   decision states.
4. **What is handed on.** Under the register's own schema, rest-on test and method, both routes included, the pass
   hands on what the user's decision names: that the slab's own dynamics select the slow limit and the far
   behaviour.
   - The requirement names the object required, not a derivation route, law, material, premise or candidate.
   - It cites the record's statements and the filed evidence it rests on, and its consumers: S14, S15 and S20.
   - Its target follows the sources' owners (the S12 record for the forms, under "Keep S12").
   - It merges into an existing entry where one covers the same requirement, rather than duplicating it.
   - Where this list and the register disagree, the register wins, and the disagreement is reported.
5. **What remains open.** Each item the step's decision list retains as OPEN is stated with its owner, as the sources
   state them. The readings the user did not choose are listed as not computed. None is closed here.
6. **Status and provenance lines.**
   - S14a part 1's record and register pass were accepted at `4a4dfaf5` (L4). In the part 1 record, the register (its
     status and its S14a part 1 pass) and `../../STATUS.md`, every statement that presents them as not yet accepted,
     or their review as still to come, states that acceptance instead, heading or not. The lines found so far are in
     L7 and N5; that list is not claimed to be complete. Statements about S12's earlier state stay as they are.
   - The two notes carried from the part 1 record's round-1 review are repaired in the part 1 record: the l.59
     citation names its actual source, and the §6 quotation of the earlier population directive is retrieved from
     that directive, `directives/S11_requirements_pass2_directive.md` (L7).
   - Nothing else in the part 1 record or the register's S14a part 1 pass changes.
   - `../../STATUS.md` states where S14a part 2 stands, as the record does.
7. **The account of edits.** The record names every file this work changes outside the record, its measurements
   file and its generator, with each change. That account matches `git diff` against the commit the work starts
   from, and the measurements file retrieves that diff.

## Exclusions
- No new CAS, derivation, witness, closure or profile choice. No rerun of any witness.
- Edit only:
  - the record, its measurements file and generator;
  - `SUBSTRATE_REQUIREMENTS.md`;
  - `steps/S14a_1_projection.md`, for item 6 only;
  - `../../STATUS.md`.
- No commit and no push.

## STOP report (at most one page)
- The record's headline claims, each with its source.
- The evidence each claim rests on, with its pass, leg and standing.
- Register entry counts before and after, new entries by target, and existing entries that gained a source.
- Disagreements between sources, with both sides quoted.
- Open questions that need computation beyond retrieval.
- The account of edits (item 7).
