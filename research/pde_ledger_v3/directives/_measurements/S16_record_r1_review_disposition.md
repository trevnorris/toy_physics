# S16 record review, round 1: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_record_r1_review_lookups.md`, beside this file (generator
`_scratch/s16/gen/rec_r1_lookups.sh`). "Z-n" cites a section of that file.

## The round

- **Artifact.** The record author's repair of the round 0 version, which is preserved at `4ae7e6cf`. The repair is
  uncommitted on top of that commit. It changed the record, its measurements file and generator, the register and
  `../../STATUS.md`; `steps/S14a_2_slow_limit.md` is unchanged (Z0).
- **Prompt and legs.** `_scratch/s16/rec/s16_record_review_prompt_r1.md` (Z0). It is the round 0 prompt, plus an item 0
  that tests Y1–Y4 and the diff against `4ae7e6cf`, with the round 0 dispositions handed. Codex wrote the artifact, so
  the legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`, high) (G1). They had the identical prompt and ran in
  parallel. Grok exited 0.
- **Baseline.** All six baseline files matched after the round. The working tree was unchanged, apart from the lookups
  file this disposition adds (Z0).
- **Verdicts.** Claude: NEEDS REVISION, one finding. Grok: CLEAR (Z0).
- **Both legs** report that Y1–Y4 now hold. Both reran the generator and reproduced the filed measurements file byte
  for byte.

## Findings

**W1. The list author's filings are admitted on the record list's direction, not on the register's test, and the
disagreement is not reported.** (Claude 1.) **CONFIRMED.**
- The register's method makes each "*this is deferred to*" statement an entry (Z1). Its earlier passes applied it more
  narrowly (Z1):
  - pass 1: "A requirement is not everything a later step would find useful. … ⛔ not an obligation on which a banked
    result rests. It got no entry.";
  - pass 2: "an OPEN item sources no entry under the rest-on test";
  - the O2 pass: "An unspecified response form needed only to select or solve a later closed material state is an
    OPEN handoff … A future closure owner alone does not satisfy this test." Its row for the density response, the O2
    `:138` object, reads "No new density-response entry".
  - The one forward object the register has admitted is S14a part 2's, which "is explicitly required by the user, not
    entered merely because a future solve could use it."
- The S16 pass admits compactness (`R-S16-02`), face momentum (`-03`) and density response (`-05`) on "Record contract
  item 7 directs these named filings". `R-S22-02`'s admission ground also cites item 7 (Z2).
  - `R-S16-05` is the O2 `:138` object the O2 pass declined. Its entry says that decision "remains unchanged" (Z2).
  - After Y3, none of the four has a user requirement behind it (round 0 disposition, Y3).
- The record list, item 9: "Where this list and the register disagree, the register wins, and the disagreement is
  reported" (Z3). Neither the record nor the pass reports such a disagreement for these filings. The only "register
  wins" in either is about the population route, at register l.1782 (Z2).
- **The other side,** which the leg also states:
  - The method's "Each becomes an entry keyed by the step that must deliver" can be read to admit them.
  - The plan puts compactness among the "supplied inputs that belong in S22's debt table".
  - Pass-2 finding E3 has the record file the body's smallness as its own requirement.

  So the register text alone does not settle it. The pass rests the admission on item 7 and does not take this up.
- **Cause: the record list.** Item 7 directs the record to file each of the step list's requirements "as that section
  states" (Z3). Item 9 makes the register's test govern. For these four filings the two items point different ways,
  and the record followed item 7. The orchestrator wrote both items.
- **What must be true.** Each list-author filing (`R-S16-02`, `-03`, `-05`, `R-S22-02`) is admitted by the register's
  rest-on test, as its earlier passes applied it, or by an explicit user requirement.
  - Where neither admits a filing, it is handed on as an OPEN item, not an entry. The disagreement with the list is
    reported with both sides quoted.
  - The census, the record and `../../STATUS.md` follow from the result.

## Notes below the filter (Claude)

- **The step list's `:20` sub-heading,** "Where the user's decisions differ from that sentence", also groups the
  compactness filing under the user's decisions (Z4). It is a second instance of the heading attribution in Y3. The
  entries already name the right filer.
- **The U1 row in record §1** quotes S12's words as part 1 quotes them. The standing is unchanged.
- **The leg's generator rerun** created or refreshed the git-ignored
  `scripts/__pycache__/S16_record_measurements.cpython-312.pyc`. No tracked file changed.

## Outcome

One finding is confirmed. It decides whether up to four of the six new entries are register entries or OPEN items. That
changes the census, the by-target rows and `../../STATUS.md`. Whether a forward need with no banked result resting on
it goes on the register has so far been the user's decision (S14a part 2). So the user is asked whether to require
these filings. The reviewed version is preserved with these dispositions. The record's author repairs it to what must
be true, given the user's answer. The repair gets a round 2 review by a fresh Claude agent and Grok.
