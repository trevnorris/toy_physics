# S16 comparator scoped review, round 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S16_comparator_review_r2_lookups.md`, beside this file
(generator `_scratch/s16/cmp_review/gen_r2_lookups.sh`). "R-n" cites a section of that file. Rounds 0 and 1 are in
`S16_comparator_review_r{0,1}_disposition.md`.

## The round

- **Artifact.** The comparator, its tests, the run record and the builder's report, after repair 2. Codex wrote repair
  2 in its repair 1 session, resumed. Its prompt is `_scratch/s16/cmp_repair2/resume_prompt_repair2.md`. It points at
  K5's "What must be true" and extends it to every SymPy citation the comparator prints. It asks for nothing else.
  - The production output is `_scratch/s16/cmp_work/repair_r2_full_01.jsonl`, with the test output
    `repair_r2_tests_01.json` and the builder's own citation audit `repair_r2_audit_01.json`.
  - The builder's report is the repair's final message, copied verbatim from its transcript.
- **Scope.** The user chose a citation-only repair with a scoped check (round 1 disposition, "Outcome").
- **Prompt and legs.** `_scratch/s16/cmp_review/s16_cmp_build_review_prompt_r2.md` (R1). It covers four things:
  - every printed citation, in both engines, against the engines at `1d990ee6`;
  - that the code diff against `de28e059` touches only citations;
  - that the output differs from `repair_r1_full_01.jsonl` only in citations and runtime fields;
  - the run record.

  The legs are a fresh Claude agent (Opus) and Grok (`grok-4.7`, high) (G1). They had the identical prompt and ran in
  parallel. Grok exited 0.
- **Baseline.** All seven baseline files matched after the round. The working tree was unchanged (R1).
- **Verdicts.** Claude: CLEAR, no findings. Grok: NEEDS REVISION, three findings (R2).
- **Method deviation.** Claude reports one: it read the comparator before writing its engine index. Its index script
  reads only the two engines, and its containment check rests on that index.

## Scoped checks

- **K5 is repaired.** The rules round 1 named now cite their build or print lines (R4):
  - `definitions` 201–204, with the names at 101–103;
  - `momentumReceived` and `momentumRate` 211–214;
  - `momentum_surface_motion` 215–216;
  - `force_companion` 193;
  - `partition_routes` 289;
  - `coverage` 435–437;
  - `dimensions` 359–407 and 425–428.
  
  The three wrong reason strings no longer occur in the output (R4).
- **No citation misses its object, and none leaves its file** (Claude, R5):
  - the highest cited lines are SymPy 437 of 439 and Wolfram 246 of 249;
  - no cited line is outside its file;
  - none of the 94 rule citations lacks a print or direct-build line for a tag it covers;
  - every name-table citation contains the line that introduces the name.

  Grok's check finds no cited line out of range (R6).
- **The code changed only in citations** (R3).
  - Grok's token comparison against `de28e059`: 38 string literals changed, all of them citations. No other token
    changed, and no comment changed.
  - The tests are identical to `de28e059`, and the engines to `1d990ee6`.
- **The output changed only in citations and runtime fields** (R3). Grok's comparison of all records against
  `repair_r1_full_01.jsonl`: citation spans and line specs differ, the two runtime and memory fields differ, and
  nothing else does. Claude finds the same, and that the printed tables regenerate exactly from the final code.
- **The run record** matches the guard records and the files (both legs).

## Grok's findings: ranges that also include other objects' lines

**Grok 1–3. Nineteen citations name a range that also includes lines building a neighbouring object.**
**FACTS CONFIRMED; NOT A FINDING.**
- Grok's check marks a citation as failing when it names a code line outside Grok's own block for the object. The
  failing kinds are "foreign" lines, never "out of range", and the check has no case for a citation that misses its
  object (R6).
- Of the 19, 13 are Wolfram citations. Those were not part of K5 and have been unchanged since round 0. The other 6
  are SymPy: the rate and alias ranges, and three reason strings (R6).
- An example is the Wolfram rates range, 99–125 (R6). It contains the rate construction, its emits and the calls. It
  also contains line 113, the force-companion emit inside the same function, and the moment and definition lines
  115–122.
- Claude saw the same ranges, and recorded them as notes, not findings.
- **Why not a finding.**
  - The ranges change nothing the comparator computes (R3).
  - They change nothing that may be claimed from them. Each citation contains the lines that build or print its
    object, inside the engine (R5). That is what decision-list items 2 and 3 ask for: each cites "the line in each
    engine that builds it" or "introduces it". It is also what K5 required, and what it found violated: citations
    that missed their object, or left the file.
  - The bindings themselves rest on the round 1 measurements: the independent comparisons, the form ablations and
    the one-sided corruption.
- **What may be claimed.** Each printed citation contains the lines that build or print its object, inside the engine
  as accepted at `1d990ee6`.
- **What may not be claimed.** That every line a citation names belongs to that object. Some ranges also span
  neighbouring objects' lines.

## Not findings

- **Claude: six rules print their first row's object text as the rule's object.** One example is `dimensions`,
  "item 5 dimension of supplied material equation". This is unchanged since round 1, and each binding row prints its
  own correct object text. It changes nothing computed. A record cites a binding's own object text, not the rule's.

## Outcome

Nothing outstanding changes what the comparator computes or what may be claimed from its output or its citations,
with the claim stated above. K1–K5 are repaired.

This was the third review round, and one leg's verdict was NEEDS REVISION. So the finish was put to the user, with
two options:
- accept now, with the claim stated as containment;
- trim every range to its object's lines, then check again.

The user chose to accept (2026-10-24).
- **Accepted.** The production output is this round's reviewed run, `repair_r2_full_01`, itself, not a rerun. It is
  filed as `research/pde_ledger_v3/scripts/out/S16_cross_engine_comparator.out`. Its identity and provenance are in
  `S16_comparator_production_run.md`.
- **Carried to the S16 step record.**
  - A citation contains its object's lines, and a range may also span neighbouring objects' lines.
  - The bindings rest on the round 1 measurements, not on the citations.
  - Six rules print their first row's object text; cite each binding's own object text.
