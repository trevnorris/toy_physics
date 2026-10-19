# S9b record, review round 2: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S9b_record_r2_review_lookups.md` beside this file (generator
`_scratch/s9b_build/gen/s9b_record_r2_lookups.sh`).

## The round
- **Artifacts.** The record, its measurements file and generator, and the S9b register pass, as repaired in round 2
  by the same author (brief `_scratch/s9b_build/record/resume_record_r2.md`, against the round-1 disposition's "What
  must be true"). The reviewed previous version is `d1511b20`. Four files in the baseline
  `_scratch/s9b_build/record/s9b_record_review_baseline_r2.sha256`, all OK at review.
- **Prompt and legs.** The prompt is `_scratch/s9b_build/record/s9b_record_review_prompt_r2.md` (sha256 `ac430675…`):
  round 1's prompt with the baseline and the previous version advanced. The legs were a fresh Claude agent (opus) and
  Grok (`grok-4.7`, effort high), launched identically, with scratch under
  `/var/projects/_review_scratch/s9b-record-r2-{claude,grok}/`. Reports:
  `_scratch/s9b_build/record/s9b_record_review_r2_claude.md` and `…_r2_grok.txt`. Grok returned **CLEAR**; Claude
  returned **NEEDS REVISION** with two findings.
- **Review scope.** This is the third review of the record. The user chose a full review for it, from the options
  full, scoped or stop.
- **A user decision this round.** Asked whether a check of the flow-only condition against the river/acoustic-flow
  prior art should be queued, the user chose "Queue it after record". The check is new computation, done after the
  record is accepted and outside it: Codex writes the check and it is reviewed under G1.

## Findings

| # | Leg | Finding | Disposition |
|---|---|---|---|
| 1 | Claude 1 | **The unperformed oracle applicability checks are not handed on as open items with an owner.** The record applies no oracle (692) and says "Applicability checks of these candidate oracles remain unperformed, not negative results" (709–710). The OPEN list (736–744) and the §7 owner table carry no such item. The register's method says "Run this route for every sector with identified prior art" (835, 1099). The decision list excludes an oracle check from the record (152). | **Confirmed.** The lookups show no oracle or applicability item in §7 (602–662) or in the OPEN list (736–745). An oracle the record names but does not apply is an M3 obligation, and as written it is handed on to no one. Under the user's decision above, the river/acoustic-flow check now has an owner and a place: after the record, outside it. So the decision list's exclusion and the register's route do not conflict. The record runs no oracle; the check runs after it. The leg's request to report a source tension therefore becomes a request to report the decision. |
| 2 | Claude 2 | **The comparator's printed stratum inventory is not reported.** The comparator prints 660 `stratum_inventory` lines, one of its six line kinds. 658 of them carry the coverage label "unresolved coverage unless sampled in paired comparison", and the leg counts 19,534 strata and 2 empty inventories. The record, the measurements file, its generator and the register never mention the line kind or the label. The record says "Every remaining coverage reason and every affected object is retained in M7, verbatim from the accepted coverage appendix" (553–554). | **Confirmed.** The lookups count 660 lines in the comparator output, 658 with that label, and zero mentions of either in the four files. The sentence at 553–554 is true of the accepted coverage appendix, but not of the comparator's printed output. A coverage class of the comparator is not reported, so what the record hands on as uncovered is incomplete. Grok's enumeration counted the same 660 lines (`recount_coverage_stdout.txt`) and did not report it. The comparator's dispositions mention the inventory only twice, both in r4: its run-to-run order (item c), and as the place one reader ablation's change showed (item b). |

**Sound, both legs, each with its own scripts and stdout.**
- Both legs reran the generator, and it reproduces the filed measurements file byte for byte (sha256 `ff14384c…`).
- Both legs' comparator enumerations match the record's row classes, sample counts, `not_formed` reasons, symbolic
  classes and Appendices B and C.
- Round 1's repairs hold:
  - The item-5 count follows the printed γ-difference objects, and the Wolfram field-equivalence question stays open
    on its own.
  - Each mechanism alone is reported with both engines' printed objects.
  - The engine-specific branch bindings and same-named row differences are listed.
  - S9b is a source on `R-S8-03`, and the census follows the rest-on test.
- The forward case is labelled conditional. The premises, D5's obligation and the dependences are carried. No
  expected value, sign or density-link effect is stated.
- The Claude leg reports 14 paragraphs after `R-S21-01`'s last schema field, each copied verbatim from §3. It
  considered that and did not report it, since nothing handed on changes. I agree.

**Not clear.** Both findings change what the record hands on. Neither is a defect the round-2 repair introduced.
Finding 2 has been there since round 0 and five earlier legs passed it. Finding 1 is in the second-route text that
round 1 rewrote, and the user's decision this round resolves it. The reviewed round-2 record is preserved before any
repair overwrites it.

## What must be true after the repair
1. Every identified prior-art oracle whose applicability check is unperformed is handed on as an open item, in the
   record and in the S9b register pass. Each has its owner, or "no owner named" where there is none. The
   river/acoustic-flow check carries the user's decision above. The second-route question about each such oracle
   remains open with it. No oracle is applied, and no result of one is stated.
2. The record reports every line kind the comparator prints, including the stratum inventory, with each kind's
   printed coverage label, its size and the objects it leaves uncompared. Each of these rests on a retrieval in the
   measurements file. The record says what M7 does and does not contain.
3. Nothing else changes. The measurements file is regenerated by the final generator.
