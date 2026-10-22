# Gravity prior-art section, review round 2 (scoped): dispositions (orchestrator)

The mechanical lookups behind each disposition are in `gravity_prior_art_review_r2_lookups.md`, beside this file
(generator `_scratch/prior_art_gravity/gen_r2_lookups.sh`). "P-n" cites a section of that file.

## The round

- **Scope.** After round 1 the user chose "Repair, then scoped review (Recommended)" (10-22). The round reviews
  repair 2 only: S1 from round 1, and whether anything else changed.
- **Artifact.** The section after repair 2, which the authoring subagent made. The round 1 version is preserved at
  `d469f0dc`.
  - The repair adds Wang's presupposition (4), no external body forces, to G7's conditions, with the results that
    assume it (P3).
  - In the sources file, it extends the Theorem 3 quotation to the sentence that inherits that presupposition (P3).
  - Its diff against `d469f0dc` is `_scratch/prior_art_gravity/doc_edit_r2.diff`: 3 lines added, 2 removed (P3).
- **Prompt and legs.** `_scratch/prior_art_gravity/review_prompt_r2.md` (P1). Codex (`gpt-6.1-sol`, xhigh) and Grok
  (`grok-4.7`, high), in parallel, with web access. Both exited 0. Both files matched the baseline after the round,
  and the working tree was unchanged apart from the lookups file (P1).
- **Verdicts.** Codex: CLEAR. Grok: CLEAR (P2).
  - Each opened Wang's paper (`arXiv:physics/0506062`, v5).
  - Each found presupposition (4) in Theorem 2, and Theorem 3 and Corollary 4 assuming presuppositions (1)–(5).
  - Each found the extended quotation accurate, and the repair the only change against `d469f0dc`.

## Disposition: cleared

Nothing outstanding changes what the section may attribute or claim. The gravity attribution section, its
cross-reference and its sources file are accepted at the commit that adds this file.

## Carried

- The legs could open some works only as abstracts or bibliographic records (rounds 0 and 1). The section and its
  sources file mark what was opened for each citation.
- The S16 record cites this section's rows as the section labels them. It checks its computed objects against a prior
  result only where that row's conditions match S16's.

## The section's review history

| Round | Codex | Grok | Outcome |
|---|---|---|---|
| 0 | 4 findings | 5 findings | R1–R7 confirmed; preserved `eb5e95db`; repair 1 |
| 1 | 1 finding | CLEAR | S1 confirmed; preserved `d469f0dc`; user: repair, then scoped review |
| 2 (scoped) | CLEAR | CLEAR | accepted |
