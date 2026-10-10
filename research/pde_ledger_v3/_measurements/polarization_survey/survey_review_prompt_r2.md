# Review: the polarization survey, repair round 2 (Codex-written)

## Artifact
`/var/projects/toy_physics/_scratch/polarization/POLARIZATION_SURVEY.md`. Its sha256 is in
`/var/projects/toy_physics/_scratch/polarization/survey_review_baseline_r2.sha256`. Codex wrote it from the prompt
`/var/projects/toy_physics/_scratch/polarization/survey_prompt_r0.md`; read the prompt to see what was asked.

This is the survey after two repairs. The previous rounds' dispositions, with what must be true after each repair,
are in `/var/projects/toy_physics/_scratch/polarization/survey_r0_review_disposition.md` and
`survey_r1_review_disposition.md` in the same directory. The reviewed previous versions are
`POLARIZATION_SURVEY_reviewed_r0.md` and `POLARIZATION_SURVEY_reviewed_r1.md`. The class definitions now in force are
in the "The classes" section of `survey_repair2_prompt.md`; they replace the task's class list. Check that each "what must be true"
holds, and that the repair introduced nothing wrong. Review the whole document, not only the changes.

The survey is for the user, a programmer and not a physicist. It explains light's polarization, lists what
experiments show, inventories what this repository says, and classifies each experimental fact against the toy
model. In the model, light is a transverse shear wave of a brane: an ordered, finite-thickness slab centred at
`w = 0` in a four-dimensional superfluid bulk. The classification will be used to decide what the model must explain,
so an error in it changes what may be claimed.

## What to check
1. **Part 1, the explainer.** Is each statement correct? Is each one correctly labelled as a classical-wave fact or
   a quantum fact? Would a careful non-physicist come away with a wrong picture anywhere?
2. **Part 2, the experimental record.** For each entry, open the cited source where you can. Check:
   - that the source says what the survey says;
   - the year, the authors, and the precision or bound;
   - that unestablished signals are marked as such.

   Name any load-bearing experiment or bound that is missing.
3. **Part 3, the repository.** Check the quoted `file:line` sources and the status each one is given. Name any
   source in the v3 ledger (`research/pde_ledger_v3/`) that bears on polarization, the number of transverse modes, or
   the brane's displacement directions, and that the survey missed or misread.
4. **Part 4, the classification.** Are the class definitions in force clear and mutually exclusive, and do they
   sort the facts in a way the user can act on? Is each row in the class those definitions give? In particular:
   - Is anything marked "not addressed" or "required" actually in conflict with a model source?
   - Is anything marked "reproduced" stronger than its source supports?
5. **Scope.** Does the survey propose mechanisms or resolve conflicts, which it was asked not to do?

## Method
- Every claim about a repository file needs the command you ran and its literal output.
- Every claim about an external source needs the source: a URL, DOI, or file you opened, and the quoted passage.
  Say whether you opened it or relied on a secondary source.
- If a finding rests on a calculation, write a script, run it, and save the script and its literal stdout to named
  absolute paths under `/tmp`. A prose derivation is worth nothing.
- Quote both sides for every finding: the survey's text and the source's.
- Do not modify any file under `/var/projects/toy_physics`. Do not read `~/.claude/projects/`. Do not read the
  review reports of either round.

## Physics filter
Report a finding only if it would change what the user understands, what an experiment is said to show, or how a
fact is classified against the model. Do not report wording or formatting that changes none of these.

## Output
Your final message is your report:
- **Verdict:** CLEAR or NEEDS REVISION.
- **Findings,** numbered. For each: the location in the survey, what is wrong, the evidence, and the fix.
- A short list of what you checked and found sound.

Keep it under 1,500 words. Write your report and exit. Do not spawn agents.
