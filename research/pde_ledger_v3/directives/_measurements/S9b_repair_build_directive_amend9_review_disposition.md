# S9b repair build directive, amendment 9: review dispositions (orchestrator)

**Artifact.** `directives/S9b_repair_build_directive.md` with amendment 9 (header paragraph, item 14, and a "Repair
round 6" line in Parts 2 and 3), orchestrator-written, for the build's repair round 6 after review r5
(`S9b_repair_build_r5_review_disposition.md`, row C1). Reviewed at sha256 `f6574006…`
(`_scratch/s9b_build/s9b_repair_build_directive_amend9_review_baseline_r0.sha256`). It changes what two ablation
harnesses print, so it is reviewed as physics until clear.

**Legs.** Orchestrator-written, so Codex (`gpt-6.1-sol`, xhigh) and Grok (`CLAUDE.md` G1), on the identical prompt
`_scratch/s9b_build/s9b_repair_build_directive_amend9_review_prompt_r0.md` (sha256 `99af0fac…`). Both reported before
I adjudicated. Commands and literal output: `S9b_repair_build_directive_amend9_review_disposition_lookups.md`
(generator `_scratch/s9b_build/gen/s9b_amend9_r0_lookups.sh`). The legs' scripts and stdout are copied to
`_scratch/s9b_build/s9b_amend9_review_r0_{codex,grok}_evidence/`.

## Round 0

- **Grok: CLEAR, no findings.** `_scratch/s9b_build/s9b_amend9_review_r0_grok.txt` (the report starts at
  "**Verdict: CLEAR**"). Its SymPy and Wolfram scripts, run through the guard, classify C1's cases under the amended
  text. The coefficient pair differs at a point where both are defined. The inert integral of zero evaluates to zero.
  A condition inside a set differs at a membership point. Two quantified statements with the same truth are the
  same value. On "admissible profile", it reads the class as the spec's and item 13's general radial functions, and
  finds no knife whose reading depends on a narrower class.
- **Codex: verdict withheld, then CLEAR, no findings.** Its first pass could not run its physics script. The guard
  refused it twice ("aggregate 16 GiB reservation budget exceeded"), because my S9b production runs held the main
  tree's 16 GiB budget at the time. Its first report (`s9b_amend9_review_r0_codex_final.txt`) withholds a verdict
  and establishes no finding. I resumed the same session after the production runs ended, with a status note only
  (`s9b_amend9_review_r0_codex_resume.md`); it carries nothing from the other leg. Its scripts then ran through the
  guard, and its report (`s9b_amend9_review_r0b_codex_final.txt`, full text in the evidence directory's
  `review_report.md`) is CLEAR. It reads "admissible" as bounded by the spec and the directive's general-radial-profile
  requirement, not by a family a builder chooses.

**Clearance.** Both legs CLEAR in round 0, each with its own scripts, so nothing outstanding changes what is
computed or may be claimed (G4). Amendment 9 is accepted.

**Process record.** The withheld first pass was an orchestrator scheduling error: a guarded job of mine ran in the
same checkout while a leg needed the guard. It cost one resume, and no finding.
