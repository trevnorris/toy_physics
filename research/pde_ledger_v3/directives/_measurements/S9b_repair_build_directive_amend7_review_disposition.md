# S9b repair build directive, amendment 7: review dispositions (orchestrator)

**Artifact.** Amendment 7 to `directives/S9b_repair_build_directive.md`, written by the orchestrator for the build's
repair round 4 (`S9b_repair_build_r3_review_disposition.md`, findings C1 and C2). It changes what the engines print,
so it is reviewed as physics until clear, by Codex and Grok (`CLAUDE.md` G1, orchestrator-written).

## Round 0

**Reviewed version.** The working-tree directive, sha256 `fd0d8931…`
(`_scratch/s9b_build/s9b_repair_build_directive_amend7_review_baseline_r0.sha256`), frozen as
`_scratch/s9b_build/S9b_repair_build_directive_amend7_reviewed_r0.md`. Identical prompt
`_scratch/s9b_build/s9b_repair_build_directive_amend7_review_prompt_r0.md` (sha256 `f090caad…`) to both legs; both
reported before I adjudicated.
- **Codex (`gpt-6.1-sol`, xhigh): CLEAR, no findings.** Final message
  `_scratch/s9b_build/s9b_amend7_review_r0_codex_final.txt`; full report, scripts and stdout copied to
  `_scratch/s9b_build/s9b_amend7_review_r0_codex_evidence/`.
- **Grok (`grok-4.7`): CLEAR, no findings.** Report `_scratch/s9b_build/s9b_amend7_review_r0_grok.txt` (the report
  starts at "**Verdict: CLEAR**"; the text before it is progress narration); script and stdout copied to
  `_scratch/s9b_build/s9b_amend7_review_r0_grok_evidence/`.

**Verification.** Mechanical lookups in `S9b_repair_build_directive_amend7_review_disposition_lookups.md`, generated
by `_scratch/s9b_build/gen/s9b_amend7_r0_lookups.sh`: the committed directive is the reviewed version, and both
verdicts are CLEAR with no findings.

**What the legs checked.** Both checked, each with its own script, that:
- item 5's new bullet names the radii the flyby and both round-trip legs traverse, with each turning radius fixed by
  that ray's construction; it is the domain of the spec's "Branch existence" paragraph, and neither engine is handed
  a root-solving recipe;
- the bullet holds together with item 6's every-`b` quantifier, the item 10 restrictions and forward copies, and K10;
- item 8 names one predicate by its meaning (true for path dependence, false for endpoint dependence) with the gate
  of observable 3, and states no value of it;
- the amendment carries no expected value, sign, limit or reviewer measurement, and drops nothing already required.

Both also found that the r3 disposition's G2 reason holds against item 14 for the cited K1 instance. Codex adds that
this does not certify future harness differences; item 14 binds the repaired harness in review r4.

**Codex's note on the repair route.** C1 also broke the spec's "everywhere along the ray" requirement, which already
bound the SymPy builder. The amendment states that requirement for both engines; it does not replace the spec
obligation. I agree. The r3 disposition's "comes from the directive" names why the SymPy builder missed it, not a
waiver.

**Clearance.** Nothing outstanding changes what is computed or what may be claimed (G4). Amendment 7 is accepted as
the round-0 version.

**Process.** Codex left an auxiliary Wolfram diagnostic running at its task cap (unit `s11c-guard-db76ab8a59db`),
reported it, and excluded it from its verdict. Its launcher was dead. I stopped the unit and reconciled the stale
reservation; the record is `_scratch/s11c/s9b-amend7-review-r0/independent-wolfram-mylKeM/RECONCILIATION.md`. Grok
ran with its cwd outside the repository, and its prompt told it that a turn ending without a tool call ends the run;
it finished in one session.
