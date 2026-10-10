# S9b spec amendment 2 and build directive amendment 5: review dispositions (orchestrator)

**Artifacts.** Two orchestrator-written amendments, reviewed together:
- amendment 2 to `directives/S9b_SHARED_PHYSICS.md`, on top of `ede8aa21`: `c_γ` is the non-negative root of `c_γ²`,
  and branch existence is stated through `c_γ²` without presupposing its sign;
- amendment 5 to `directives/S9b_repair_build_directive.md`, on top of `ede8aa21`: items 5, 6 and 9, and one pointer
  line in each of Parts 2 and 3.

**Why.** Build review r1 accepted four findings
(`S9b_repair_build_r1_review_disposition.md`; engines preserved at `f4dc2d3d`): the implied `j_n` was not computed
from its condition; Wolfram's Part C domains carried an unsupplied amplitude cut; branch existence presupposed the
sign of `c_γ²` and the spec left `c_γ`'s root open; the `γ` difference was printed unreduced. The amendments state
what must be true after the repair. They change what the engines compute, so they were reviewed as physics until
clear, with Codex and Grok (G1, orchestrator-written). This is a technical repair under the user's standing approval
(2026-10-09).

Commands and literal output: `S9b_SHARED_PHYSICS_amend2_review_lookups.md`. Generator:
`_scratch/s9b_build/gen/s9b_spec_amend2_r0_lookups.sh`.

## Round 0

Both legs used the identical prompt `_scratch/s9b_build/s9b_spec_amend2_review_prompt_r0.md`. They reviewed the
working-tree versions at `_scratch/s9b_build/s9b_spec_amend2_review_baseline_r0.sha256`: spec `f4c9d351…`, directive
`b3274adb…`. Frozen copies: `_scratch/s9b_build/S9b_SHARED_PHYSICS_amend2_reviewed_r0.md` and
`…/S9b_repair_build_directive_amend5_reviewed_r0.md`.
- **Codex** (gpt-6.1-sol, xhigh): CLEAR, no findings. Final report `_scratch/s9b_build/s9b_spec_amend2_review_r0_codex_final.txt`;
  full report and derivations `…_r0_codex_evidence/` (`report.md`, `derive_review_v2.py`, `derive_root_live.py`
  and their stdout).
- **Grok** (grok-4.7): CLEAR, no findings. Report `_scratch/s9b_build/s9b_spec_amend2_review_r0_grok.txt`;
  derivation `…_r0_grok_evidence/derive_s9b_amend2.py` and its stdout.

Both reported before I adjudicated. Both checked that a repair meeting the text, and no weaker repair, closes each
of the four findings.

**Recorded, not findings.**
- **The r1 disposition's knife sentence.** My r1 disposition listed "A knife that changes a condition changes its
  implied `j_n`" among what must be true. Codex notes this does not hold as a blanket rule: where a restriction
  fixes `V ≡ 0`, a changed condition need not change the implied `j_n`. The sentence was not carried into the
  directive (lookup: zero occurrences). The build review checks dependence by FORM ablation instead.
- **The root and the optical chart.** Both legs find that the non-negative root agrees with the chart the optical
  counting already used near `δ = 0`. The root fixes the speed convention, not the sign of the live stiffness ratio
  or of `V`.
- **The legs' constructions,** including their witness profiles, the branch structure they derived for the implied
  `j_n`, and their reduced form of the `γ` difference, are reviewer measurements. They stay out of builder packets
  and leg prompts.

**Outcome.** Nothing outstanding changes what is computed or what may be claimed. Both amendments are accepted at the
reviewed versions. Next: both builders repair their engines (repair round 2), then build review r2.
