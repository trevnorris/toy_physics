# S15 comparator build review, round 0: dispositions (orchestrator)

The mechanical lookups behind each disposition are in `S15_comparator_build_review_r0_lookups.md`, beside this file
(generator `_scratch/s15/gen/s15_cmp_r0_lookups.sh`). "M-n" cites a section of that file. Every computed result quoted
here is a leg's saved stdout; the legs' scripts were not re-run.

## The round

- **Artifact.** These files, uncommitted on `f16725a9`:
  - the comparator `scripts/S15_cross_engine_comparator.py`;
  - its tests `scripts/test_S15_cross_engine_comparator.py`;
  - the run record `directives/_measurements/S15_comparator_run_record.md`;
  - the builder's final report, `directives/_measurements/S15_comparator_builder_report.md`, copied verbatim from its
    transcript;
  - the production output `_scratch/s15/cmp_work/full_07.jsonl`, with its audits;
  - the builder's change to the guard, `scripts/s11c_guarded_run.py`.

  All matched the baseline after the round. The working tree was unchanged apart from this file and the lookups file
  (M0).
- **Gate.** G1: a Codex-written script, so its legs are fresh Claude + Grok, reviewed until clear.
- **Prompt and legs.** `_scratch/s15/cmp_review/s15_cmp_build_review_prompt_r0.md`, sha256 `28ef471f…` (M0). Both legs
  received it identically, in parallel:
  - Grok (`grok-4.7`, effort high): `/var/projects/_review_scratch/s15-cmp-r0-grok/REVIEW.md`;
  - fresh Claude (Opus agent): `/var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md`.
- **Verdicts.** Grok: CLEAR. Claude: NEEDS REVISION, with two findings (M1).

## Findings

**C1 (Claude): the comparison-disable check switches off the whole comparison. Confirmed.**
- *List, item 7:* "Each control fails when the comparison step it targets is disabled."
- *Source:* the comparator has one disable switch, at line 1180, and it returns `not_compared` for every comparison
  (M2). The tests set it for each control in turn (lines 770–777). `detects()` first asserts that the baseline
  comparison is active (line 98). With the switch on, 15 of the control groups fail on that baseline assertion, and the
  rest fail because nothing is sampled (M2, the leg's guarded run `tests_base_disable_01`).
- *Effect:* the run record's "twenty-five comparison-disable groups rejected their reductions" (line 15, M2) shows only
  that a fully disabled comparison is noticed. It does not show that any control fails with the step it targets
  disabled, which is the property item 7 asks for. Grok reports the same global switch rejecting four groups, and did
  not examine which step each control targets.
- *What must be true:* each control fails when the comparison step it targets is disabled or altered, with the rest of
  the comparison left on. The run record states, for each control, the step and how it was disabled.

**C2 (Claude): the sparse-moment control's oracle uses the comparator's own moment function. Confirmed.**
- *Source:* `sparse_integral_oracle` is commented "Independent dense expansion" (test line 697). It builds its
  reference from `c.angular_moment` (line 706), the function the sparse path uses (comparator line 853) (M3).
- *Evidence:* the leg's ablation A10 halved the full-range moment (M3). Under it, all 26 synthetic groups passed,
  `sparse_integral_oracle` included. Production values moved in 9 of 17 selected bindings, while every sample residual
  stayed 0, because one evaluator serves both operands (M3).
- *The output stands.* The leg's own quadrature agrees with the printed integral values (M4). Grok's monomial-drop
  ablation moved the output (its report). So the moment table in the accepted code is not in question. What fails is
  the claim that the test checks it independently.
- *What must be true:* the oracle obtains its moments without the comparator's moment function, so that an altered
  moment fails the control.

## Examined and not filed

- **Ablations that move production and no test.** These are A01–A03 (surfaces swapped, drains exchanged, `Π_o`
  transposed in the tables), A06 (one stratum dropped), A14 (positional pairing) and A16 (a supplied equation applied
  to one operand inside the comparator).
  - The synthetic repoint tests build their fixtures from the loaded tables, so they test the comparison machinery,
    not the tables' content. Each table error drives the real-stream residuals nonzero.
  - Both legs' compare-first scripts found each selected object bound to the right ends (their reports).
  - The strata are printed in the output: 45 samples per binding, 3 forms × 5 strata × 3 points (M8).
  - A14 cannot move on these streams, whose key order is identical.
- **The test suite needs eight CPUs.** It asserts 8 CPUs and runs 8 workers. Under the guard's default one-CPU pooled
  job it refuses, so the legs ran their copies with a switch. The tests and the full run were made under the guard with
  `--all-cpus-exclusive` (M6).
- **Guard settings against the list.** Item 8 names `--pool s15-comparator --memory-gib 50`. The final runs used
  `--all-cpus-exclusive --memory-gib 55` (M6), under my resume notes 2 and 3: use the machine's cores, and a 55 GiB cap
  after a run was killed at 50 GiB. The deviation is mine and authorized.
- **The guard change (item 9).** Both legs found that it keeps the memory cap, zero swap, task cap and admission in
  every mode, and that a run without the new option behaves as before. Claude's comparison found two CPU-list spellings
  (`'6-6'`, `'06'`) now accepted. They name the same CPU set (M6). The builder's report does not mention the change (M6).
  It is reviewed here and is committed with the preserved baseline.
- **Binding `region_momentum_groups:1373`.** The comparator prints an exact residual `0` and both operands `0` at all 45
  samples (M8).
  - The Claude leg's quadrature of that integral gives `1.9e-07` and `1.0e-04` on its two configurations (M8). Its
    coarse and fine grids disagree, and the line reports `rel 1.05e+296` against zero. Its summary, "22 integrals × 2
    configurations agree … at ≤5.5e-15 relative" (M8), does not describe that line.
  - The leg did not print the integrand's scale, so its filed stdout cannot say whether these are rounding noise around
    an exact zero.
  - This is not a finding. The question of whether an operand that evaluates to zero at every sample is the printed
    object's value goes to round 1, for both legs to settle by computation.
- **A08 (integral compared by its integrand alone) went long-silent** in the Claude leg (stderr unchanged after 13:38).
  The leg finished that check from its test failures and partial samples. After the leg's report was in, I stopped the
  unit, `s11c-guard-1b2527c930bf`, at 16:17 (M7).

## Decision

**Not accepted at this version.**
- The comparator's computed output stands on both legs' evidence. C1 and C2 change what may be claimed for its
  controls (list item 7), and both are confirmed.
- This version, the guard change and the production output's hashes are committed as the reviewed baseline before
  repair. That commit preserves the version; it does not accept it.
- Repair goes to the builder: C1 and C2 as stated above. If the comparator's source changes, the complete comparison is
  rerun with the final code, as the brief says.
- Round 1 is then reviewed by fresh Claude + Grok. Its prompt adds the zero-operand question above.
