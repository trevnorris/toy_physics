# Measurements — S9b directive amendment 10, review round 0 (generated 2026-10-12 23:30)

Generator: `_scratch/s9b_build/gen/s9b_amend10_r0_lookups.sh`. The legs' scripts and stdout are copied (files under
2 MB) to `_scratch/s9b_build/s9b_amend10_review_r0_codex_evidence/` (from `/tmp/s9b_amend10_codex_review_r0_P9t3gy/`) and
`_scratch/s9b_build/s9b_amend10_review_r0_grok_evidence/` (from `/tmp/s9b-amend10-review-r0/`).

## The reviewed text is preserved byte-exact
```
$ sha256sum _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
8dbf60257ad5cc2aea35d469ba86d1631f29d9d817b0612a5614f2ab536e5c2c  _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
```

```
$ grep 'S9b_repair_build_directive.md' _scratch/s9b_build/s9b_repair_build_directive_amend10_review_baseline_r0.sha256
8dbf60257ad5cc2aea35d469ba86d1631f29d9d817b0612a5614f2ab536e5c2c  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

## Both verdicts
```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend10_review_r0_codex_final.txt
Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend10_review_r0_grok.txt
Verdict: NEEDS REVISION
```

```
$ cat _scratch/s9b_build/s9b_amend10_review_r0_codex.DONE _scratch/s9b_build/s9b_amend10_review_r0_grok.DONE
EXIT=0
EXIT=0
```

## The round-0 sentences the findings cite
```
$ grep -n 'by more than their printed error bound' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
291:      choice under item 13. When the two evaluations differ, exactly or by more than their printed error bound, the
```

```
$ grep -n "that meets the conditions the spec and the part's domain place on it" _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
288:      that meets the conditions the spec and the part's domain place on it. To show a change, the harness may
```

```
$ grep -n 'has looked for such a point and profile and found none' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
294:      has looked for such a point and profile and found none, its difference prints that it is undecided, so it
```

```
$ grep -n 'so an unchanged part prints as zero and a' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
280:      condition. Each part's difference is computed by the engine, so an unchanged part prints as zero and a
```

```
$ grep -n 'No profile family is used for a shared-vocabulary object' _scratch/s9b_build/preserved/S9b_repair_build_directive_amend10_r0.md
267:    - No profile family is used for a shared-vocabulary object. An engine may also print family-restricted
```

## The legs' evidence lines
```
$ grep -n 'activate_sameq_zero\|substituted_sameq_zero\|closed_evaluations_equal\|matrix_p_at_3_2\|matrix_q_at_3_2\|exists_profile_meeting_both' _scratch/s9b_build/s9b_amend10_review_r0_grok_evidence/cases_executed.wl.stdout
4:activate_sameq_zero True
6:substituted_sameq_zero False
7:substituted_activate_sameq_zero True
28:closed_evaluations_equal True
29:matrix_p_at_3_2 True
30:matrix_q_at_3_2 False
32:exists_profile_meeting_both False
```

```
$ ls _scratch/s9b_build/s9b_amend10_review_r0_codex_evidence/
assemble_report.py
baseline_final.stderr
baseline_final.stdout
baseline.stderr
baseline.stdout
check_amendment.py
check_equalities.wl
collect_excerpts.py
collect_files.py
diff.stderr
diff.stdout
directive.stderr
directive.stdout
disposition.stderr
disposition.stdout
file_commands.json
file_evidence.md
finish_checks.command.json
finish_checks.py
guard.stderr
guard.stdout
lookups.stderr
lookups.stdout
py_engine.stderr
py_engine.stdout
py_harness.stderr
py_harness.stdout
report_body.md
report.md
run_guard.py
spec_diff.stderr
spec_diff.stdout
spec.stderr
spec.stdout
symbolic.attempts.jsonl
symbolic.launcher-1.stderr
symbolic.launcher-1.stdout
symbolic.launcher-2.stderr
symbolic.launcher-2.stdout
symbolic.launcher-3.stderr
symbolic.launcher-3.stdout
symbolic.launcher-4.stderr
symbolic.launcher-4.stdout
symbolic.launcher-5.stderr
symbolic.launcher-5.stdout
symbolic.launcher-6.stderr
symbolic.launcher-6.stdout
symbolic.launcher-7.stderr
symbolic.launcher-7.stdout
symbolic.outcome.json
symbolic.stderr
symbolic.stdout
wl_engine.stderr
wl_engine.stdout
wl_harness.stderr
wl_harness.stdout
wolfram.attempts.jsonl
wolfram_fixed.attempts.jsonl
wolfram_fixed.launcher-1.stderr
wolfram_fixed.launcher-1.stdout
wolfram_fixed.outcome.json
wolfram_fixed.stderr
wolfram_fixed.stdout
wolfram.launcher-1.stderr
wolfram.launcher-1.stdout
wolfram.outcome.json
wolfram.stderr
wolfram.stdout
```

## The round-1 text
```
$ grep -n 'prints by this rule, in this order' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
289:      Each part's difference is computed by the engine and prints by this rule, in this order:
```

```
$ grep -n 'compared as the condition it states on its free variables' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
287:      compared as the condition it states on its free variables, so its point gives values to those only.
```

```
$ grep -n 'a bound that covers the error of both' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
295:         a bound that covers the error of both. The change prints with the point, the profile when there is one,
```

```
$ grep -n 'It need not meet the condition being compared' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
302:      derivative the part takes, at the witness point. It need not meet the condition being compared, nor the
```

```
$ grep -n 'neither a profile choice nor a profile family' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
304:      this test only. They enter no engine object, and they are neither a profile choice nor a profile family
```

