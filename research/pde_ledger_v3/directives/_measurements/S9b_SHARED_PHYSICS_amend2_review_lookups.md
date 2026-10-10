# Measurements — S9b spec amendment 2 and directive amendment 5, review round 0 (generated 2026-10-09 18:55)

Generator: `_scratch/s9b_build/gen/s9b_spec_amend2_r0_lookups.sh` (sha256sum/grep only).

## The reviewed versions are the working-tree files
```
$ sha256sum -c _scratch/s9b_build/s9b_spec_amend2_review_baseline_r0.sha256
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

```
$ sha256sum _scratch/s9b_build/S9b_SHARED_PHYSICS_amend2_reviewed_r0.md _scratch/s9b_build/S9b_repair_build_directive_amend5_reviewed_r0.md
f4c9d3516a50c091a56ef4dd09366b9213b0c3a38e463daf0f43dbbaf550daaa  _scratch/s9b_build/S9b_SHARED_PHYSICS_amend2_reviewed_r0.md
b3274adb086b46f2965ff33a92682da4253b3c71fa83ce1c13f6aff0d5102653  _scratch/s9b_build/S9b_repair_build_directive_amend5_reviewed_r0.md
```

## Both verdicts
```
$ grep -n -o 'Verdict:\*\* [A-Z]*\|Verdict: [A-Z]*' _scratch/s9b_build/s9b_spec_amend2_review_r0_codex_final.txt _scratch/s9b_build/s9b_spec_amend2_review_r0_grok.txt
_scratch/s9b_build/s9b_spec_amend2_review_r0_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_spec_amend2_review_r0_grok.txt:1:Verdict: CLEAR
```

```
$ grep -n '^\*\*Findings' _scratch/s9b_build/s9b_spec_amend2_review_r0_codex_evidence/report.md _scratch/s9b_build/s9b_spec_amend2_review_r0_grok.txt
_scratch/s9b_build/s9b_spec_amend2_review_r0_codex_evidence/report.md:3:**Findings: None within the requested amendment scope.** This clears the governing text, not the repaired engines.
_scratch/s9b_build/s9b_spec_amend2_review_r0_grok.txt:29:**Findings:** none.
```

## The r1 disposition's knife sentence is not in the directive
```
$ grep -n 'A knife that changes a condition' research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r1_review_disposition.md
37:| C1 | The implied `j_n` is not computed from its condition. Each `*_IMPLIED_JN_*` payload pairs the reduced condition with the supplied balance under the row's profile replacement. `V` is left as a free profile, so the `j_n` part is the same whatever the condition says. All 18 tags, both engines. (Claude) | **ACCEPT, both en
```

```
$ grep -c 'A knife that changes a condition' research/pde_ledger_v3/directives/S9b_repair_build_directive.md
0
```

```
$ grep -n -o 'the disposition.s blanket knife statement is not a physics requirement' _scratch/s9b_build/s9b_spec_amend2_review_r0_codex_evidence/report.md
17:the disposition's blanket knife statement is not a physics requirement
```

