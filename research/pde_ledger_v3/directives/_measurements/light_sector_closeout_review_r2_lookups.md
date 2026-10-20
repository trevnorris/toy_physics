# Measurements: light-closeout review round 2 (generated 2026-10-19 19:32)

Generator: `_scratch/light_closeout/gen/review_r2_lookups.sh`. Literal retrieval only. The reviewed revision is
commit `53e0f2ea` plus the r2 working-tree repair, committed unchanged with this file.

## F1. The round: prompt, baseline, exits and verdict lines
```
$ sha256sum /var/projects/toy_physics/_scratch/light_closeout/review_prompt_r2.md | cut -c1-64
87ac5764a2e0ced63655844b0c7aacc4c0f82948cf8e849cee46dab83b3f54bf
```

```
$ cat /var/projects/toy_physics/_scratch/light_closeout/packet_baseline_r2.sha256 | sed 's#/var/projects/toy_physics/##'
4df5da71e12f855cd4724faf57c0172332caaefebdc1959798f3ea196d577544  research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md
50a4dcc54e85e786e64935f435fae28d57283eb9edbb00902ee924482cfcde04  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_lookups.md
1ae4b63efffd9f52c43a00a653a311059178e860ef3ed53b656d1950e320eec0  research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_outcome_lookups.md
b42dc759d6973ddfb562bc5d79a7796b06e8ff7815846c0ce9d436741a5dd968  research/pde_ledger_v3/directives/_measurements/light_sector_reassessment_lookups.md
313bc3d3452d6ff931eb7e2fdcf915b2a8068350d817b7cc088d4d11e0d3bb18  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r0_disposition.md
0067a5e204ea57bbe4f4a572442ae53375e2c3f92b4e411857e1f4b55596065f  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r0_lookups.md
3ee8bb0850b07419744308bb3dd711ce80a93d00c86b48b757eac075726440ea  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r1_disposition.md
6f5c6bc25b45ba07e3e55c3e6ab6e6a4bfdf4e6e825eedfa8b532fd84bccd296  research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r1_lookups.md
6b8b7ee0ee812953178682471e507d1f22deda7ddda2bb490bfb56b0dac5ad42  STATUS.md
a17567baafe68ef8571c80cb0735f58c2a0be97d896a1a0fb622890c2b1b2e76  research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
197252d6ee4549e599123eb086da8075949fa160ca9793fe87f1a0b754a9bee7  research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md
```

```
$ cat /var/projects/_review_scratch/light-closeout-r2-codex/DONE /var/projects/_review_scratch/light-closeout-r2-grok/DONE
EXIT=0
EXIT=0
```

```
$ grep -n '^CLEAR\|^\*\*CLEAR\*\*\|^FINDINGS\|^## FINDINGS' /var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt /var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt | cut -c1-200
/var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt:95:CLEAR
/var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt:88:**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
```

```
$ grep -n '^R2-[1-4] \|^## R2-[1-4]' /var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt /var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt | cut -c1-120
/var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt:15:R2-1 — HOLDS. The repaired R-S22-01 source 
/var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt:41:R2-2 — HOLDS. Packet quote (opening, lines 
/var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt:60:R2-3 — HOLDS. E15 prints the complete zero-
/var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt:67:R2-4 — HOLDS with explicit LC_ALL=en_US.utf
/var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt:3:## R2-1. Holds
/var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt:28:## R2-2. Holds
/var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt:62:## R2-3. Holds
/var/projects/_review_scratch/light-closeout-r2-grok/grok_review.txt:75:## R2-4. Holds
```

## F2. Both repositories and the packet unchanged by the legs
The first command drops this round's own disposition and lookups, written after both legs exited.
```
$ git -C /var/projects/toy_physics status --porcelain | grep -v 'light_sector_closeout_review_r2_' | diff /var/projects/toy_physics/_scratch/light_closeout/status_before_r2.txt - && echo IDENTICAL
IDENTICAL
```

```
$ git -C /var/projects/toy_physics-p2b status --porcelain | diff /var/projects/toy_physics/_scratch/light_closeout/p2b_status_before_r2.txt - && echo IDENTICAL
IDENTICAL
```

```
$ sha256sum -c /var/projects/toy_physics/_scratch/light_closeout/packet_baseline_r2.sha256 | sed 's#/var/projects/toy_physics/##'
research/pde_ledger_v3/steps/S9b_flow_oracle_comparison.md: OK
research/pde_ledger_v3/directives/_measurements/light_sector_closeout_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/S9b_flow_oracle_check_outcome_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/light_sector_reassessment_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r0_disposition.md: OK
research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r0_lookups.md: OK
research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r1_disposition.md: OK
research/pde_ledger_v3/directives/_measurements/light_sector_closeout_review_r1_lookups.md: OK
STATUS.md: OK
research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md: OK
research/pde_ledger_v3/DEFERRED_HEAVY_RUNS.md: OK
```

## F3. Codex's memory-named scratch scripts: the paths they open
```
$ grep -c '\.claude' /var/projects/_review_scratch/light-closeout-r2-codex/03_memory_source_pointers.py /var/projects/_review_scratch/light-closeout-r2-codex/14_memory_evidence.py
/var/projects/_review_scratch/light-closeout-r2-codex/03_memory_source_pointers.py:0
/var/projects/_review_scratch/light-closeout-r2-codex/14_memory_evidence.py:0
```

```
$ grep -n "Path('" /var/projects/_review_scratch/light-closeout-r2-codex/03_memory_source_pointers.py /var/projects/_review_scratch/light-closeout-r2-codex/14_memory_evidence.py
/var/projects/_review_scratch/light-closeout-r2-codex/03_memory_source_pointers.py:4:root = Path('/var/projects/toy_physics')
/var/projects/_review_scratch/light-closeout-r2-codex/14_memory_evidence.py:4:root = Path('/var/projects/toy_physics')
/var/projects/_review_scratch/light-closeout-r2-codex/14_memory_evidence.py:13:p2b = Path('/var/projects/toy_physics-p2b/_scratch/s11c')
```

## F4. Codex's locale note on reproducing the lookups
```
$ grep -n 'LC_ALL\|locale' /var/projects/_review_scratch/light-closeout-r2-codex/review_report.txt
11:E13–E16 are the primary checks; E11–E12 document the locale qualification.
67:R2-4 — HOLDS with explicit LC_ALL=en_US.utf8. E16 copies the generator without
79:Locale qualification: under the inherited LC_ALL=C.UTF-8, byte reproduction
82:The differences in the inherited-locale run are ordering only; they change no
```

```
$ grep -n 'LC_ALL\|sort' /var/projects/toy_physics/_scratch/s9b_oracle/gen/light_sector_closeout_lookups.sh | head -5
2:# Mechanical fact-lookups (grep/sed/awk field retrieval/sort/uniq/cut/wc/git log/git show only; no algebra) behind
25:  run "grep -n -F '\"residual_tag\" -> \"FLOW_SQUARED_RESIDUAL\", \"domain\" -> \"$d\", \"value\" -> True' $OUT | cut -d: -f1 | while read n; do awk -v n=\$n 'NR<n && /^FLOW_ORACLE_OPERAND/ {fo=\$0} NR<n && /^RESIDUAL_OPERANDS/ {op=\$0} NR==n {match(fo, /\"equations\" -> \\{[^}]*\\}/); print substr(fo, RSTART, RLENGTH) \" | \" op; exit}' $OUT; done | sort | uniq -c"
29:run "grep '^RESIDUAL_ZERO_TEST: <|\"residual_tag\" -> \"FLOW_SIGNED_RESIDUAL\"' $OUT | grep -o '\"value\" -> [A-Za-z]*|>' | sort | uniq -c"
31:run "grep '^RESIDUAL_ZERO_TEST: <|\"residual_tag\" -> \"EXCHANGE_RESIDUAL\"' $OUT | grep -o '\"value\" -> .*' | sort | uniq -c"
33:run "grep '^OPTICAL_IDENTIFICATION_MAP' $OUT | grep -o '\"oracle\" -> \"[^\"]*\"' | sort | uniq -c"
```

