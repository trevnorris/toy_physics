# Measurements: S15 reading build repair-4b resume note (restart interruption), leak and lane check

Literal-match counts and hashes only. The note is operational: it reports that the session was interrupted and points back at `resume_s15_reading_repair_4.md` (gated at `5878692b`). It adds no physics, no values, no methods and no conditions, so it changes nothing that is computed or may be claimed.

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_4b.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
b92971154b95d2b1a6b12dfb957d20cffd47f862dea277fd82d35eb1c3c14a10  _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
1069
```

```
$ /usr/bin/grep -c -i -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'beta' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'aa[' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'bb[' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'exponent' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'amplitude' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'condition' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Reduce' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Solve' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'SolveAlways' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
1
```

```
$ /usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
0
```

```
$ grep -c -F 'resume_s15_reading_repair_4.md' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md; grep -c -F 'Do not read the orchestrator' _scratch/s15/reading_build/resume_s15_reading_repair_4b.md
1
1
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; systemctl --user list-units --type=service --state=active --no-legend | grep -c s11c-guard; pgrep -c codex; cat _scratch/s11c/s15/repair3_90441_native_rate_factor_list_02/outcome.json | grep -o '"exitCode": *[-0-9]*'
{}

0
0
"exitCode": 1
```
