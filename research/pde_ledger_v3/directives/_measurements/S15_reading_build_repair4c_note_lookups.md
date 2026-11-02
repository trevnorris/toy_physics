# Measurements: S15 reading build repair-4c resume note (memory stop), leak and lane check

Literal-match counts, hashes and stored facts only. The note is operational: it reports the builder's memory stop, raises the guard memory cap from 50 to 55 GiB, and points back at `resume_s15_reading_repair_4.md` (gated at `5878692b`). It adds no physics, no values, no methods and no conditions, so it changes nothing that is computed or may be claimed.

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_4c.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
6c8e1a3d1179aaf1e9b395e1be45199b84da1d1ad1e62027d5fb05fcb55d4b62  _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
1315
```

```
$ /usr/bin/grep -c -i -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'beta' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'aa[' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'bb[' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'exponent' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'amplitude' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'condition' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Reduce' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Solve' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'SolveAlways' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
1
```

```
$ /usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

```
$ /usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
0
```

The one `Claude` hit, if present, is the `~/.claude/projects` path in the do-not-read line:

```
$ /usr/bin/grep -n -i -F -- 'claude' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
11:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
```

```
$ grep -c -F 'resume_s15_reading_repair_4.md' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md; grep -c -F 'Do not read the orchestrator' _scratch/s15/reading_build/resume_s15_reading_repair_4c.md
1
1
```

Facts the note states:

```
$ cat _scratch/s11c/s15/repair3_90441_native_mixed_momentum_domain_01/outcome.json
{
  "exitCode": 1,
  "wallSeconds": 2790.3142783610383,
  "unit": "s11c-guard-c3c98b10f9d3",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 1,
    "guardReason": null,
    "wallSeconds": 2790.1777805419406,
    "childPid": 407015
  }
}
```

```
$ grep -c -F 'No more memory available' _scratch/s11c/s15/repair3_90441_native_mixed_momentum_domain_01/stdout
1
```

```
$ python3 -c 'import json;print(max(int(json.loads(l).get("memory.peak",0)) for l in open("_scratch/s11c/s15/repair3_90441_native_mixed_momentum_domain_01/resource-samples.jsonl") if l.strip()))'
52120858624
```

```
$ grep -E '^(MemTotal|MemAvailable):' /proc/meminfo
MemTotal:       65779296 kB
MemAvailable:   59313296 kB
```

```
$ grep -n -E '^AVAILABLE_MIN =|^POOL_MEMORY_MAX =' scripts/s11c_guarded_run.py
76:AVAILABLE_MIN = 0
86:POOL_MEMORY_MAX = physical_memory()
```

```
$ grep -n -F -- '--memory-gib 50' _scratch/s15/reading_build/brief_s15_reading_repair_3.md
228:`python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/<a new directory name per run> -- <command>`.
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; systemctl --user list-units --type=service --state=active --no-legend | grep -c s11c-guard; pgrep -c codex; pgrep -c WolframKernel
{}

0
0
0
```

```
$ cat _scratch/s15/reading_build/repair4.DONE
exit 0
```
