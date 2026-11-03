# Measurements: S15 reading build repair-4i note (user: 1-hour trial limit), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit
is read by the orchestrator. The run and work directories named in these commands are local scratch and disposable;
every output this note relies on is quoted here in full. The note reports that the orchestrator stopped the guarded trial
run `repair4_scalar_wl_c1_momentum_03` on 2026-11-03 at the user's decision, and carries that decision: in the trial of
a cheaper item-4 construction asked for by `resume_s15_reading_repair_4h.md` (gated at `d371f16c`), an attempt stops 1 hour
after it started; a trial stop is not an item-4 label; the final runs keep the rule of `resume_s15_reading_repair_4g.md`
(`1a703b67`). The note names no value, exponent or condition.

## The user's words and decision (2026-11-03)

AskUserQuestion, ~21:00Z. Question: "The trial's first object has taken ~3 hours across 3 tries, not minutes. How should
the trial continue?" Answer chosen: **"1-hour trial limit (Recommended)"**: "Stop object 1 now (1.6 h this try, memory
climbing). Run the other 11 trial objects with a 1-hour limit each, then report which returned and how fast. A trial stop is
not a label; the final runs keep the 6-hour rule. At most ~11 more hours before you see the numbers."

The other options offered: "Keep going as is", "Stop trial, finish old way", "Stop, consider numeric-first".

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
dd3879f7a7a1089abf629bfb1f327ca1a9935250ed6e2e6be6853bea900026e2  research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
3024 research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
5:On 2026-11-03 I stopped `repair4_scalar_wl_c1_momentum_03` through its unit (`systemctl --user stop s11c-guard-a21d67f54053`). The journal records the stop request at 21:05:48Z and the unit stopped 
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
26:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
5:On 2026-11-03 I stopped `repair4_scalar_wl_c1_momentum_03` through its unit (`systemctl --user stop s11c-guard-a21d67f54053`). The journal records the stop request at 21:05:48Z and the unit stopped 
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md)"
9:The trial measures how fast your construction is. The user decided that in this trial, an attempt is stopped 1 hour after it started and is reported as not returned within 1 hour. The Wolfram trial 
21:  - in Wolfram first, `region/component_1/rate_labels_{drain1, drain2}`, `region/component_2/momentum_labels_{drain1, drain2}`, `region/component_2/rate_labels_{drain1, drain2}`, `region/component_
count: 2
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4f.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4g.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4h.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
2
```

```
$ /usr/bin/grep -c -F -- 'method_trial_external_stop_03.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
2
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'one job at a time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- '1 hour' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
5
```

```
$ /usr/bin/grep -c -F -- '6 hours' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not start the final runs' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'poll each job until it ends' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
1
```

```
$ ls _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
```

```
$ /usr/bin/grep -n -F 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:18:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:7:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:9:Everything under "What must be true" in `resume_s15_reading_repair_4.md` stands. For an object whose attempt ran out of memory, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md` and in `resume_s15_reading_repair_4d.md` stands. For an object whose attempt was stopped at 6 hours, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md` and `resume_s15_reading_repair_4f.md` stands. This decision replaces the bullet in `resume_s15_reading_repair_4f.md` about an attempt that goes silent. Also:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md` and `resume_s15_reading_repair_4g.md` stands, with these changes and additions:
```

```
$ /usr/bin/grep -n -F -e 'only when it runs out of memory' -e '6 hours after it started' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-220
9:`resume_s15_reading_repair_4f.md` said that an item-4 attempt that goes silent for a long stretch is labelled, and it did not say how long. That sentence was mine; the user's decision was the 6-hour stop. The user has 
```

```
$ /usr/bin/grep -n -F -e 'region/component_' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-120
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_20` through its unit (`systemctl --user stop s11c-guard-1af66ae63
24:- Try the construction on the six objects `region/component_1/momentum_labels_{drain1, drain2}`, `region/component_1/
```

## R5. Facts the note states

S0, captured by the orchestrator at 2026-11-03T21:05:44Z, before the stop (commands and literal output):

````
S0 snapshot 2026-11-03T21:05:44Z
--- registry
{
  "s11c-guard-a21d67f54053.service": {
    "pool": "s15",
    "poolMemoryMax": 66571993088,
    "memoryMax": 59055800320,
    "tasksMax": 32,
    "cpu": 7,
    "owner": {
      "pid": 3544174,
      "startTicks": "100390329",
      "bootId": "c3b7b2c7-a5fc-4b45-adf0-6cbde3d609bb"
    },
    "logDirectory": "/var/projects/toy_physics/_scratch/s11c/s15/repair4_scalar_wl_c1_momentum_03"
  }
}

--- unit
ControlGroup=/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-a21d67f54053.service
ActiveState=active
SubState=running
ActiveEnterTimestamp=Tue 2026-11-03 12:22:34 MST
--- memory.stat
anon 21336997888
file 1725116416
inactive_file 1725112320
memory.peak 23157850112
memory.max 59055800320
oom 0
oom_kill 0
--- procs
3544181 3544180 99.8 20922344 01:43:10 /var/projects/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel -script /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3.wl
3505258 3505257  0.0 230916   03:04:52 codex exec resume -m gpt-6-astra -c model_reasoning_effort=high -c sandbox_mode="danger-full-access" 01a18cb6-20cc-7a40-b214-5419860e4fb6 -
3505716 3505258  0.0 26220    03:04:44 /var/projects/_home/.codex/packages/standalone/releases/0.162.1-x86_64-unknown-linux-musl/bin/codex-code-mode-host
3544180 3544178  0.0 17988    01:43:10 python3 /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/run_trial.py --source /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3
3544178    1646  0.0 15772    01:43:10 /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/toy_physics/_scratch/s11c/s15/repair4_scalar_wl_c1_momentum_03/invocation.json
3544174 3505258  0.0 15452    01:43:10 python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 55 --log-directory /var/projects/toy_physics/_scratch/s15/guard/repair4_scalar_wl_c1_momentum_03 -- python3 /var/projects/toy
3563798 3505258  0.0 12900    01:35:55 python3 _scratch/s15/reading_build_work/repair3/method_trial/sequence_trials_v3.py
3544177 3544174  0.0  7636    01:43:10 systemd-run --user --quiet --wait --pipe --collect --unit=s11c-guard-a21d67f54053 --property=MemoryMax=59055800320 --property=MemorySwapMax=0 --property=TasksMax=32 --property=RuntimeMaxSec=infinity --property=R
3505257  175046  0.0  3780    03:04:52 /bin/bash -c source /home/trevnorris/.claude/shell-snapshots/snapshot-bash-1793555799511-nktnl5.sh 2>/dev/null || true && shopt -u extglob 2>/dev/null || true && { \builtin unalias -- 'unsetenv'; \builtin unset 
--- run dir
total 32800
drwxrwxr-x   2 trevnorris trevnorris     4096 Nov  3 12:22 .
drwxrwxr-x 247 trevnorris trevnorris    20480 Nov  3 12:22 ..
-rw-rw-r--   1 trevnorris trevnorris      164 Nov  3 12:22 duration-validation.json
-rw-rw-r--   1 trevnorris trevnorris      597 Nov  3 12:22 effective-limits.json
-rw-rw-r--   1 trevnorris trevnorris     1375 Nov  3 12:22 invocation.json
-rw-rw-r--   1 trevnorris trevnorris      679 Nov  3 12:22 limit-validation.json
-rw-rw-r--   1 trevnorris trevnorris  1251274 Nov  3 14:05 resource-samples.jsonl
-rw-rw-r--   1 trevnorris trevnorris        0 Nov  3 12:22 stderr
-rw-rw-r--   1 trevnorris trevnorris 32281821 Nov  3 14:05 stdout
stdout WORKING_ON 206315
WORKING_ON: reading/region/component_1/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / scalar quotient 1 / eliminate series pivot 26
--- guard log dir
total 32800
drwxrwxr-x   2 trevnorris trevnorris     4096 Nov  3 12:22 .
drwxrwxr-x 247 trevnorris trevnorris    20480 Nov  3 12:22 ..
-rw-rw-r--   1 trevnorris trevnorris      164 Nov  3 12:22 duration-validation.json
-rw-rw-r--   1 trevnorris trevnorris      597 Nov  3 12:22 effective-limits.json
-rw-rw-r--   1 trevnorris trevnorris     1375 Nov  3 12:22 invocation.json
-rw-rw-r--   1 trevnorris trevnorris      679 Nov  3 12:22 limit-validation.json
-rw-rw-r--   1 trevnorris trevnorris  1251274 Nov  3 14:05 resource-samples.jsonl
-rw-rw-r--   1 trevnorris trevnorris        0 Nov  3 12:22 stderr
-rw-rw-r--   1 trevnorris trevnorris 32281821 Nov  3 14:05 stdout
--- work dir
total 1647960
drwxrwxr-x 2 trevnorris trevnorris       4096 Nov  3 12:22 .
drwxrwxr-x 5 trevnorris trevnorris       4096 Nov  3 12:29 ..
-rw-rw-r-- 1 trevnorris trevnorris        606 Nov  3 12:22 execution.json
-rw-rw-r-- 1 trevnorris trevnorris   32281821 Nov  3 14:05 kernel.stdout
-rw-rw-r-- 1 trevnorris trevnorris 1655205989 Nov  3 14:05 transcript.out
{
  "run": "repair4_scalar_wl_c1_momentum_03",
  "object": "region/component_1/momentum_labels_{drain1, drain2}",
  "source": "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3.wl",
  "source_sha256": "6eaf145ecf6bb609820318f58616412678f8efea49834efe182b5826b7dcb20d",
  "started": "2026-11-03T19:22:34.176956+00:00",
  "command": [
    "math",
    "-script",
    "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3.wl"
  ],
  "pid": 3544181,
  "attempt_started": "2026-11-03T19:22:48.778857+00:00"
}

--- builder log
311760 _scratch/s15/reading_build/repair4h.log
2026-11-03 13:33:36.498861113 -0700
--- trial files sha256
a866f015a45cee799c69ea9410c9a5c1a9a812ff24a150115dec85caa04bf048  _scratch/s15/reading_build_work/repair3/method_trial/add_positive_normalization.py
01e591577831a54eb01fe5bb65acc005238be03849a51653df7a8ef844efb700  _scratch/s15/reading_build_work/repair3/method_trial/add_sympy_positive_normalization.py
9ed4eb258af3d659229bc20e64e6c1e56941ce95b956b5826831ad66792873ca  _scratch/s15/reading_build_work/repair3/method_trial/build_sympy_scalar.py
0a77a0fc4db6b375f1a0cf4ac31218e932a1df90e13f5540f2e7b66d565509b5  _scratch/s15/reading_build_work/repair3/method_trial/build_wolfram_scalar.py
6c5977030c34bc8d728090baa519eaa92644b8cc978d4983346781f0c8dce69e  _scratch/s15/reading_build_work/repair3/method_trial/run_py_trial.py
55308071eb9366d6efc71e01418b01dc582ed96b3fc023256ed4a40fbd35bed2  _scratch/s15/reading_build_work/repair3/method_trial/run_trial.py
1fff7743e439f7c88e4c7d2786aad883cac6c40888e2b5fb57b6c4f601af2496  _scratch/s15/reading_build_work/repair3/method_trial/sequence_trials.py
9772f7a4a69a6f4c95d05de3f2ac901a83905fcc3fd5007009973c8d4a20268f  _scratch/s15/reading_build_work/repair3/method_trial/sequence_trials_v3.py
7504649916817f386eb8f5d0987fdc0674ed8e22846ca76f173de6453210d904  _scratch/s15/reading_build_work/repair3/method_trial/sympy_scalar_trial.py
d16b79d7983bfa155099c4249f511741c885608e8371d94321d0aa4a4b0379d2  _scratch/s15/reading_build_work/repair3/method_trial/sympy_scalar_trial_v2.py
db3c581af28e770fe65ac29bd216f3985f502721479285cdc7bfff8d951dad17  _scratch/s15/reading_build_work/repair3/method_trial/watch_trial.py
fc6118660e4119a42de14e61e5559e05a9df78e0672a820aebe4662f708d3300  _scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v2.wl
6eaf145ecf6bb609820318f58616412678f8efea49834efe182b5826b7dcb20d  _scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3.wl
f035569c75bfa95ab7c45a1ff29234c6fa0384f267375ed6589e28f9456c1a1e  _scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial.wl
--- repo
 M .claude/skills/prep-for-compact/SKILL.md
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1  research/pde_ledger_v3/mathematica/S15_separation_departure.wl
48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f  research/pde_ledger_v3/scripts/S15_separation_departure.py
b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
````

```
$ journalctl --user -u s11c-guard-a21d67f54053.service --no-pager -o short-iso-precise --since '2026-11-03 21:05:30 UTC' --until '2026-11-03 21:06:30 UTC' --utc | cut -c1-200
2026-11-03T21:05:48.899995+00:00 pop-os systemd[1646]: Stopping s11c-guard-a21d67f54053.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-03T21:05:54.058062+00:00 pop-os systemd[1646]: s11c-guard-a21d67f54053.service: State 'final-sigterm' timed out. Killing.
2026-11-03T21:05:54.058147+00:00 pop-os systemd[1646]: s11c-guard-a21d67f54053.service: Killing process 3544181 (WolframKernel) with signal SIGKILL.
2026-11-03T21:05:54.058249+00:00 pop-os systemd[1646]: s11c-guard-a21d67f54053.service: Failed with result 'timeout'.
2026-11-03T21:05:54.058278+00:00 pop-os systemd[1646]: s11c-guard-a21d67f54053.service: Unit process 3544181 (WolframKernel) remains running after unit stopped.
2026-11-03T21:05:54.058535+00:00 pop-os systemd[1646]: Stopped s11c-guard-a21d67f54053.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-03T21:05:54.058689+00:00 pop-os systemd[1646]: s11c-guard-a21d67f54053.service: Consumed 1h 43min 17.859s CPU time, 21.5G memory peak, 0B memory swap peak.
```

```
$ cat _scratch/s11c/s15/repair4_scalar_wl_c1_momentum_03/outcome.json; echo; cat _scratch/s15/guard/repair4_scalar_wl_c1_momentum_03/outcome.json; echo; /usr/bin/grep -a -c -F 'No more memory available' _scratch/s11c/s15/repair4_scalar_wl_c1_momentum_03/stdout
{
  "exitCode": 1,
  "wallSeconds": 6199.971640525968,
  "unit": "s11c-guard-a21d67f54053",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}

{
  "exitCode": 1,
  "wallSeconds": 6199.971640525968,
  "unit": "s11c-guard-a21d67f54053",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}

0
```

```
$ cat _scratch/s15/reading_build_work/repair3/method_trial/wl_c1_momentum_03/execution.json
{
  "run": "repair4_scalar_wl_c1_momentum_03",
  "object": "region/component_1/momentum_labels_{drain1, drain2}",
  "source": "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3.wl",
  "source_sha256": "6eaf145ecf6bb609820318f58616412678f8efea49834efe182b5826b7dcb20d",
  "started": "2026-11-03T19:22:34.176956+00:00",
  "command": [
    "math",
    "-script",
    "/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/method_trial/wolfram_scalar_trial_v3.wl"
  ],
  "pid": 3544181,
  "attempt_started": "2026-11-03T19:22:48.778857+00:00"
}
```

```
$ for r in 01 02 03; do echo "$r $(stat -c '%y' _scratch/s11c/s15/repair4_scalar_wl_c1_momentum_$r/invocation.json) $(cat _scratch/s11c/s15/repair4_scalar_wl_c1_momentum_$r/outcome.json | tr -d '\n ')"; done
01 2026-11-03 11:03:55.547801251 -0700 {"exitCode":1,"wallSeconds":78.0212976608891,"unit":"s11c-guard-ed6ed80e6ab7","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
02 2026-11-03 11:05:27.107316975 -0700 {"exitCode":1,"wallSeconds":4574.71407812892,"unit":"s11c-guard-a72a6f8c5032","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
03 2026-11-03 12:22:34.099253008 -0700 {"exitCode":1,"wallSeconds":6199.971640525968,"unit":"s11c-guard-a21d67f54053","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
```

```
$ ls -l _scratch/s15/reading_build_work/repair3/method_trial/method_trial_external_stop_03.md | awk '{print $5, $9}'; sha256sum _scratch/s15/reading_build_work/repair3/method_trial/method_trial_external_stop_03.md
2883 _scratch/s15/reading_build_work/repair3/method_trial/method_trial_external_stop_03.md
daf02418db194445abb047ebac07c1e837748e9524976d73cdb7bb9e2afa1924  _scratch/s15/reading_build_work/repair3/method_trial/method_trial_external_stop_03.md
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c -f 'codex exec'; ls _scratch/s15/reading_build/repair4h.DONE 2>&1; cat _scratch/s15/reading_build/repair4h.DONE 2>&1
{}

0
0
_scratch/s15/reading_build/repair4h.DONE
exit 0
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md; echo "run record diff vs 90441b89 exit $?"; git status --short -- research/pde_ledger_v3
run record diff vs 90441b89 exit 0
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4i_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null; echo "final run dirs after final_20: $(ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null | wc -l)"
final run dirs after final_20: 0
```
