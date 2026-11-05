# Measurements: S15 reading build repair-4m note (user: run the expansion route longer), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes, process and
unit state. Lines cut to a byte width can end inside a character; the record drops such partial characters. Every hit is
read by the orchestrator. The run and work directories named in these commands are local scratch and disposable; every
output this note relies on is quoted here in full. The note carries a user decision on how the one-object trial of
`resume_s15_reading_repair_4l.md` (reviewed by two legs and gated at `37aeedeb`, updated at `e5e5a220`) continues: the
same object, with the same construction in each engine, runs until each attempt ends on its own, in place of the 1-hour
trial stop. It changes no supplied equation, no object, no construction, no contract and nothing a record may claim (a
trial outcome is still not an item-4 label), so it is an operational note, gated like `4i`, `4j` and `4k`.

## The user's decision (2026-11-05)

AskUserQuestion, ~09:45Z, after the builder's trial report. Question: "The expansion route didn't finish the test object
in an hour on either engine. How should item 4 go from here?" Options offered:
- "Finish with labels (Recommended)": "Back to your 11-02 plan: final runs for all 48 objects. Any object that doesn't
  finish in 6 hours, or runs out of memory, is labelled 'not obtained' (never 'no condition') and goes on the deferred
  list for a bigger machine. Probably 2–4 days of machine time, mostly the six two-drain objects. That figure comes from
  older runs, and SymPy's pace on the finals is unknown."
- "Finish, label the six now": "Same, but the six two-drain region objects take their 'not obtained' label from the runs
  already made, with no new 6-hour attempt each. Saves up to about 3 days. A short builder note gets two review legs first."
- "Rethink item 4": "Change what item 4 asks for on the six two-drain objects, for example a numeric check first, with
  exact work only where it's needed. This changes the reading decision list, so it gets two review legs and then a new
  build. It takes longer but could give answers instead of labels. I'd need your direction on what to ask for."
- "Run the expansion longer": "Run the one test object alone on the whole machine with no time limit, to see whether the
  expansion route finishes at all. Wolfram had already used 30 GiB after one hour, so it may run out of memory."

The user selected **"Run the expansion longer"**. Before the question the orchestrator told the user that the SymPy rate
measurement keeps running unless they say otherwise; they did not. So each trial attempt is the only trial job, beside
that measurement, at the largest cap the guard admits there (R5 quotes the caps the guard admitted in the 1-hour trial).

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
f8f43c1d596417d99aa36190f666a6ac414d9797961a0d6616b83b1e1174a81e  research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
2702 research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gegenbauer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Gegenbauer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
1:Note from the orchestrator. I have read your trial report (`series_trial_report.md`). The user has made a further decision about the trial.
5:Neither engine returned the conditions of `region/component_1/momentum_labels_{drain1, drain2}` within 1 hour. The user decided to find out whether the construction that uses the supplied series ret
19:- Attempt `region/component_1/momentum_labels_{drain1, drain2}` in Wolfram and then in SymPy, with the construction of your last trial attempt in each engine (`repair4_series_wl_c1_momentum_04` and
count: 3
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[|hh[12])' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[|hh[12])' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
22:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
19:- Attempt `region/component_1/momentum_labels_{drain1, drain2}` in Wolfram and then in SymPy, with the construction of your last trial attempt in each engine (`repair4_series_wl_c1_momentum_04` and
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'GiB' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'GiB' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md)"
19:- Attempt `region/component_1/momentum_labels_{drain1, drain2}` in Wolfram and then in SymPy, with the construction of your last trial attempt in each engine (`repair4_series_wl_c1_momentum_04` and
20:- When both engines' attempts have ended, stop and write a short report. For each attempt give the run, its start and end, its `--memory-gib`, whether it returned its conditions, ran out of memory 
count: 2
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'series_trial_report.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4g.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4l.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
2
```

```
$ /usr/bin/grep -c -F -- 'except the 1-hour trial stop' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'No attempt is stopped for time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'repair4_size_py_c1_rate_01' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 's11c-guard-fd3dd4540a3c' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'repair4_series_wl_c1_momentum_04' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'repair4_series_py_c1_momentum_02' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'one trial job at a time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'report every attempt' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not attempt any other object' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'do not start the final runs' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'No production source changes' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'not an item-4 label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
1
```

```
$ ls research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
```

```
$ git log -1 --format='%h %s' -- research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
e5e5a220 Update the S15 repair-4l note: the SymPy momentum measurement attempt has ended
```

```
$ /usr/bin/grep -n -F -e 'What must be true' -e 'the 1-hour trial stop of' -e 'Each attempt stops 1 hour after it started' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-240
27:## What must be true
29:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_reading_repair_4h.md` and `resume_s15_reading_
36:- Write the construction in each engine, in your work directory, and attempt `region/component_1/momentum_labels_{drain1, drain2}` with it, Wolfram first. Run through the guard, one job at a time, with the largest `--memory-gib` the guar
```

```
$ /usr/bin/grep -n -F -e '1 hour after it started' -e 'No trial outcome is printed as an item-4 label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-240
9:The trial measures how fast your construction is. The user decided that in this trial, an attempt is stopped 1 hour after it started and is reported as not returned within 1 hour. The Wolfram trial of `region/component_1/momentum_labels_{
16:- No trial outcome is printed as an item-4 label. A trial attempt ends when it returns its conditions, when it runs out of memory, or 1 hour after it started.
23:- Run through the guard, one job at a time, with `--memory-gib 55`, and poll each job until it ends. Each attempt stops 1 hour after it started. If an attempt ends on an error in your code, you may fix it and run that object again; repor
```

```
$ sed -n '17p' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-240
- every item-4 object that carries a label carries one of the two kinds, the 6-hour stop or running out of memory, and the run record says which.
```

## R5. Facts the note states

The trial report, its end and the two timed attempts it names.

```
$ cat _scratch/s15/reading_build/repair4l.DONE; sha256sum _scratch/s15/reading_build_work/repair3/series_trial/series_trial_report.md
exit 0
b8b1da310c7e2d4840445b04ea8c060581ada64073f6771492d70a214910c515  _scratch/s15/reading_build_work/repair3/series_trial/series_trial_report.md
```

```
$ /usr/bin/grep -n -F -e 'not_returned_within_one_hour' -e 'No complete condition classification returned' -e 'Production sources were not changed' -e 'Stop here' _scratch/s15/reading_build_work/repair3/series_trial/series_trial_report.md | cut -c1-240
3:Only `region/component_1/momentum_labels_{drain1, drain2}` was attempted, independently in Wolfram and SymPy, Wolfram first. Production sources were not changed. The original SymPy rate measurement was neither stopped nor restarted. No fi
14:| wl_c1_momentum_04 | 07:24:48 | 08:24:55 | not_returned_within_one_hour | 3620.4076844048686 | 32393408512 |
16:| py_c1_momentum_02 | 08:39:31 | 09:39:35 | not_returned_within_one_hour | 4005.1941500520334 | 660975616 |
20:No complete condition classification returned in this trial. Consequently there is no returned classification whose completeness or discarded-branch inventory can be certified. Partial conditions and Q values are provisional. The Wolfram
24:Open result: this trial has not established a complete, lower-cost item-4 construction. Stop here for the orchestrator’s reply.
```

```
$ for r in repair4_series_wl_c1_momentum_04 repair4_series_py_c1_momentum_02; do python3 -c "import json;j=json.load(open('_scratch/s15/reading_build_work/repair3/series_trial/$r'+'_result.json'));print('$r',j['result'],j['cap_gib'],j['source_sha256'],j['execution']['source'].rsplit('/',1)[1],j['unit_after_stop'],j['reservation_after_stop'])"; done
repair4_series_wl_c1_momentum_04 not_returned_within_one_hour 53 c2b52ca9b9d1b74789ad84b508938af2b0f9c1a1e088e409df6e9f26c470cba8 wolfram_series_trial.wl inactive None
repair4_series_py_c1_momentum_02 not_returned_within_one_hour 54 34db16b8f4a3bd980c8e9e279899dae34a3afb6f0bd42e5c7c8388d107d1bf60 sympy_series_trial.py inactive None
```

```
$ for r in repair4_series_wl_c1_momentum_04 repair4_series_py_c1_momentum_02; do printf '%s ' $r; cat _scratch/s11c/s15/$r/outcome.json | tr -d '\n '; echo; done
repair4_series_wl_c1_momentum_04 {"exitCode":1,"wallSeconds":3620.4076844048686,"unit":"s11c-guard-c487f523cee6","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
repair4_series_py_c1_momentum_02 {"exitCode":1,"wallSeconds":4005.1941500520334,"unit":"s11c-guard-66ff6e5e34a0","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
```

```
$ sha256sum _scratch/s15/reading_build_work/repair3/series_trial/wolfram_series_trial.wl _scratch/s15/reading_build_work/repair3/series_trial/sympy_series_trial.py
c2b52ca9b9d1b74789ad84b508938af2b0f9c1a1e088e409df6e9f26c470cba8  _scratch/s15/reading_build_work/repair3/series_trial/wolfram_series_trial.wl
34db16b8f4a3bd980c8e9e279899dae34a3afb6f0bd42e5c7c8388d107d1bf60  _scratch/s15/reading_build_work/repair3/series_trial/sympy_series_trial.py
```

The SymPy measurement attempt and the guard registry.

```
$ systemctl --user is-active s11c-guard-fd3dd4540a3c.service; python3 -c "import json;r=json.load(open('/home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json'));[print(k,v['memoryMax'],v['owner']['pid'],v['logDirectory'].rsplit('/',1)[1]) for k,v in r.items()]"
active
s11c-guard-fd3dd4540a3c.service 2147483648 237635 repair4_size_py_c1_rate_01
```

```
$ test -e _scratch/s11c/s15/repair4_size_py_c1_rate_01/outcome.json && echo 'outcome.json present' || echo 'outcome.json absent'
outcome.json absent
```

```
$ echo "WolframKernel processes: $(pgrep -c WolframKernel)"; echo "processes named codex: $(pgrep -c -x codex)"
WolframKernel processes: 0
processes named codex: 0
```

Production sources, final runs.

```
$ sha256sum research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f  research/pde_ledger_v3/scripts/S15_separation_departure.py
4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1  research/pde_ledger_v3/mathematica/S15_separation_departure.wl
b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ git status --short -- research/pde_ledger_v3
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4m_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ echo "final run dirs after final_20: $(ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null | wc -l)"
final run dirs after final_20: 0
```
