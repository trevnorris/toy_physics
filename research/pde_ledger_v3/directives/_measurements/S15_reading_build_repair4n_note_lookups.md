# Measurements: S15 reading build repair-4n note (user: stop the SymPy attempts, label the six two-drain objects, finish), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes, process and
unit state. Lines cut to a byte width can end inside a character; the record drops such partial characters. Every hit is
read by the orchestrator. The run and work directories named in these commands are local scratch and disposable; every
output this note relies on is quoted here in full. The note carries user decisions that change what the S15 run record may
claim for six item-4 objects (a third label kind) and how the final runs are made, so it is physics-bearing and gets two
review legs before commit.

## The user's decisions (2026-11-06)

Asked during the 4m long run, about 11.5 h into the SymPy attempt. The user's words before the question: "What I'm getting
at is that if this script will also run out of memory, just a lot slower. Seems like we need to rethink this. Perhaps this
goes onto the list of runs we need more memory for."

AskUserQuestion: "What should happen to the SymPy attempt that's running now (11.5 h, 1.8 GB, silent since 11:23 am)?"
The user selected **"Stop it, finish S15 (Recommended)"**: "Stop it and record it as not returned. The six two-drain objects
get 'not obtained' labels and go on the bigger-machine list; the other 42 objects get their final runs. A short builder note
gets two review legs first."

AskUserQuestion: "The SymPy rate measurement is the only SymPy job still running (about 40 h, one of the six two-drain
objects that will now be labelled). Stop it too?" The user selected **"Stop it too (Recommended)"**.

Parallel final runs follow the user's standing rule for this machine (2026-10-17, their words): "These jobs take hours to
run, so running multiple in parallel will save us a lot of time."; and (2026-10-18) "We set limits so we can run things in
parallel. If we're not running it in parallel then it can use the whole box."

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
8f144178b5bfe9bb562b2cd4209b329ff71dc1e8e318ccddc582c2219e2755e7  research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
5057 research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
5:On 2026-11-06, on the user's decision:
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
6:- `repair4_series_unlimited_py_c1_momentum_01` (unit `s11c-guard-350cfe7325e2`). At about 05:40Z I ended your session. That also ended `launch_trial.py` and the guard launcher of this run, while its
count: 1
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
6:- `repair4_series_unlimited_py_c1_momentum_01` (unit `s11c-guard-350cfe7325e2`). At about 05:40Z I ended your session. That also ended `launch_trial.py` and the guard launcher of this run, while its
count: 1
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gegenbauer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Gegenbauer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
6:- `repair4_series_unlimited_py_c1_momentum_01` (unit `s11c-guard-350cfe7325e2`). At about 05:40Z I ended your session. That also ended `launch_trial.py` and the guard launcher of this run, while its
count: 1
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[|hh[12])' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[|hh[12])' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
32:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Singular' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Singular' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Maple' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Maple' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md)"
count: 0
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4f.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4g.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4h.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4i.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4l.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4m.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
2
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
2
```

```
$ /usr/bin/grep -c -F -- 'third kind of label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'not attempted in that run' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'listed as a trial' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- '6 hours after it started' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'has not run out of memory' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'when no other guarded job is running' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'the kind of each label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'For these six objects it replaces the rules' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'For every other object those rules stand' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'whether another guarded job was running during it' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'none of them gives any of the six' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'machine with more memory' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'final guarded run of each engine' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'at the same time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'Stop there' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```
$ ls research/pde_ledger_v3/directives/resume_s15_reading_repair_4{f,g,h,i,l,m}.md _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md _scratch/s15/reading_build/brief_s15_reading_repair_3.md
_scratch/s15/reading_build/brief_s15_reading_repair_3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md
```

```
$ for f in 4g 4h 4i 4l 4m; do git log -1 --format="$f %h %s" -- research/pde_ledger_v3/directives/resume_s15_reading_repair_$f.md; done
4g 1a703b67 Gate the S15 reading build repair-4g note (user: 6 hours only)
4h d371f16c Gate the S15 reading build repair-4h note (user: try a faster method first)
4i 4d2c29fc Gate the S15 reading build repair-4i note (user: 1-hour trial limit)
4l e5e5a220 Update the S15 repair-4l note: the SymPy momentum measurement attempt has ended
4m 8c81d55a Gate the S15 reading build repair-4m note (user: run the expansion route longer)
```

```
$ /usr/bin/grep -n -F -e 'one of the two kinds' -e 'Run the jobs one at a time' -e 'Finish as the brief says' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-240
17:- every item-4 object that carries a label carries one of the two kinds, the 6-hour stop or running out of memory, and the run record says which.
22:- The cap stays `--memory-gib 55`. Run the jobs one at a time.
23:- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.
```

```
$ /usr/bin/grep -n -F -e 'its label says what was attempted' -e 'the run record lists, for each engine' research/pde_ledger_v3/../../_scratch/s15/reading_build/resume_s15_reading_repair_4d.md 2>/dev/null | cut -c1-240; /usr/bin/grep -n -F -e 'its label says what was attempted' -e 'the run record lists, for each engine' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-240
10:- its label says what was attempted, that the attempt ran out of memory, and where: the last step the attempt reported and the peak memory of the run;
13:- the run record lists, for each engine, every object that carries such a label.
10:- its label says what was attempted, that the attempt ran out of memory, and where: the last step the attempt reported and the peak memory of the run;
13:- the run record lists, for each engine, every object that carries such a label.
```

```
$ /usr/bin/grep -n -F -e 'which of the two it is' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-200; /usr/bin/grep -n -F -e 'no label is printed for any other reason' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; /usr/bin/grep -n -F -e '## What must be true' -e '## Going on' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
20:- the run record lists, for each engine, every object that carries such a label or the label for running out of memory, and which of the two it is.
19:- An item-4 object either prints its complete conditions or carries one of the two labels; no label is printed for any other reason. The assessment reports that each engine can print a not-obtained
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:22:## Going on
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:19:## Going on
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md:13:## What must be true
```

```
$ /usr/bin/grep -n -F -e 'No production source changes' -e 'do not start the final runs' -e 'No attempt is stopped for time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4m.md | cut -c1-200
14:- No attempt is stopped for time or for a period without output. An attempt ends when it returns its conditions, when it runs out of memory at its cap, or on an error in your code.
15:- No production source changes.
20:- When both engines' attempts have ended, stop and write a short report. For each attempt give the run, its start and end, its `--memory-gib`, whether it returned its conditions, ran out of memory 
```

## R5. Facts the note states

The six objects: each name as the note gives it, counted in the note and in `resume_s15_reading_repair_4h.md`.

```
$ printf 'note %s  4h %s\n' $(/usr/bin/grep -c -F -- 'region/component_1/momentum_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md) $(/usr/bin/grep -c -F -- 'region/component_1/momentum_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)
note 1  4h 1
```

```
$ printf 'note %s  4h %s\n' $(/usr/bin/grep -c -F -- 'region/component_1/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md) $(/usr/bin/grep -c -F -- 'region/component_1/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)
note 1  4h 1
```

```
$ printf 'note %s  4h %s\n' $(/usr/bin/grep -c -F -- 'region/component_2/momentum_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md) $(/usr/bin/grep -c -F -- 'region/component_2/momentum_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)
note 1  4h 1
```

```
$ printf 'note %s  4h %s\n' $(/usr/bin/grep -c -F -- 'region/component_2/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md) $(/usr/bin/grep -c -F -- 'region/component_2/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)
note 1  4h 2
```

```
$ printf 'note %s  4h %s\n' $(/usr/bin/grep -c -F -- 'region/component_3/momentum_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md) $(/usr/bin/grep -c -F -- 'region/component_3/momentum_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)
note 1  4h 2
```

```
$ printf 'note %s  4h %s\n' $(/usr/bin/grep -c -F -- 'region/component_3/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md) $(/usr/bin/grep -c -F -- 'region/component_3/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)
note 1  4h 2
```

```
$ /usr/bin/grep -n -F 'None of the six region objects that carry both drains has returned its conditions' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200
9:None of the six region objects that carry both drains has returned its conditions in any run of this repair: two ran out of memory within minutes, three were stopped at 6 hours or more, and I stoppe
```

Runs on these objects after `resume_s15_reading_repair_4h.md`: the method trial, the scope measurements, the series trial and the long run (process keys only).

```
$ /usr/bin/grep -n -F 'did not return complete conditions within one hour for any of the twelve engine/object pairs' _scratch/s15/reading_build_work/repair3/method_trial/method_trial_report.md | cut -c1-200
3:**The construction did not return complete conditions within one hour for any of the twelve engine/object pairs.** None of these trial runs exhausted memory. The five newly timed Wolfram attempts an
```

```
$ for d in _scratch/s11c/s15/repair4_size_*; do printf '%s ' $(basename $d); test -e $d/outcome.json && python3 -c "import json,sys;o=json.load(open(sys.argv[1]));print('exitCode',o.get('exitCode'),'wallSeconds',round(o.get('wallSeconds',0),1))" $d/outcome.json || echo 'no outcome.json'; done
repair4_size_py_c1_momentum_01 exitCode 1 wallSeconds 60837.2
repair4_size_py_c1_rate_01 exitCode 1 wallSeconds 143174.4
repair4_size_wl_c1_momentum_01 exitCode 0 wallSeconds 17420.5
repair4_size_wl_c1_momentum_02 no outcome.json
repair4_size_wl_c1_rate_01 exitCode 0 wallSeconds 7436.3
```

```
$ /usr/bin/grep -n -F -e 'No complete condition classification returned' -e 'not_returned_within_one_hour' _scratch/s15/reading_build_work/repair3/series_trial/series_trial_report.md | cut -c1-200
14:| wl_c1_momentum_04 | 07:24:48 | 08:24:55 | not_returned_within_one_hour | 3620.4076844048686 | 32393408512 |
16:| py_c1_momentum_02 | 08:39:31 | 09:39:35 | not_returned_within_one_hour | 4005.1941500520334 | 660975616 |
20:No complete condition classification returned in this trial. Consequently there is no returned classification whose completeness or discarded-branch inventory can be certified. Partial conditions a
```

```
$ python3 -c "import json;j=json.load(open('_scratch/s15/reading_build_work/repair3/series_unlimited/repair4_series_unlimited_wl_c1_momentum_01_result.json'));print({k:j.get(k) for k in ('result','stop_kind','cap_gib','started','ended','peak_memory_bytes')})"
{'result': 'ran_out_of_memory', 'stop_kind': None, 'cap_gib': 54, 'started': '2026-11-05T14:48:45.126183+00:00', 'ended': '2026-11-05T18:02:10.265255+00:00', 'peak_memory_bytes': None}
```

The SymPy long attempt and its stop.

```
$ ls _scratch/s11c/s15/repair4_series_unlimited_py_c1_momentum_01; test -e _scratch/s11c/s15/repair4_series_unlimited_py_c1_momentum_01/outcome.json && echo 'outcome.json present' || echo 'outcome.json absent'
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
resource-samples.jsonl
stderr
stdout
outcome.json absent
```

```
$ ls _scratch/s15/reading_build_work/repair3/series_unlimited/ | /usr/bin/grep -F 'py_c1_momentum_01'
repair4_series_unlimited_py_c1_momentum_01
repair4_series_unlimited_py_c1_momentum_01_health.json
repair4_series_unlimited_py_c1_momentum_01_launch.json
repair4_series_unlimited_py_c1_momentum_01_launcher.log
```

```
$ tail -1 _scratch/s11c/s15/repair4_series_unlimited_py_c1_momentum_01/resource-samples.jsonl | python3 -c "import json,sys;j=json.loads(sys.stdin.read());print({k:j.get(k) for k in ('utc','memory.peak')})"
{'utc': '2026-11-06T05:41:07.391364+00:00', 'memory.peak': '1916407808'}
```

```
$ journalctl --user -u s11c-guard-350cfe7325e2.service -u s11c-guard-fd3dd4540a3c.service --since '2026-11-06 05:35:00 UTC' --until '2026-11-06 05:45:00 UTC' -o short-iso-precise --no-pager | cut -c1-200
2026-11-05T22:41:09.096152-07:00 pop-os systemd[1646]: Stopping s11c-guard-350cfe7325e2.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-05T22:41:09.106445-07:00 pop-os systemd[1646]: Stopped s11c-guard-350cfe7325e2.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-05T22:41:09.106605-07:00 pop-os systemd[1646]: s11c-guard-350cfe7325e2.service: Consumed 11h 38min 36.537s CPU time, 1.7G memory peak, 0B memory swap peak.
2026-11-05T22:42:51.396769-07:00 pop-os systemd[1646]: Stopping s11c-guard-fd3dd4540a3c.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-05T22:42:51.418999-07:00 pop-os systemd[1646]: Stopped s11c-guard-fd3dd4540a3c.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-05T22:42:51.419263-07:00 pop-os systemd[1646]: s11c-guard-fd3dd4540a3c.service: Consumed 1d 15h 46min 4.258s CPU time, 375.2M memory peak, 0B memory swap peak.
```

```
$ systemctl --user show -p LoadState,ActiveState s11c-guard-350cfe7325e2.service s11c-guard-fd3dd4540a3c.service
LoadState=not-found
ActiveState=inactive

LoadState=not-found
ActiveState=inactive
```

```
$ python3 -c "import json,sys;o=json.load(open('_scratch/s11c/s15/repair4_size_py_c1_rate_01/outcome.json'));print({k:o.get(k) for k in ('exitCode','wallSeconds')})"
{'exitCode': 1, 'wallSeconds': 143174.44159308192}
```

The guard registry and the kept copy.

```
$ cat /home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json; ls -l --time-style=+%Y-%m-%dT%H:%M:%SZ /home/trevnorris/.local/state/s11c-guard/resource-pool/; python3 -c "import json;r=json.load(open('/home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json.stale-20261106T054326Z'));[print(k,v['memoryMax'],v['owner']['pid'],v['logDirectory'].rsplit('/',1)[1]) for k,v in r.items()]"
{}
total 12
-rw-rw-r-- 1 trevnorris trevnorris   0 2026-10-15T21:09:59Z admission.lock
-rw-rw-r-- 1 trevnorris trevnorris   3 2026-11-05T22:43:26Z reservations.json
-rw-rw-r-- 1 trevnorris trevnorris 413 2026-10-28T18:28:58Z reservations.json.stale-20261029T005431Z
-rw-rw-r-- 1 trevnorris trevnorris 405 2026-11-05T22:42:51Z reservations.json.stale-20261106T054326Z
s11c-guard-350cfe7325e2.service 57982058496 1856816 repair4_series_unlimited_py_c1_momentum_01
```

```
$ echo "WolframKernel processes: $(pgrep -c WolframKernel)"; echo "processes named codex: $(pgrep -c -x codex)"; echo "guard units: $(systemctl --user list-units --no-legend 's11c-guard-*' | wc -l)"
WolframKernel processes: 0
processes named codex: 1
guard units: 0
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
 M research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4n_note_lookups.md
 M research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/S15_item4_singular_side_test.md
?? research/pde_ledger_v3/directives/_measurements/S15_item4_singular_side_test_lookups.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ echo "final run dirs after final_20: $(ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null | wc -l)"
final run dirs after final_20: 0
```
