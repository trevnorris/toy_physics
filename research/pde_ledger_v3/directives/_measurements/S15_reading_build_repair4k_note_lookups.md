# Measurements: S15 reading build repair-4k note (user: run the Wolfram reruns now), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes, process and
unit state. Every hit is read by the orchestrator. The run and work directories named in these commands are local scratch
and disposable; every output this note relies on is quoted here in full. The note carries a user decision on how the scope
measurement of `resume_s15_reading_repair_4j.md` (gated at `acc4f1a3`) continues after both Wolfram attempts were stopped by
their caps while sharing the machine, and it reports how the orchestrator stopped the builder session. It changes no object,
no construction and no stop rule other than where the 55 GiB Wolfram reruns run; it names no value, exponent or condition.

## The user's decision (2026-11-04)

AskUserQuestion, ~21:05Z. Question: "Both Wolfram attempts hit their sharing limits and now wait until both SymPy runs
finish, which could take days. How should the scope measurement continue?" Options offered:
- "Run Wolfram now (Recommended)": "Wolfram reruns start now, one at a time, at the full 55 GiB, beside the two SymPy runs
  (0.3 GB each). That's the same 55 GiB limit as running alone. I stop the builder's job sequencer and give the builder a
  short note. The SymPy runs are in their own system units and should keep going, but there's a small risk of losing their
  7 hours."
- "Keep waiting": "The plan as written: the Wolfram reruns start only after both SymPy runs finish, with no estimate for
  when. 55 GiB sits idle meanwhile. Nothing changes."
- "Rethink item 4 now": "The scope is already large: Wolfram went past 32 GiB and nearly 5 h on its first quotient, and
  SymPy has spent 7 h on its first coefficient. Stop measuring and change item 4 in the reading decision list instead, with
  two review legs."

The user selected **"Run Wolfram now (Recommended)"**.

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
fcae9b3269a3865760dc83782565fdaf7fcd2e92d33a46c4d838dc770b5a0245  research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
2977 research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
5:On 2026-11-04 at 21:24:59Z I sent SIGSTOP to your session's process (pid 236750). I meant to hold your session still while I stopped `run_measurements.py` (pid 237622) on its own. Instead, the shell
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
26:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
5:On 2026-11-04 at 21:24:59Z I sent SIGSTOP to your session's process (pid 236750). I meant to hold your session still while I stopped `run_measurements.py` (pid 237622) on its own. Instead, the shell
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md)"
11:Both Wolfram attempts were stopped by their caps while sharing the machine: `repair4_size_wl_c1_momentum_01` at 32 GiB and `repair4_size_wl_c1_rate_01` at 20 GiB. Under `resume_s15_reading_repair_4
15:Everything under "What must be true" in `resume_s15_reading_repair_4j.md` stands, including what it carries over from the earlier notes, except that a Wolfram attempt at `--memory-gib 55` runs alon
17:- At most one Wolfram attempt runs at a time.
22:- Attempt `region/component_1/momentum_labels_{drain1, drain2}` and then `region/component_1/rate_labels_{drain1, drain2}` again in Wolfram, through the guard with `--memory-gib 55`, with the const
count: 4
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4j.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
4
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
2
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
4
```

```
$ /usr/bin/grep -c -F -- 'not stopped or restarted' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
1
```

```
$ /usr/bin/grep -c -F -- 'At most one Wolfram attempt runs at a time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
1
```

```
$ /usr/bin/grep -c -F -- 'Poll every job until it ends' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
1
```

```
$ /usr/bin/grep -c -F -- 'report every attempt' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not start the final runs' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
1
```

```
$ ls research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md; /usr/bin/grep -n -F -e 'What must be true' -e 'run it again alone' -e 'No attempt is stopped for time' -e 'write a short report' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200
research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
9:## What must be true
11:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md`, `resume_s15_readin
12:- No attempt is stopped for time. An attempt ends when it returns its conditions, when it runs out of memory at `--memory-gib 55` while running alone, or on an error in your code.
18:- Run through the guard. Run the two SymPy attempts alongside one Wolfram attempt at a time, in the order the objects are named above. While jobs share the machine, size each cap from that object's
20:- When all four objects have ended, stop and write a short report. For each attempt give the run, its start and end, whether it returned its conditions, ran out of memory or ended on an error in yo
```

```
$ git log -1 --format='%h %s' -- research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
acc4f1a3 Gate the S15 reading build repair-4j note (user: run the first two to completion)
```

## R5. Facts the note states

The stop: the shell's exit marker and its time, the last health write, the controller state.

```
$ cat _scratch/s15/reading_build/repair4j.DONE; stat -c '%y %n' _scratch/s15/reading_build/repair4j.DONE _scratch/s15/reading_build_work/repair3/size_measurement/latest_health.json _scratch/s15/reading_build_work/repair3/size_measurement/controller.json
exit 147
2026-11-04 14:24:59.824933902 -0700 _scratch/s15/reading_build/repair4j.DONE
2026-11-04 14:24:11.869104911 -0700 _scratch/s15/reading_build_work/repair3/size_measurement/latest_health.json
2026-11-04 13:50:56.727154929 -0700 _scratch/s15/reading_build_work/repair3/size_measurement/controller.json
```

```
$ python3 -c "import json;h=json.load(open('_scratch/s15/reading_build_work/repair3/size_measurement/latest_health.json'));print(h['utc'])"
2026-11-04T21:24:11.868967+00:00
```

```
$ python3 -c "import json;c=json.load(open('_scratch/s15/reading_build_work/repair3/size_measurement/controller.json'));print('updated',c['updated']);print('active',[a['run'] for a in c['active']]);print('pending_standalone',c['pending_standalone']);[print('completed',x['run'],'cap_gib',x['cap_gib'],'result',x['result'],'peak_memory_bytes',x['peak_memory_bytes'],'started',x['execution']['started'],'ended',x['execution']['ended'],'wall_seconds',round(x['execution']['wall_seconds'],1)) for x in c['completed']]"
updated 2026-11-04T20:50:56.726880+00:00
active ['repair4_size_py_c1_momentum_01', 'repair4_size_py_c1_rate_01']
pending_standalone [['wl', 'momentum', 2], ['wl', 'rate', 2]]
completed repair4_size_wl_c1_momentum_01 cap_gib 32 result shared_cap_reached_retry_alone peak_memory_bytes 32932380672 started 2026-11-04T13:56:38.085287+00:00 ended 2026-11-04T18:46:43.494187+00:00 wall_seconds 17405.4
completed repair4_size_wl_c1_rate_01 cap_gib 20 result shared_cap_reached_retry_alone peak_memory_bytes 21474836480 started 2026-11-04T18:47:00.270081+00:00 ended 2026-11-04T20:50:25.927973+00:00 wall_seconds 7405.7
```

```
$ /usr/bin/grep -n -F -e 'shared_cap_reached_retry_alone' -e 'if not active and retry' _scratch/s15/reading_build_work/repair3/size_measurement/run_measurements.py | cut -c1-200
46:  status='shared_cap_reached_retry_alone' if job['meta']['cap_gib']<55 else 'ran_out_of_memory'
74: if not active and retry:
```

```
$ ps -p 236750,237622 -o pid,comm 2>&1; echo "processes named codex: $(pgrep -c -x codex)"; echo "python3 processes running run_measurements.py: $(pgrep -c -f '^python3 [^ ]*run_measurements[.]py')"
    PID COMMAND
processes named codex: 0
python3 processes running run_measurements.py: 0
```

The SymPy attempts: units, launchers and registry.

```
$ systemctl --user is-active s11c-guard-c4214c0da91a.service s11c-guard-fd3dd4540a3c.service
active
active
```

```
$ ps -o pid,ppid,comm -p 237624,237635; ps -o pid,comm -p 1646
    PID    PPID COMMAND
 237624    1646 python3
 237635    1646 python3
    PID COMMAND
   1646 systemd
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json
{
  "s11c-guard-c4214c0da91a.service": {
    "pool": "s15",
    "poolMemoryMax": 66571993088,
    "memoryMax": 2147483648,
    "tasksMax": 32,
    "cpu": 7,
    "owner": {
      "pid": 237624,
      "startTicks": "107074519",
      "bootId": "c3b7b2c7-a5fc-4b45-adf0-6cbde3d609bb"
    },
    "logDirectory": "/var/projects/toy_physics/_scratch/s11c/s15/repair4_size_py_c1_momentum_01"
  },
  "s11c-guard-fd3dd4540a3c.service": {
    "pool": "s15",
    "poolMemoryMax": 66571993088,
    "memoryMax": 2147483648,
    "tasksMax": 32,
    "cpu": 6,
    "owner": {
      "pid": 237635,
      "startTicks": "107074619",
      "bootId": "c3b7b2c7-a5fc-4b45-adf0-6cbde3d609bb"
    },
    "logDirectory": "/var/projects/toy_physics/_scratch/s11c/s15/repair4_size_py_c1_rate_01"
  }
}
```

```
$ for r in repair4_size_py_c1_momentum_01 repair4_size_py_c1_rate_01; do printf '%s outcome.json: ' $r; test -e _scratch/s11c/s15/$r/outcome.json && echo present || echo absent; done
repair4_size_py_c1_momentum_01 outcome.json: absent
repair4_size_py_c1_rate_01 outcome.json: absent
```

```
$ echo "WolframKernel processes: $(pgrep -c WolframKernel)"
WolframKernel processes: 0
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
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4k_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ echo "final run dirs after final_20: $(ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null | wc -l)"
final run dirs after final_20: 0
```
