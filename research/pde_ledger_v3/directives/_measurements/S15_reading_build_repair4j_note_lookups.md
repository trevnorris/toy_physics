# Measurements: S15 reading build repair-4j note (user: run the first two to completion), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit
is read by the orchestrator. The run and work directories named in these commands are local scratch and disposable;
every output this note relies on is quoted here in full. The note carries the user's decision after the trial of a
cheaper item-4 construction (`resume_s15_reading_repair_4h.md`, gated at `d371f16c`; `resume_s15_reading_repair_4i.md`,
gated at `4d2c29fc`) ended with no attempt returning within 1 hour: two objects are attempted in both engines with no time
stop, as a measurement, not a final run. The final runs keep the rule of `resume_s15_reading_repair_4g.md` (`1a703b67`).
The note names no value, exponent or condition.

## The user's words and decision (2026-11-04)

AskUserQuestion, ~11:00Z. Question: "How should the final runs treat the six whole-region objects that carry both drains?
None returned in 1 hour with the new method, and Wolfram mostly didn't return within 6 hours with the old one." Options
offered: "Shorter stop for the six (Recommended)", "Full 6-hour rule", "Rethink item 4 first". The user answered in their
own words: **"Allow the first two to run to completion, regardless of how long. We should understand the scope of what
we're dealing with."**

The orchestrator reads "the first two" as the first two of the six objects in the order the trial named them,
`region/component_1/momentum_labels_{drain1, drain2}` and `region/component_1/rate_labels_{drain1, drain2}`, in both
engines, ending only on a return, running out of memory, or an error in the builder's code.

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
94a862c71a6cb10375e6841b86e4a7f1a421e9f7e2cb6e530823aee329e85a59  research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
2708 research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
22:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md)"
18:- Run through the guard. Run the two SymPy attempts alongside one Wolfram attempt at a time, in the order the objects are named above. While jobs share the machine, size each cap from that object's
count: 1
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4f.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4g.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4h.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4i.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'method_trial_report.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'method_trial_attempt_inventory.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
2
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
2
```

```
$ /usr/bin/grep -c -F -- 'No attempt is stopped for time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'No production source changes' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'Poll each job until it ends' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not start the final runs' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
1
```

```
$ ls _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
```

```
$ /usr/bin/grep -n -F 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:18:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:7:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:9:Everything under "What must be true" in `resume_s15_reading_repair_4.md` stands. For an object whose attempt ran out of memory, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md` and in `resume_s15_reading_repair_4d.md` stands. For an object whose attempt was stopped at 6 hours, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md` and `resume_s15_reading_repair_4f.md` stands. This decision replaces the bullet in `resume_s15_reading_repair_4f.md` about an attempt that goes silent. Also:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md` and `resume_s15_reading_repair_4g.md` stands, with these changes and additions:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md`, `resume_s15_reading_repair_4f.md`, `resume_s15_reading_repair_4g.md` and `resume_s15_reading_repair_4h.md` stands. Also:
```

```
$ /usr/bin/grep -n -F -e '1 hour after it started' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-220
9:The trial measures how fast your construction is. The user decided that in this trial, an attempt is stopped 1 hour after it started and is reported as not returned within 1 hour. The Wolfram trial of `region/component
16:- No trial outcome is printed as an item-4 label. A trial attempt ends when it returns its conditions, when it runs out of memory, or 1 hour after it started.
23:- Run through the guard, one job at a time, with `--memory-gib 55`, and poll each job until it ends. Each attempt stops 1 hour after it started. If an attempt ends on an error in your code, you may fix it and run that
```

```
$ /usr/bin/grep -n -F -e 'only when it runs out of memory' -e '6 hours after it started' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-220
9:`resume_s15_reading_repair_4f.md` said that an item-4 attempt that goes silent for a long stretch is labelled, and it did not say how long. That sentence was mine; the user's decision was the 6-hour stop. The user has 
```

```
$ /usr/bin/grep -n -F -e 'region/component_1/momentum_labels_{drain1, drain2}' -e 'region/component_1/rate_labels_{drain1, drain2}' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-140
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md:24:- Try the construction on the six objects `region/component_1/momentum_
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md:5:On 2026-11-03 I stopped `repair4_scalar_wl_c1_momentum_03` through its u
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md:9:The trial measures how fast your construction is. The user decided that 
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md:21:  - in Wolfram first, `region/component_1/rate_labels_{drain1, drain2}`
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md:22:  - then all six in SymPy, including `region/component_1/momentum_label
```

## R5. Facts the note states

```
$ cat _scratch/s15/reading_build/repair4i.DONE; ls -l _scratch/s15/reading_build_work/repair3/method_trial/method_trial_report.md _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md | awk '{print $5, $9}'; sha256sum _scratch/s15/reading_build_work/repair3/method_trial/method_trial_report.md _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md _scratch/s15/reading_build_work/repair3/method_trial/method_trial_manifest.json
exit 0
14410 _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md
3774 _scratch/s15/reading_build_work/repair3/method_trial/method_trial_report.md
540db289b181f3069e97091c30baebd70bedb33e95ab6030dfc137a9be75e7d3  _scratch/s15/reading_build_work/repair3/method_trial/method_trial_report.md
649ff683b1a1cf9bebea87330672d9a338cd925ec09d8b298b4d0c0bca917bdd  _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md
ee8a6f86e791dc5559cb91b7a88bb52742d20f6e7bebc09171acf6590c42d81d  _scratch/s15/reading_build_work/repair3/method_trial/method_trial_manifest.json
```

```
$ /usr/bin/grep -n -F -e 'did not return complete conditions within one hour' -e 'None of these trial runs exhausted memory' -e 'No final runs were started' _scratch/s15/reading_build_work/repair3/method_trial/method_trial_report.md | cut -c1-260
3:**The construction did not return complete conditions within one hour for any of the twelve engine/object pairs.** None of these trial runs exhausted memory. The five newly timed Wolfram attempts and six corrected SymPy attempts reached their one-hour cutoff
16:All created/changed trial scripts and their SHA256 values are listed in the inventory; [artifact checksums](method_trial_artifacts.sha256) also cover this report and its manifests. Production engine/launcher hashes match the captured pre-trial hashes. Revie
```

```
$ /usr/bin/grep -c -F '| not_returned_within_one_hour |' _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md; /usr/bin/grep -c -F -e '| conditions |' -e '| returned |' _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md
12
0
```

```
$ /usr/bin/grep -F -e 'c1_momentum' -e 'c1_rate' _scratch/s15/reading_build_work/repair3/method_trial/method_trial_attempt_inventory.md | /usr/bin/grep -F '| 20' | cut -d'|' -f2,3,5,6,7,8
 `repair4_scalar_wl_c1_momentum_01` | Wolfram c1_momentum | 58.115654 | 78.0212976608891 | 405708800 | code_error_author_stop 
 `repair4_scalar_wl_c1_momentum_02` | Wolfram c1_momentum | 4554.828131 | 4574.71407812892 | 3492872192 | construction_error_author_stop 
 `repair4_scalar_wl_c1_momentum_03` | Wolfram c1_momentum | 6180.121138 | 6199.971640525968 | 23163592704 | not_returned_within_one_hour 
 `repair4_hour_wl_c1_rate_01` | Wolfram c1_rate | 3601.4597240389558 | 3673.0518820240395 | 13472899072 | not_returned_within_one_hour 
 `repair4_hour_py_c1_momentum_01` | SymPy c1_momentum | not started | 397.49822601804044 | 185503744 | code_error_no_selected_object 
 `repair4_hour_py_c1_rate_01` | SymPy c1_rate | not started | 391.64498175296467 | 164569088 | code_error_no_selected_object 
 `repair4_hour_py_c1_momentum_02` | SymPy c1_momentum | 2749.2573616049485 | 3152.7048533839406 | 650948608 | construction_error 
 `repair4_hour_py_c1_momentum_03` | SymPy c1_momentum | 3600.3420180459507 | 4002.123119862983 | 596459520 | not_returned_within_one_hour 
 `repair4_hour_py_c1_rate_02` | SymPy c1_rate | 3600.343907054048 | 4594.864452120033 | 321032192 | not_returned_within_one_hour 
 `repair4_scalar_wl_c1_momentum_01` | 2026-11-03T18:04:10.286865+00:00 | WORKING_ON: "reading/region/component_1/momentum_labels_{drain1, drain2} / closed geometric differentiation" |
 `repair4_scalar_wl_c1_momentum_02` | 2026-11-03T18:05:41.746030+00:00 | WORKING_ON: "reading/region/component_1/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / scalar quotient 1 / FullSimplify substituted functional range" |
 `repair4_scalar_wl_c1_momentum_03` | 2026-11-03T19:22:48.778857+00:00 | WORKING_ON: "reading/region/component_1/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / scalar quotient 1 / eliminate series pivot 29" |
 `repair4_hour_wl_c1_rate_01` | 2026-11-03T21:11:02.481843+00:00 | WORKING_ON: reading/region/component_1/rate_labels_{drain1, drain2} / item4 exponent-sensitive / scalar quotient 6 / branch_41 PiecewiseExpand Q and departure strata |
 `repair4_hour_py_c1_momentum_01` | None | WORKING_ON: reading family and derivative operator / solve symbolic regular reference at source line 2279 |
 `repair4_hour_py_c1_rate_01` | None | WORKING_ON: reading family and derivative operator / solve symbolic regular reference at source line 2279 |
 `repair4_hour_py_c1_momentum_02` | 2026-11-04T03:03:06.494878+00:00 | reading/region/component_1/momentum_labels_drain1+drain2/item4 / sparse_scalar/coefficient_0/polynomial_branch_44/solve_quadratic_parameter |
 `repair4_hour_py_c1_momentum_03` | 2026-11-04T03:56:19.063370+00:00 | reading/region/component_1/momentum_labels_drain1+drain2/item4 / sparse_scalar/coefficient_0/polynomial_branch_45/solve_radius_ratio_algebraic_roots |
 `repair4_hour_py_c1_rate_02` | 2026-11-04T05:12:57.219650+00:00 | reading/region/component_1/rate_labels_drain1+drain2/item4 / sparse_scalar/coefficient_0/primitive_stratum_0/node_1/coefficient_row_rank_16 |
```

```
$ for r in repair4_scalar_wl_c1_momentum_0{1,2,3} repair4_hour_wl_{c1_rate,c2_momentum,c2_rate,c3_momentum,c3_rate}_01 repair4_hour_py_{c1_momentum,c1_rate,c2_momentum,c2_rate,c3_momentum,c3_rate}_01 repair4_hour_py_c1_momentum_0{2,3} repair4_hour_py_{c1_rate,c2_momentum,c2_rate,c3_momentum,c3_rate}_02; do printf '%s ' $r; tr -d '\n ' < _scratch/s11c/s15/$r/outcome.json | cut -c1-120; echo; done
repair4_scalar_wl_c1_momentum_01 {"exitCode":1,"wallSeconds":78.0212976608891,"unit":"s11c-guard-ed6ed80e6ab7","stderrBytes":0,"limitsVerified":true,"chi

repair4_scalar_wl_c1_momentum_02 {"exitCode":1,"wallSeconds":4574.71407812892,"unit":"s11c-guard-a72a6f8c5032","stderrBytes":0,"limitsVerified":true,"chi

repair4_scalar_wl_c1_momentum_03 {"exitCode":1,"wallSeconds":6199.971640525968,"unit":"s11c-guard-a21d67f54053","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_wl_c1_rate_01 {"exitCode":0,"wallSeconds":3673.0518820240395,"unit":"s11c-guard-83f334261db2","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_wl_c2_momentum_01 {"exitCode":0,"wallSeconds":3619.3081173480023,"unit":"s11c-guard-4591b5d05c73","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_wl_c2_rate_01 {"exitCode":0,"wallSeconds":3670.0726244989783,"unit":"s11c-guard-5e4dc31c4ffd","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_wl_c3_momentum_01 {"exitCode":0,"wallSeconds":3619.275505815982,"unit":"s11c-guard-efd6367ddf05","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_wl_c3_rate_01 {"exitCode":0,"wallSeconds":3670.0597738939105,"unit":"s11c-guard-ac3b95a3359b","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_py_c1_momentum_01 {"exitCode":0,"wallSeconds":397.49822601804044,"unit":"s11c-guard-42e0a4538ba5","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_py_c1_rate_01 {"exitCode":0,"wallSeconds":391.64498175296467,"unit":"s11c-guard-df8f34eee679","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_py_c2_momentum_01 {"exitCode":0,"wallSeconds":392.444408598938,"unit":"s11c-guard-4f542667bea6","stderrBytes":0,"limitsVerified":true,"chi

repair4_hour_py_c2_rate_01 {"exitCode":0,"wallSeconds":394.0546806290513,"unit":"s11c-guard-588a188b83d7","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c3_momentum_01 {"exitCode":0,"wallSeconds":392.6414507990703,"unit":"s11c-guard-2ca231a49ad6","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c3_rate_01 {"exitCode":0,"wallSeconds":393.9424918759614,"unit":"s11c-guard-fcb61760700c","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c1_momentum_02 {"exitCode":2,"wallSeconds":3152.7048533839406,"unit":"s11c-guard-7e262183a11d","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_py_c1_momentum_03 {"exitCode":0,"wallSeconds":4002.123119862983,"unit":"s11c-guard-aa47dd5ebef3","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c1_rate_02 {"exitCode":0,"wallSeconds":4594.864452120033,"unit":"s11c-guard-bc23ba10ea25","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c2_momentum_02 {"exitCode":0,"wallSeconds":4008.262891040067,"unit":"s11c-guard-0f71c29c37ea","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c2_rate_02 {"exitCode":0,"wallSeconds":4460.512710707146,"unit":"s11c-guard-6f020549008e","stderrBytes":0,"limitsVerified":true,"ch

repair4_hour_py_c3_momentum_02 {"exitCode":0,"wallSeconds":4003.0028028779197,"unit":"s11c-guard-32fa2795faff","stderrBytes":0,"limitsVerified":true,"c

repair4_hour_py_c3_rate_02 {"exitCode":0,"wallSeconds":4449.721750638913,"unit":"s11c-guard-1cf3e8cde009","stderrBytes":0,"limitsVerified":true,"ch

```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c -f 'codex exec'
{}

0
1
```

```
$ sha256sum research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f  research/pde_ledger_v3/scripts/S15_separation_departure.py
4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1  research/pde_ledger_v3/mathematica/S15_separation_departure.wl
b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ python3 -c "import json;d=json.load(open('_scratch/s15/reading_build_work/repair3/method_trial/method_trial_manifest.json'));print(json.dumps(d['production_source_hash_check'],indent=1))"
{
 "research/pde_ledger_v3/mathematica/S15_separation_departure.wl": {
  "expected": "4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1",
  "actual": "4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1",
  "unchanged": true
 },
 "research/pde_ledger_v3/mathematica/S15_separation_departure_run.py": {
  "expected": "b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302",
  "actual": "b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302",
  "unchanged": true
 },
 "research/pde_ledger_v3/scripts/S15_separation_departure.py": {
  "expected": "48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f",
  "actual": "48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f",
  "unchanged": true
 }
}
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md; echo "run record diff vs 90441b89 exit $?"; git status --short -- research/pde_ledger_v3
run record diff vs 90441b89 exit 0
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4j_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null; echo "final run dirs after final_20: $(ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null | wc -l)"
final run dirs after final_20: 0
```

```
$ cat _scratch/s11c/s15/repair4_hour_py_c3_rate_02/duration-validation.json
{
  "verified": true,
  "actual": {
    "RuntimeMaxUSec": "infinity",
    "Restart": "no"
  },
  "wallDeadlineSeconds": null,
  "inactivityDeadlineSeconds": null
}
```
