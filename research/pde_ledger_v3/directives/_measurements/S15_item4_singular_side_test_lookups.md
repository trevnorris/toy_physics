# Measurements: S15 item-4 Singular side-test directive, leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts and hashes. Lines cut
to a byte width can end inside a character; the record drops such partial characters. Every hit is read by the
orchestrator. The run directories named here are local scratch and disposable; every output relied on is quoted here.
The directive is orchestrator-written and goes to a fresh Codex builder, so it gets two review legs before launch.

## The user's decision (2026-11-06)

After the orchestrator suggested a Singular trial, the AskUserQuestion "Which engine should we try for heavy algebra
alongside Mathematica?" got the answer "Not yet". The user then wrote, in their words: "Doing a side test with Singular
using codex doesn't sound like a bad idea as well to get another data point". Singular was installed at the user's
request (their message was the install command) under /var/projects/_tools, without system changes.

## R0. The directive

```
$ sha256sum research/pde_ledger_v3/directives/S15_item4_singular_side_test.md; wc -c research/pde_ledger_v3/directives/S15_item4_singular_side_test.md
fbc7cbd689a3e75f0847c887610ac00adc1f84543c947f922db3671b11e4294b  research/pde_ledger_v3/directives/S15_item4_singular_side_test.md
4481 research/pde_ledger_v3/directives/S15_item4_singular_side_test.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the directive

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
7:Item 4 of `research/pde_ledger_v3/directives/S15_reading_decision_list.md` asks, for each departure object of item 3 that depends on the separation, for the exact and complete conditions under which
11:`region/component_1/momentum_labels_{drain1, drain2}`, as items 3 and 4 of the reading decision list define it, on the region and domain of `research/pde_ledger_v3/directives/S15_decision_list.md`.
15:The section "Supplied (unfalsifiable in this trial)" of `research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md`: the distance from drain 1 over the region, and its series in `z = r/d`. 
count: 3
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gegenbauer' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Gegenbauer' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[|hh[12])' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[|hh[12])' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

## R3. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
41:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
39:This task ends at your report; it includes no review.
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
39:This task ends at your report; it includes no review.
count: 1
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Maple' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Maple' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'expected' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'expected' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'should equal' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'should equal' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'match' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'match' research/pde_ledger_v3/directives/S15_item4_singular_side_test.md)"
count: 0
```

## R4. Pointers

```
$ ls research/pde_ledger_v3/directives/S15_reading_decision_list.md research/pde_ledger_v3/directives/S15_decision_list.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl scripts/s11c_guarded_run.py
research/pde_ledger_v3/directives/S15_decision_list.md
research/pde_ledger_v3/directives/S15_reading_decision_list.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
research/pde_ledger_v3/mathematica/S15_separation_departure.wl
research/pde_ledger_v3/scripts/S15_separation_departure.py
scripts/s11c_guarded_run.py
```

```
$ /usr/bin/grep -n -F 'Supplied (unfalsifiable in this trial)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
15:## Supplied (unfalsifiable in this trial)
```

```
$ git log -1 --format='%h %s' -- research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
e5e5a220 Update the S15 repair-4l note: the SymPy momentum measurement attempt has ended
```

```
$ ls -l /var/projects/_tools/singular/bin/Singular; /var/projects/_tools/singular/bin/Singular --version </dev/null 2>&1 | head -3
-rwxrwxr-x 2 trevnorris trevnorris 22272 Feb  9  2026 /var/projects/_tools/singular/bin/Singular
Singular for x86_64-Linux version 4.4.1 (44100, 64 bit) Jan 2025
with
	GMP(6.3.0),NTL(11.6.0),FLINT(3.3.1),
```

```
$ ls /var/projects/_tools/singular/share/singular/LIB/ | wc -l; ls /var/projects/_tools/singular/share/singular/LIB/grobcov.lib /var/projects/_tools/singular/share/singular/LIB/primdec.lib
219
/var/projects/_tools/singular/share/singular/LIB/grobcov.lib
/var/projects/_tools/singular/share/singular/LIB/primdec.lib
```

```
$ python3 -c 'import sympy; print(sympy.__version__)'
1.14.0
```

```
$ /usr/bin/grep -n -F -e '--pool' scripts/s11c_guarded_run.py | head -5 | cut -c1-200
8:Defaults remain 2 GiB and 32 tasks. Opt-in --pool jobs share one machine-wide
452:    parser.add_argument('--pool', help='opt in to shared machine-wide admission; named purpose')
468:        parser.error('--pool must be a short alphanumeric purpose name')
470:        parser.error('--all-cpus-exclusive requires --pool')
```

```
$ ls -d _scratch/singular_side_test 2>&1
ls: cannot access '_scratch/singular_side_test': No such file or directory
```

## R5. Facts the directive states

The engines' output name of the object (progress lines of the last SymPy and Wolfram attempts on it).

```
$ echo 'SymPy:'; /usr/bin/grep -a -m1 -o -F 'region/component_1/momentum_labels_drain1+drain2' _scratch/s11c/s15/repair4_series_unlimited_py_c1_momentum_01/stdout; echo 'Wolfram:'; /usr/bin/grep -a -m1 -o -F 'region/component_1/momentum_labels_{drain1, drain2}' _scratch/s11c/s15/repair4_series_unlimited_wl_c1_momentum_01/stdout; echo 'Wolfram, the SymPy spelling, count:'; /usr/bin/grep -a -c -F 'region/component_1/momentum_labels_drain1+drain2' _scratch/s11c/s15/repair4_series_unlimited_wl_c1_momentum_01/stdout
SymPy:
region/component_1/momentum_labels_drain1+drain2
Wolfram:
region/component_1/momentum_labels_{drain1, drain2}
Wolfram, the SymPy spelling, count:
0
```

No exact construction in either engine has returned the conditions of the six objects: the evidence is quoted in
`S15_reading_build_repair4n_note_lookups.md` R5 (note 4h line 9, the method trial report line 3, the scope measurements,
the series trial report lines 14, 16 and 20, the long-run Wolfram result, and the stopped SymPy long attempt).

```
$ sha256sum research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4n_note_lookups.md
f7475dca1ebb514024cfe25d58300c2f5aaf7d0b4b93c1cba275807ff9720d98  research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4n_note_lookups.md
```
