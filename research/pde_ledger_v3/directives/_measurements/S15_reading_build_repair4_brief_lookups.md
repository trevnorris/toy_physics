# Measurements: S15 reading build repair-4 note, leak gate and citation check

Generator: `_scratch/s15/gen/s15_reading_repair4_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by the
orchestrator. The note resumes the repair-3 builder session, whose brief (gated at `0e51c6de`) stands.

## R0. The note

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_4.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_4.md; git log -1 --format='%h %s'
f9df2f12a63ca5fb8da9e4ebd803e13636a90e133bb7f2e0a0b1ab5efdb63eda  _scratch/s15/reading_build/resume_s15_reading_repair_4.md
5261
90441b89 Preserve the S15 reading build repair 3 as reviewed in round 3 (not accepted)
```

## R1. Values, exponents and methods from the review legs (fixed strings)

```
$ /usr/bin/grep -n -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'λ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
11:- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre tran
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
11:- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre tran
count: 1
```

```
$ /usr/bin/grep -n -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
7:Review confirmed that every condition your engines print is true, that no condition restates the question, that your run record's claims are accurate, and that items 1,
12:- Independent review reached exponent conditions, by exact routes checked by a second route, for representative objects of these classes, including the region momentum
15:- Wolfram: the identity attempt is `SolveAlways[Implies[dist > …, …], dist]` (`mathematica/S15_separation_departure.wl:519–522`). All 48 returns are the unevalua
21:- the conditions under which the item-3 object vanishes for every `D` on the domain, computed exactly from the printed closed form, complete, with every special expone
count: 4
```

```
$ /usr/bin/grep -n -F -- '+2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '+2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 0' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
11:- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre tran
15:- Wolfram: the identity attempt is `SolveAlways[Implies[dist > …, …], dist]` (`mathematica/S15_separation_departure.wl:519–522`). All 48 returns are the unevalua
count: 2
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
11:- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre tran
count: 1
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
11:- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre tran
count: 1
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b1' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '2a' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a+b' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a +' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a +' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q=' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q =' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'w=' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'S_k' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'S_k' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ẋ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ẋ' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
16:- SymPy: when a value holds `appellf1` or `polylog`, or more than one exponent form, the construction goes to `direct_attempt(coefficients,'all_parameters')` (`scripts
count: 1
```

```
$ /usr/bin/grep -n -F -- 'Log' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Log' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gauss' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Gauss' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'quadrature' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'quadrature' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pochhammer' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pochhammer' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'hypergeom' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'hypergeom' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'strata' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'strata' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'co-moving' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'co-moving' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'comoving' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'comoving' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'continuity' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'continuity' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'flux' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'flux' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-220
10:- Wolfram prints 210 numbered conditions on the 48 objects, and none of them names an exponent (`aa[·]` or `bb[·]`). Each is an amplitude, rate or velocity condition.
11:- SymPy prints 782 conditions on the 48 objects. Conditions naming an exponent appear in 18 objects: the inner and outer drain-1 groups of carried current, centre transport and radius transport. The other 30 objects h
16:- SymPy: when a value holds `appellf1` or `polylog`, or more than one exponent form, the construction goes to `direct_attempt(coefficients,'all_parameters')` (`scripts/S15_separation_departure.py:306–331`). That att
```

## R3. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
40:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
5:Your engines, transcripts and run record are committed as reviewed at 90441b89. The two transcripts are git-annex links again, with read-only targets: before a run writ
7:Review confirmed that every condition your engines print is true, that no condition restates the question, that your run record's claims are accurate, and that items 1,
24:- no printed line or run-record sentence states a property that no computation in that engine produced;
38:Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.
count: 4
```

```
$ /usr/bin/grep -n -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
1:Note from the orchestrator. Your repair has been reviewed. The review found item 4 still not delivered, and the user has chosen a further repair, to be done by you in t
5:Your engines, transcripts and run record are committed as reviewed at 90441b89. The two transcripts are git-annex links again, with read-only targets: before a run writ
7:Review confirmed that every condition your engines print is true, that no condition restates the question, that your run record's claims are accurate, and that items 1,
12:- Independent review reached exponent conditions, by exact routes checked by a second route, for representative objects of these classes, including the region momentum
count: 4
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
1:Note from the orchestrator. Your repair has been reviewed. The review found item 4 still not delivered, and the user has chosen a further repair, to be done by you in t
5:Your engines, transcripts and run record are committed as reviewed at 90441b89. The two transcripts are git-annex links again, with read-only targets: before a run writ
7:Review confirmed that every condition your engines print is true, that no condition restates the question, that your run record's claims are accurate, and that items 1,
12:- Independent review reached exponent conditions, by exact routes checked by a second route, for representative objects of these classes, including the region momentum
count: 4
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' _scratch/s15/reading_build/resume_s15_reading_repair_4.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' _scratch/s15/reading_build/resume_s15_reading_repair_4.md)"
count: 0
```

## R4. Operational sentences carried from the repair-3 notes

```
$ /usr/bin/grep -c -F -- 'polling' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:0
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:1
_scratch/s15/reading_build/foot_repair3.md:0
```

```
$ /usr/bin/grep -c -F -- 'queued notices' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:0
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:1
_scratch/s15/reading_build/foot_repair3.md:0
```

```
$ /usr/bin/grep -c -F -- 'long sleeps' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:0
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:1
_scratch/s15/reading_build/foot_repair3.md:0
```

```
$ /usr/bin/grep -c -F -- 'does not need a new go-ahead' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:0
```

```
$ /usr/bin/grep -c -F -- 'guard refuses admission' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:1
```

```
$ /usr/bin/grep -c -F -- 'outside item 4' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:0
```

```
$ /usr/bin/grep -c -F -- 'delete the link itself' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:0
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:1
```

```
$ /usr/bin/grep -c -F -- 'final guarded run' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:1
```

```
$ /usr/bin/grep -c -F -- 'own construction' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:2
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:1
```

```
$ /usr/bin/grep -c -F -- 'not-obtained label' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:0
_scratch/s15/reading_build/foot_repair3.md:1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_3b.md _scratch/s15/reading_build/resume_s15_reading_repair_3c.md _scratch/s15/reading_build/foot_repair3.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3b.md:1
_scratch/s15/reading_build/resume_s15_reading_repair_3c.md:1
_scratch/s15/reading_build/foot_repair3.md:1
```

## R5. Every line and count the note cites, at the reviewed commit

```
$ git show 90441b89:research/pde_ledger_v3/mathematica/S15_separation_departure.wl | sed -n '519,522p;530p;546p' | cut -c1-170
 query=Implies[dist>(rout[t]/.rules),And@@Thread[coefficients==0]];
 solo[prefix<>"/identity_solver_operand",query/.inverse];
 Block[{$MessageList={}},
  answer=Quiet[trackedSolveAlways[query,dist]/.SolveAlways->Inactive[SolveAlways]];messages=$MessageList];
 polynomialSlots=Pick[slots,Map[!MemberQ[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}],{rin[t],rout[t]}],#]&,physical]];
 Block[{$MessageList={}},monomialAnswer=Quiet[trackedReduce[monomialQuery,polynomialSlots,Complexes]/.Reduce->Inactive[Reduce]];messages=$MessageList];
```

```
$ git show 90441b89:research/pde_ledger_v3/scripts/S15_separation_departure.py | sed -n '306,308p;331p' | cut -c1-170
    all_values=[*coefficients,q]
    if any(v.has(S.appellf1,S.polylog) for v in all_values):
        direct_attempt(coefficients,'all_parameters');return answers,gaps
    if any(f.has(*remaining) for f in varying_forms):
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md && echo 'engines and run record unchanged since 90441b89'; ls -la research/pde_ledger_v3/scripts/out/S15_separation_departure.out research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | awk '{print $1, $9, $10, $11}' | cut -c1-140
engines and run record unchanged since 90441b89
lrwxrwxrwx research/pde_ledger_v3/mathematica/out/S15_separation_departure.out -> ../../../../.git/annex/objects/WQ/zZ/MD5E-s472539680--3ee0
lrwxrwxrwx research/pde_ledger_v3/scripts/out/S15_separation_departure.out -> ../../../../.git/annex/objects/gZ/F4/MD5E-s427868598--72ecbc78
```

```
$ sed -n '111743p' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-200
reading/region/component_1/momentum_labels_drain1/item4/all_parameters/solve_undetermined_coeffs/return: [{A_1(t): 0}, {B_1(t): 0}]
```

```
$ echo "wl conditions $(grep -c '/item4/condition_[0-9]*: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out); with aa[ or bb[ $(grep '/item4/condition_[0-9]*: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -c -E 'aa\[|bb\[')"; echo "wl Inactive SolveAlways $(grep -c 'Inactive\[SolveAlways\]' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out)"; echo "py all_parameters $(grep -c "unresolved_parameter_branches: \[('all_parameters'" research/pde_ledger_v3/scripts/out/S15_separation_departure.out)"
wl conditions 210; with aa[ or bb[ 0
wl Inactive SolveAlways 48
py all_parameters 30
```

```
$ grep -n -E 'conditions 782|conditions 210|with a condition' research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r3_lookups.md | cut -c1-200
204:$ echo "numbered item-4 condition lines in the transcript: $(grep -c '/item4/condition_[0-9]*: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out)"; awk '{c+=$1; e+=$2; if ($1>0
206:48 objects: conditions 210 objects with a condition 48 conditions naming aa[ or bb[ 0 objects with such a condition 0
321:$ echo "numbered item-4 condition lines in the transcript: $(grep -c '/item4/condition_[0-9]*: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out)"; awk '{c+=$1; e+=$2; if ($1>0) n+
323:48 objects: conditions 782 objects with a condition 48 conditions naming an exponent symbol 596 objects with such a condition 18
```
