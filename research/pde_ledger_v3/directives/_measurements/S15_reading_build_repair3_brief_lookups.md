# Measurements: S15 reading build repair-3 brief, leak gate and citation check

Generator: `_scratch/s15/gen/s15_reading_repair3_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by the
orchestrator. Probes run over the whole brief and over its new parts (the head and foot), since the contract between
the rules is copied byte for byte from the cleared lists.

## R0. The brief

```
$ sha256sum _scratch/s15/reading_build/brief_s15_reading_repair_3.md _scratch/s15/reading_build/head_repair3.md _scratch/s15/reading_build/foot_repair3.md _scratch/s15/reading_build/assemble_repair3.sh; wc -c < _scratch/s15/reading_build/brief_s15_reading_repair_3.md; git log -1 --format='%h %s'
141c4f818420129692feb09fb52a15b064d89d819fce3b545526bab3caa33fcf  _scratch/s15/reading_build/brief_s15_reading_repair_3.md
bdf86b8281d003bbae9b08dba49311ecc918a59ae91d427ec2b35efba3bbdb0f  _scratch/s15/reading_build/head_repair3.md
92b6d9eb84183d95963a617d60e0e84241d3f611e87db4aa87a19aebabdb2645  _scratch/s15/reading_build/foot_repair3.md
c96eb48a8e7e2f8623294d361840dd07c51cdf20b2e531af1f29b1fd31a41365  _scratch/s15/reading_build/assemble_repair3.sh
22059
4c4e2558 Preserve the S15 reading build repair 2 as reviewed in round 2 (not accepted)
```

```
$ diff <(sed -n '/^# Part A/,/^-\{100\}$/p' _scratch/s15/reading_build/brief_s15_reading_build.md) <(sed -n '/^# Part A/,/^-\{100\}$/p' _scratch/s15/reading_build/brief_s15_reading_repair_3.md) && echo 'contract parts byte-identical to the original gated brief (b6e6f025)'
contract parts byte-identical to the original gated brief (b6e6f025)
```

```
$ git show b6e6f025 --stat --format='%h %s' | head -3
b6e6f025 Gate the S15 reading build brief (leak probe and contract check)

 .../S15_reading_build_directive_lookups.md         | 366 +++++++++++++++++++++
```

## R1. Values, exponents and methods from the review legs (fixed strings, in the new parts)

```
$ /usr/bin/grep -n -F -- 'alpha' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'λ' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lam' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lam' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
18:- There are 48 objects whose item-3 closed form depends on the separation `d`: three components each of the drain-1 and mixed groups of the carried current, the centre
count: 1
```

```
$ /usr/bin/grep -n -F -- '= -1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= -1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' _scratch/s15/reading_build/r3_gate_new_parts.md)"
18:- There are 48 objects whose item-3 closed form depends on the separation `d`: three components each of the drain-1 and mixed groups of the carried current, the centre
27:- For each of the 48 objects, in each engine, the conditions under which the item-3 object vanishes for every `D` on the domain:
count: 2
```

```
$ /usr/bin/grep -n -F -- '−4' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−4' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-4' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-4' _scratch/s15/reading_build/r3_gate_new_parts.md)"
18:- There are 48 objects whose item-3 closed form depends on the separation `d`: three components each of the drain-1 and mixed groups of the carried current, the centre
count: 1
```

```
$ /usr/bin/grep -n -F -- '−5' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−5' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-5' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-5' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' _scratch/s15/reading_build/r3_gate_new_parts.md)"
4:- Part A is from research/pde_ledger_v3/directives/S15_decision_list.md at 91248b67: the sections "Supplied equations", "The drains, their fields and the region", "What
5:- Part B is from research/pde_ledger_v3/directives/S15_reading_decision_list.md at 7f52b4d2: the sections "Supplied (unfalsifiable in this build)", "The separation and 
14:Their run record is research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md. All five files are committed at 4c4e2558 as reviewed. Line 
58:Finish with the final guarded run of each engine, from the final code. Then update research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.
count: 4
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' _scratch/s15/reading_build/r3_gate_new_parts.md)"
19:  - Each engine prints, as the condition, a finite sum of squared terms set to zero. In SymPy this is `finite_identity_condition` (`scripts/S15_separation_departure.py
count: 1
```

```
$ /usr/bin/grep -n -F -- 'a_1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '2a' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a+b' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a +' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a +' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'A_1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'A_1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'B_1' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'B_1' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'A₁' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'A₁' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ẋ' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ẋ' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'X1dot' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'X1dot' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' _scratch/s15/reading_build/r3_gate_new_parts.md)"
51:`python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 50 --log-directory /var/projects/toy_physics/_scratch/s15/guard/<a new directory
count: 1
```

```
$ /usr/bin/grep -n -F -- 'Log' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Log' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Laurent' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Laurent' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' _scratch/s15/reading_build/r3_gate_new_parts.md)"
40:The only place the physical symbols may be combined by hand is in writing: the supplied balance with P5's momentum density and carried current and P6's stress, as Part
count: 1
```

```
$ /usr/bin/grep -n -F -- 'Pochhammer' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pochhammer' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'hypergeom' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'hypergeom' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Appell' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Appell' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'appell' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'appell' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'quadrature' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'quadrature' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'stratum' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'stratum' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'strata' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'strata' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q=' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '900' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '900' _scratch/s15/reading_build/r3_gate_new_parts.md)"
count: 0
```

## R2. Symbols in proximity (new parts): every line naming an exponent or amplitude symbol, and every line naming d with a number

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_|b_|A_k|B_k|aa\[|bb\[)' _scratch/s15/reading_build/r3_gate_new_parts.md | cut -c1-200
```

```
$ /usr/bin/grep -n -E '`d[^a-z]|\bd\b' _scratch/s15/reading_build/r3_gate_new_parts.md | /usr/bin/grep -E '[0-9]' | cut -c1-200
18:- There are 48 objects whose item-3 closed form depends on the separation `d`: three components each of the drain-1 and mixed groups of the carried current, the centre transport and the radius tran
19:  - Each engine prints, as the condition, a finite sum of squared terms set to zero. In SymPy this is `finite_identity_condition` (`scripts/S15_separation_departure.py:293–303`, used at `:1093–
20:  - `Q_under_condition` is one line per object: the closed `Q` at `d = 2R`, carried by an Euler template (`.py:1105–1108`, `.wl:663–666`). It is the same whichever condition holds.
```

## R3. The builder stays in its lane (whole brief)

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
243:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
3:You are the new author of the S15 reading build, working in /var/projects/toy_physics. Read AGENTS.md first. The build's two engines already exist and have been reviewe
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
113:  printed them (`3192123b`), and the accepted comparator compared them (`5a019394`) (L3).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
191:Their run record is research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md. All five files are committed at 4c4e2558 as reviewed. Line
198:  - The equivalence of the finite equations to the identity is stated in prose, not computed (`.py:30–36`, `.wl:613–624`, and a printed certificate line in each t
211:- No printed line, and no run-record sentence, states a property that no computation in that engine produced: an equivalence, a bound, completeness, or that something
235:Finish with the final guarded run of each engine, from the final code. Then update research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record
count: 4
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
1:Status: the S15 reading decision list has passed its review (commit 7f52b4d2), so you are cleared to PROCEED.
3:You are the new author of the S15 reading build, working in /var/projects/toy_physics. Read AGENTS.md first. The build's two engines already exist and have been reviewe
191:Their run record is research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md. All five files are committed at 4c4e2558 as reviewed. Line
193:Review found Part B's items 1, 2, 3, 5 and 6 sound in both engines. It confirmed one defect, in item 4, in both engines.
196:  - Each engine prints, as the condition, a finite sum of squared terms set to zero. In SymPy this is `finite_identity_condition` (`scripts/S15_separation_departure.p
199:  - Independent review reached explicit conditions, and `Q` under them, for representative objects of these classes by exact routes, and checked them against a second
count: 6
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' _scratch/s15/reading_build/brief_s15_reading_repair_3.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' _scratch/s15/reading_build/brief_s15_reading_repair_3.md)"
count: 0
```

## R4. Every engine and run-record line the brief cites, at the preserved commit 4c4e2558

```
$ git show 4c4e2558:research/pde_ledger_v3/scripts/S15_separation_departure.py | sed -n '30,36p;293,303p' | cut -c1-170
  The scalar output solution space has dimension r. Its Wronskians differ by
  constant factors. A nonzero output of zero order k makes their order at
  least k-r+1. Thus k<=r*m+g*r*(r-1)/2+r-1, bounded by the same expression
  with the full system dimension N. This includes rank drops and the zero
  output without dividing by any parameter-dependent minor.
  Consequently the first B+1 exact coefficient equations at any regular
  point are necessary AND sufficient for the output to vanish identically.
    T=IdentityTransfer(transfer_name,index,size)
    R=S.MatrixSymbol(projection_name,1,size)
    Y=S.MatrixSymbol(initial_name,size,1)
    product=OrderedMatrixProduct(T,index,0,power-1)
    scalar=MatrixElement(S.MatMul(R,product,Y,evaluate=False),0,0)
    # Positive reference powers give every summand the same dimension.
    equation=S.Eq(S.Sum(S.Abs(reference**power*scalar)**2,(power,0,bound)),0,evaluate=False)
    output('matrix_aliases',{transfer_name:dict(shape=(size,size),entry_table=label+'/transfer_entries',default_entry=0,index=index),projection_name:dict(shape=(1,size),e
    output('ordered_product_definition','OrderedMatrixProduct(T,j,l,u)=T(u) ... T(l), with the identity matrix for u<l. The index is an integer; all sums and products her
    equation=S.And(regularity,equation)
    output('condition',equation)
```

```
$ git show 4c4e2558:research/pde_ledger_v3/scripts/S15_separation_departure.py | sed -n '1093,1108p' | cut -c1-170
        if reduced.has(d):
            condition_parts=[]
            for powers,coefficient in reduced_coefficients.items():
                scalar=coefficient.xreplace(backward).xreplace(back_jets)
                scalar_tag=tag+'/finite_identity_'+''.join(map(str,powers))
                def identity_progress(step):
                    print('WORKING_ON: '+step,flush=True)
                    emit('WORKING_ON',step)
                condition,certificate=finite_identity_condition(scalar,d,identity_reference_point,scalar_tag,emit,identity_progress)
                condition_parts.append(condition)
            condition=S.And(rout>rin,rin>0,*condition_parts)
            solo(tag+'/condition',condition)
            reference_value=qclosed.subs(d,identity_reference_point)
            solo(tag+'/restriction_reference_operand',reference_value)
            restricted=identity_euler_template.subs({identity_distance:d,identity_reference:identity_reference_point,identity_value:reference_value},simultaneous=True)
            solo(tag+'/Q_under_condition',restricted)
```

```
$ git show 4c4e2558:research/pde_ledger_v3/mathematica/S15_separation_departure.wl | sed -n '1,2p;610,611p;613,624p' | cut -c1-170
SetSystemOptions["ParallelOptions" -> {"ParallelThreadNumber" -> 1, "MKLThreadNumber" -> 1}];
SetSystemOptions["InequalitySolvingOptions" -> "ThreadOr" -> False];
 condition=Inactive[Sum][Abs[rv^ord matrixEntry[Inactive[Dot][projectionName,orderedMatrixProduct[productName,idx,0,ord-1],startName],1,1]]^2,{ord,0,bound}]==0;
 solo[prefix<>"/condition",condition];And[regularity,condition]];
(* Multiplicity certificate. For H Y'=M Y and f=P Y put
 P_0=P, P_(j+1)=H P_j'+P_j M-j H' P_j. Then H^j f^(j)=P_j Y.
 With m=deg P and g=max(0,deg H-1,deg M), deg P_j<=m+j g.
 If their rational rank is r, the first r rows span all subsequent rows.
 Every minor of their wedge has degree <=r m+g r(r-1)/2. Multiplication
 by an invertible analytic fundamental matrix preserves the minimum order
 of vanishing of that wedge. Its nonzero Wronskians differ by constants.
 A scalar output with zero order k therefore has k<=r m+g r(r-1)/2+r-1.
 The printed bound uses size>=r and remains valid when any parameter causes
 a rank drop. Thus the finite coefficient equations are equivalent to the
 identity, with no division by a parameter-dependent rank minor. *)
emit["item4_kernel/finite_identity_certificate","For the printed polynomial system H Y'=M Y and scalar numerator P Y, the finite multiplicity bound is N deg(P)+max(0,deg(
```

```
$ git show 4c4e2558:research/pde_ledger_v3/mathematica/S15_separation_departure.wl | sed -n '657,666p' | cut -c1-170
 Do[result=Block[{currentObject=prefix<>"/scalar_"<>ToString[k]},Catch[finiteIdentity[coefficients[[k]],dist,identityReferencePoint,prefix<>"/scalar_"<>ToString[k]],ident
  If[Head[result]===List,solo[prefix<>"/attempt_return",result];
   emit[prefix<>"/status","NOT OBTAINED EXACTLY: differential closure and finite multiplicity elimination were attempted; the displayed unsupported primitive or algebraic
   AppendTo[exactGaps,{name,"item4",result}];Return[Null]];
  AppendTo[conditions,result],{k,Length[coefficients]}];
 solo[prefix<>"/condition",rout[t]>rin[t]>0 && And@@conditions];
 referenceValue=q/.dist->identityReferencePoint;
 solo[prefix<>"/restriction_reference_operand",referenceValue];
 restricted=identityRestrictionTemplate/.{identityD->dist,identityReference->identityReferencePoint,identityValue->referenceValue};
 solo[prefix<>"/Q_under_condition",restricted];solo[prefix<>"/Q_on_family",q];
```

```
$ git show 4c4e2558:research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md | sed -n '15,18p;30,33p;39p' | cut -c1-200
- Inventory: 336 Q sections; 336 item-3 sections; 336 item-4 statuses; 48 separation-dependent objects with finite-identity conditions. Transcript end marker present: `True`.
- Finite-identity objects missing a condition: `[]`; missing restricted Q: `[]`. Remaining separation occurrences in the printed condition definitions: `[]`.

Requested-object gaps indexed by this final engine: none. General-field family integrals remain unevaluated where item 3 expressly permits them.
- Inventory: 336 Q sections; 336 item-3 sections; 336 item-4 statuses; 48 separation-dependent objects with finite-identity conditions. Transcript end marker present: `True`.
- Finite-identity objects missing a condition: `[]`; missing restricted Q: `[]`. Remaining separation occurrences in the printed condition definitions: `[]`.

Requested-object gaps indexed by this final engine: none. General-field family integrals remain unevaluated where item 3 expressly permits them.
Each engine solves its computed geometric scalar equation and inserts its closed Q at a regular symbolic reference, with the same instantaneous coefficients, to obtain Q under each condition. The fami
```

```
$ git diff --quiet 4c4e2558 -- research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md && echo 'engines and run record unchanged since 4c4e2558'; ls -la research/pde_ledger_v3/scripts/out/S15_separation_departure.out research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | awk '{print $1, $9, $10, $11}' | cut -c1-140
engines and run record unchanged since 4c4e2558
lrwxrwxrwx research/pde_ledger_v3/mathematica/out/S15_separation_departure.out -> ../../../../.git/annex/objects/p4/GJ/MD5E-s1359889560--484
lrwxrwxrwx research/pde_ledger_v3/scripts/out/S15_separation_departure.out -> ../../../../.git/annex/objects/g5/Q6/MD5E-s1197632542--88d78c1
```

```
$ echo "item4 statuses beginning 'Finite parameter equations': py $(grep -c '/item4/status: Finite parameter equations' research/pde_ledger_v3/scripts/out/S15_separation_departure.out), wl $(grep -c '/item4/status: "Finite parameter equations' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out)"
item4 statuses beginning 'Finite parameter equations': py 48, wl 48
```

```
$ test -d _scratch/s15/reading_build_work/repair3 && echo 'repair3 work dir exists' || echo 'repair3 work dir absent'
repair3 work dir absent
```
