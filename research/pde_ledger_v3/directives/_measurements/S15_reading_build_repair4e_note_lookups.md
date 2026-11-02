# Measurements: S15 reading build repair-4e note (orchestrator stop of a silent run), leak gate and fact check

Generator: `_scratch/s15/gen/s15_reading_repair4e_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by the
orchestrator. The note reports that the orchestrator stopped the silent guarded run `repair3_90441_wolfram_final_17` on
2026-11-02. It names no value, method, exponent or condition, and points back at the brief and at
`resume_s15_reading_repair_4.md` (gated at `5878692b`) and `resume_s15_reading_repair_4d.md` (gated at `6f46583b`).

## R0. The note

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_4e.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
e9e829e2a8605410c147aed9aafd4de727b00536c3c4597bd1ca90bb5f916adc  _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1409
```

## R1. Values, exponents and methods (fixed strings, from the repair-4 gate)

```
$ /usr/bin/grep -n -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'λ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
5:At 2026-11-02T07:08:26Z I stopped `repair3_90441_wolfram_final_17` through its unit (`systemctl --user stop s11c-guard-4d6e832b5128`), with its launcher alive. Its stdo
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '+2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 0' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b1' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '2a' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a+b' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a +' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a +' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q=' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q =' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'w=' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'S_k' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'S_k' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ẋ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ẋ' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Log' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Log' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gauss' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Gauss' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'quadrature' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'quadrature' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pochhammer' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pochhammer' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'hypergeom' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'hypergeom' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'strata' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'strata' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'co-moving' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'co-moving' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'comoving' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'comoving' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'continuity' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'continuity' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'flux' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'flux' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Reduce' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Reduce' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Solve' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Solve' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'MemoryConstrained' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'MemoryConstrained' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

## R3. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
17:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
5:At 2026-11-02T07:08:26Z I stopped `repair3_90441_wolfram_final_17` through its unit (`systemctl --user stop s11c-guard-4d6e832b5128`), with its launcher alive. Its stdo
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md)"
count: 0
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
2
```

```
$ /usr/bin/grep -c -F -- 'brief_s15_reading_repair_3.md' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ /usr/bin/grep -c -F -- 'line 221' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ /usr/bin/grep -c -F -- 'line 36' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ /usr/bin/grep -c -F -- 'one job at a time' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' _scratch/s15/reading_build/resume_s15_reading_repair_4e.md
1
```

```
$ sed -n '221p' _scratch/s15/reading_build/brief_s15_reading_repair_3.md
Every engine run, probe or final, prints and flushes a progress line before each integration or solve it starts, naming the object and the step. A run that stops printing for a long stretch has failed: find from its last lines what it was doing and change the construction; do not wait on a silent run. A run you stop yourself is part of the repair: change the construction and continue in this session.
```

```
$ sed -n '36p' _scratch/s15/reading_build/resume_s15_reading_repair_4.md
Wait for each guarded job by polling it until it ends: check its process and its outcome.json at short intervals. Do not rely on queued notices or detached watchers, and do not use long sleeps. Keep each run's progress lines flowing as the brief says; a run that goes silent for a long stretch has failed, and you may stop your own run through its live guard as the brief describes.
```

## R5. Facts the note states

```
$ cat _scratch/s15/reading_build/final17_stop_snapshot.txt
2026-11-02T07:08:14Z
unit s11c-guard-4d6e832b5128 launcher 1829009
stdout 9939547 bytes, last write 2026-11-01 22:48:44.624552394 -0700
memory.current 40567287808 memory.peak 40569319424
WORKING_ON: reading/region/component_2/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / linear_parameter_branch_108 Reduce pivot exponent and radius stratum
    PID   RSS %CPU     ELAPSED     TIME COMMAND
1833312 38663252 99.8 01:57:23 01:57:14 WolframKernel
stop at 2026-11-02T07:08:26Z
```

```
$ ls -la _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final17.md | awk '{print $5, $9}'; stat -c '%s %y' _scratch/s11c/s15/repair3_90441_wolfram_final_17/stdout; /usr/bin/grep -c -F 'No more memory available' _scratch/s11c/s15/repair3_90441_wolfram_final_17/stdout
3234 _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final17.md
9939547 2026-11-01 22:48:44.624552394 -0700
0
```

```
$ python3 -c 'import json;o=json.load(open("_scratch/s11c/s15/repair3_90441_wolfram_final_17/outcome.json"));print(o["exitCode"],o["wallSeconds"],o["childOutcome"],o["unit"]);rows=[json.loads(l) for l in open("_scratch/s11c/s15/repair3_90441_wolfram_final_17/resource-samples.jsonl") if l.strip()];print(max(int(r.get("memory.peak",0)) for r in rows));print(repr(rows[-1]["memory.events"]))'
1 7331.1611233869335 None s11c-guard-4d6e832b5128
40569319424
'low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\noom_group_kill 0\nsock_throttled 0'
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c codex
{}

0
0
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md; echo "diff vs 90441b89 exit $?"; git status --short -- research/pde_ledger_v3
diff vs 90441b89 exit 1
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4e_note_lookups.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```
