# Measurements: S15 reading build repair-4d note (user decision: finish with gaps), leak gate and fact check

Generator: `_scratch/s15/gen/s15_reading_repair4d_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by the
orchestrator. The note carries the user decision of 2026-11-01 ("Finish with gaps"): an exact item-4 attempt that runs
out of memory counts as failed and the object is labelled. It names no value, method, exponent or condition, and points
back at `resume_s15_reading_repair_4.md` (gated at `5878692b`), whose "What must be true" stands.

## R0. The note

```
$ sha256sum _scratch/s15/reading_build/resume_s15_reading_repair_4d.md; wc -c < _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
75fec8d983557f6ccae4e11f25bc55311de419ae7d41642fbac2b51241b83f0f  _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
2028
```

## R1. Values, exponents and methods (fixed strings, from the repair-4 gate)

```
$ /usr/bin/grep -n -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'λ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '+2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 0' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b1' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '2a' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a+b' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a +' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a +' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q=' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q =' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'w=' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'S_k' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'S_k' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ẋ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ẋ' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Log' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Log' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Laurent' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Taylor' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gauss' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Gauss' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'quadrature' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'quadrature' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pochhammer' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pochhammer' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'hypergeom' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'hypergeom' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'stratum' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'strata' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'strata' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'co-moving' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'co-moving' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'comoving' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'comoving' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'continuity' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'continuity' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'flux' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'flux' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Reduce' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Reduce' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Solve' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Solve' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'MemoryConstrained' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'MemoryConstrained' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

## R3. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
21:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
13:- the run record lists, for each engine, every object that carries such a label.
19:- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.
count: 2
```

```
$ /usr/bin/grep -n -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md)"
count: 0
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
3
```

```
$ /usr/bin/grep -c -F -- 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
2
```

```
$ /usr/bin/grep -c -F -- 'can return conditions on the exponents' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
1
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
1
```

```
$ /usr/bin/grep -c -F -- 'one at a time' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
1
```

```
$ /usr/bin/grep -c -F -- 'final guarded run' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
1
```

```
$ /usr/bin/grep -c -F -- 'one that can return conditions on the exponents' _scratch/s15/reading_build/resume_s15_reading_repair_4.md; /usr/bin/grep -c -F -- 'a construction that does not reach the conditions is a reason to change the construction' _scratch/s15/reading_build/resume_s15_reading_repair_4.md; /usr/bin/grep -c -F -- 'A construction that does not reach the conditions is a reason to change the construction' _scratch/s15/reading_build/resume_s15_reading_repair_4.md
1
0
1
```

## R5. Facts the note states

```
$ ls -la _scratch/s15/reading_build_work/repair3/round90441/repair4_resource_stop_55g.md | awk '{print $5, $9}'; /usr/bin/grep -c -F 'No more memory available' _scratch/s11c/s15/repair3_90441_rule_closure_probe_01/stdout
2152 _scratch/s15/reading_build_work/repair3/round90441/repair4_resource_stop_55g.md
1
```

```
$ python3 -c 'import json;o=json.load(open("_scratch/s11c/s15/repair3_90441_rule_closure_probe_01/outcome.json"));print(o["exitCode"],o["wallSeconds"],o["childOutcome"]);print(max(int(json.loads(l).get("memory.peak",0)) for l in open("_scratch/s11c/s15/repair3_90441_rule_closure_probe_01/resource-samples.jsonl") if l.strip()))'
1 4347.966684971936 {'exitCode': 1, 'guardReason': None, 'wallSeconds': 4347.9045590509195, 'childPid': 1319806}
55644131328
```

```
$ python3 -c 'import json;print(json.load(open("_scratch/s11c/s15/repair3_90441_rule_closure_probe_01/effective-limits.json")))' | cut -c1-300
{'memory.max': '59055800320', 'memory.swap.max': '0', 'pids.max': '32', 'group': '/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-b9bcb8929773.service', 'nice': 0, 'affinity': [7], 'threads': {'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1
```

```
$ grep -E '^(MemTotal|MemAvailable):' /proc/meminfo; cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c codex
MemTotal:       65779296 kB
MemAvailable:   59395664 kB
{}

0
0
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md; echo "diff vs 90441b89 exit $?"; git status --short -- research/pde_ledger_v3
diff vs 90441b89 exit 1
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 D research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4d_note_lookups.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```
