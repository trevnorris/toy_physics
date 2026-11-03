# Measurements: S15 reading build repair-4f note (user: label and finish, 6-hour stop point), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit
is read by the orchestrator. The run and work directories named in these commands are local scratch and disposable;
every output this note relies on is quoted here in full. The note reports that the orchestrator stopped the guarded run
`repair3_90441_wolfram_final_18` on 2026-11-03 at the user's decision, and carries that decision: an exact item-4 attempt
still running 6 hours after it started counts as failed and is labelled. It names no value, method, exponent or condition,
and points back at the brief and at `resume_s15_reading_repair_4.md` (gated at `5878692b`), `resume_s15_reading_repair_4d.md`
(gated at `6f46583b`) and `resume_s15_reading_repair_4e.md` (gated at `eceafc57`).

The user decision (AskUserQuestion, 2026-11-02, after the user asked "This feels like a lot of time for something that might
fail. Are you certain we're taking the best course of action to solve this problem?"). Question: "How should we handle the
item-4 objects that won't finish?" Answer chosen: **"Label and finish (Recommended)"**, whose option text read: "Stop final_18
now and label object 3 honestly: what was tried, how long it ran, the last step, peak memory. Any later exact attempt still
running after 6 hours gets the same kind of label. Finish the other objects, run SymPy, then review. The gaps go on the
more-RAM list for a cheaper numeric check later. Roughly a day or two to review."

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md; wc -c < research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
7563293f8774de67fa854ce3daadb13627fa7b5471381033f2de8ef7308e9cfd  research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
3170
```

## R1. Values, exponents and methods (fixed strings, from the repair-4 gate)

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
5:At 2026-11-03T03:26:17Z I stopped `repair3_90441_wolfram_final_18` through its unit (`systemctl --user stop s11c-guard-b084eb25cac6`), with its launcher alive. The run 
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a +' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a +' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'S_k' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'S_k' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ẋ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ẋ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Laurent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Laurent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Taylor' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Taylor' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gauss' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Gauss' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'quadrature' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'quadrature' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pochhammer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pochhammer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'hypergeom' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'hypergeom' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'stratum' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'stratum' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'strata' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'strata' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'co-moving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'co-moving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'comoving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'comoving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'continuity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'continuity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'flux' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'flux' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Reduce' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Reduce' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Solve' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Solve' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'MemoryConstrained' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'MemoryConstrained' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'matrix' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'matrix' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pivot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pivot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'eliminate' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'eliminate' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

## R3. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
29:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
5:At 2026-11-03T03:26:17Z I stopped `repair3_90441_wolfram_final_18` through its unit (`systemctl --user stop s11c-guard-b084eb25cac6`), with its launcher alive. The run 
20:- the run record lists, for each engine, every object that carries such a label or the label for running out of memory, and which of the two it is.
27:- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.
count: 3
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md)"
1:Note from the orchestrator. I stopped `repair3_90441_wolfram_final_18`, and I have read your report of that stop (`repair4_external_stop_final18.md`). The user has made
5:At 2026-11-03T03:26:17Z I stopped `repair3_90441_wolfram_final_18` through its unit (`systemctl --user stop s11c-guard-b084eb25cac6`), with its launcher alive. The run 
11:The attempt on `region/component_2/momentum_labels_{drain1, drain2}` in `repair3_90441_wolfram_final_18` is such an attempt. Label that object from it, and do not atte
count: 3
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
2
```

```
$ /usr/bin/grep -c -F -- 'brief_s15_reading_repair_3.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
1
```

```
$ /usr/bin/grep -c -F -- 'line 221' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
2
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
1
```

```
$ /usr/bin/grep -c -F -- 'one at a time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
1
```

```
$ /usr/bin/grep -c -F -- '6 hours' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
4
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
1
```

```
$ /usr/bin/grep -c -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
0
```

```
$ sed -n '221p' _scratch/s15/reading_build/brief_s15_reading_repair_3.md
Every engine run, probe or final, prints and flushes a progress line before each integration or solve it starts, naming the object and the step. A run that stops printing for a long stretch has failed: find from its last lines what it was doing and change the construction; do not wait on a silent run. A run you stop yourself is part of the repair: change the construction and continue in this session.
```

```
$ /usr/bin/grep -n -F 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:18:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:7:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:9:Everything under "What must be true" in `resume_s15_reading_repair_4.md` stands. For an object whose attempt ran out of memory, these are also true:
```

## R5. Facts the note states

S0, captured by the orchestrator at 2026-11-03T03:26:08Z, before the stop (commands and literal output):

````
```
$ date -u '+%Y-%m-%dT%H:%M:%SZ'
2026-11-03T03:26:08Z
```

```
$ systemctl --user show -p ActiveState -p SubState -p ControlGroup s11c-guard-b084eb25cac6.service
ControlGroup=/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-b084eb25cac6.service
ActiveState=active
SubState=running
```

```
$ stat -c '%s %y' _scratch/s11c/s15/repair3_90441_wolfram_final_18/stdout
17724789 2026-11-02 20:26:07.938873209 -0700
```

```
$ LC_ALL=C grep -a -c WORKING_ON _scratch/s11c/s15/repair3_90441_wolfram_final_18/stdout
122883
```

```
$ tail -c 3000000 _scratch/s11c/s15/repair3_90441_wolfram_final_18/stdout | LC_ALL=C grep -a WORKING_ON | tail -1
WORKING_ON: reading/region/component_2/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / SeriesCoefficient 99 column 130
```

```
$ cat /sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-b084eb25cac6.service/memory.peak /sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-b084eb25cac6.service/memory.current
42904752128
42901696512
```

```
$ grep -E '^(anon|inactive_file) ' /sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-b084eb25cac6.service/memory.stat
anon 42788458496
inactive_file 26865664
```

```
$ cat /sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-b084eb25cac6.service/memory.events
low 0
high 0
max 0
oom 0
oom_kill 0
oom_group_kill 0
sock_throttled 0
```

```
$ python3 -c "import json;s=json.load(open('_scratch/s15/reading_build_work/repair3/round90441/wolfram_final18_execution/execution.json'));print(len(s['item4_worker_order']), s['item4_worker_order'].index('region/component_2/momentum_labels_{drain1, drain2}')+1); [print(k['mode'],k.get('object'),k.get('started'),k.get('ended'),k.get('pid')) for k in s['kernels']]"
48 3
base  2026-11-02T17:29:44.285187+00:00 2026-11-02T17:32:56.769773+00:00 2514963
item4 region/component_1/momentum_labels_{drain1, drain2} 2026-11-02T17:32:56.925543+00:00 2026-11-02T17:33:11.681130+00:00 2518113
item4 region/component_1/rate_labels_{drain1, drain2} 2026-11-02T17:33:11.682407+00:00 2026-11-02T17:34:20.169959+00:00 2518436
item4 region/component_2/momentum_labels_{drain1, drain2} 2026-11-02T17:34:20.170696+00:00 None 2518961
```

```
$ ls -l _scratch/s15/reading_build_work/repair3/round90441/wolfram_final18_execution
total 17332
-rw-rw-r-- 1 trevnorris trevnorris    12777 Nov  2 10:34 execution.json
-rw-rw-r-- 1 trevnorris trevnorris   304963 Nov  2 10:32 kernel_1.stdout
-rw-rw-r-- 1 trevnorris trevnorris    33971 Nov  2 10:33 kernel_2.stdout
-rw-rw-r-- 1 trevnorris trevnorris    46663 Nov  2 10:34 kernel_3.stdout
-rw-rw-r-- 1 trevnorris trevnorris 17334667 Nov  2 20:26 kernel_4.stdout
```

```
$ ps -eo pid,ppid,rss,pcpu,etime,comm --sort=-rss | head -5
    PID    PPID   RSS %CPU     ELAPSED COMMAND
2518961 2514962 41866052 99.9 09:51:48 WolframKernel
 175046    3478 376856 0.8  1-09:29:57 claude
2223810 2223809 231544 0.1    20:14:39 codex
   1923    1646 99484  0.1 10-22:54:54 gnome-shell
```

````

```
$ journalctl --user -u s11c-guard-b084eb25cac6.service --no-pager -o short-iso-precise --since '2026-11-03 03:26:00 UTC' --until '2026-11-03 03:27:00 UTC' --utc | cut -c1-200
2026-11-03T03:26:17.240998+00:00 pop-os systemd[1646]: Stopping s11c-guard-b084eb25cac6.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-03T03:26:22.308078+00:00 pop-os systemd[1646]: s11c-guard-b084eb25cac6.service: State 'final-sigterm' timed out. Killing.
2026-11-03T03:26:22.308145+00:00 pop-os systemd[1646]: s11c-guard-b084eb25cac6.service: Killing process 2518961 (WolframKernel) with signal SIGKILL.
2026-11-03T03:26:22.308240+00:00 pop-os systemd[1646]: s11c-guard-b084eb25cac6.service: Failed with result 'timeout'.
2026-11-03T03:26:22.308269+00:00 pop-os systemd[1646]: s11c-guard-b084eb25cac6.service: Unit process 2518961 (WolframKernel) remains running after unit stopped.
2026-11-03T03:26:22.308504+00:00 pop-os systemd[1646]: Stopped s11c-guard-b084eb25cac6.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-03T03:26:22.308672+00:00 pop-os systemd[1646]: s11c-guard-b084eb25cac6.service: Consumed 9h 56min 29.481s CPU time, 39.9G memory peak, 0B memory swap peak.
```

```
$ cat _scratch/s11c/s15/repair3_90441_wolfram_final_18/outcome.json; echo; /usr/bin/grep -c -F 'No more memory available' _scratch/s11c/s15/repair3_90441_wolfram_final_18/stdout
{
  "exitCode": 1,
  "wallSeconds": 35798.11544563505,
  "unit": "s11c-guard-b084eb25cac6",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}

0
```

```
$ ls -l _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final18.md | awk '{print $5, $9}'; /usr/bin/grep -n -F -e 'recorded peak' -e 'WORKING_ON' _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final18.md | cut -c1-200
3138 _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final18.md
13:Guard exit `1`, child outcome `null`, wall time `35798.11544563505` s, recorded peak `42911830016` bytes. The last sample has `max 0`, `oom 0`, `oom_kill 0`, and `oom_group_kill 0`. This is not an 
18:WORKING_ON: reading/region/component_2/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / eliminate series pivot 1
```

```
$ for r in final_16 final_17 final18; do f=_scratch/s15/reading_build_work/repair3/round90441/wolfram_${r}_execution/execution.json; python3 -c "import json,sys;s=json.load(open(sys.argv[1]));o=s['item4_worker_order'];ks=[k for k in s['kernels'] if k['mode']=='item4'];print(sys.argv[1].split('/')[-2], len(o), o.index('region/component_2/momentum_labels_{drain1, drain2}')+1, [(k['object'],k.get('ended')) for k in ks if k['object']=='region/component_2/momentum_labels_{drain1, drain2}'])" $f; done
wolfram_final_16_execution 48 3 [('region/component_2/momentum_labels_{drain1, drain2}', None)]
wolfram_final_17_execution 48 3 [('region/component_2/momentum_labels_{drain1, drain2}', None)]
wolfram_final18_execution 48 3 [('region/component_2/momentum_labels_{drain1, drain2}', None)]
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c codex; cat _scratch/s15/reading_build/repair4e.DONE
{}

0
0
exit 0
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md; echo "run record diff vs 90441b89 exit $?"; git status --short -- research/pde_ledger_v3
run record diff vs 90441b89 exit 0
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4f_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```
