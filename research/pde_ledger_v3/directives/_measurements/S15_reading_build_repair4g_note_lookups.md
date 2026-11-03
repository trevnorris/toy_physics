# Measurements: S15 reading build repair-4g note (user: 6 hours only), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit
is read by the orchestrator. The run and work directories named in these commands are local scratch and disposable;
every output this note relies on is quoted here in full. The note reports that the orchestrator stopped the guarded run
`repair3_90441_wolfram_final_19` on 2026-11-03 at the user's decision, and carries that decision: an exact item-4 attempt
is stopped and labelled only when it runs out of memory or is still running 6 hours after it started; silence alone is not
a stop. It names no value, method, exponent or condition, and points back at the brief and at `resume_s15_reading_repair_4.md`
(gated at `5878692b`), `resume_s15_reading_repair_4d.md` (gated at `6f46583b`) and `resume_s15_reading_repair_4f.md` (gated at
`74fb5700`).

The user decision (AskUserQuestion, 2026-11-02, local time ~21:35). Question: "How should a quiet item-4 attempt be stopped?
The harness now labels an object 'not obtained' after 5 minutes with no output. Object 4 got its label after 16 minutes at
3.6 GB. A kernel that is quiet at 100% CPU is usually in one long solve step." Answer chosen: **"6 hours only (Recommended)"**,
whose option text read: "Matches what you chose: an object is stopped only when it runs out of memory or reaches 6 hours.
I send a short note: drop the 5-minute rule, retry object 4, and rerun the final from the changed harness (objects 1–3 are
instant labels, so that costs little). This gives fewer gaps but may take longer: up to 6 h for each hard object, and how many
are hard is unknown. I'll keep counting and ask again if the estimate gets large."

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md; wc -c < research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
295285875b23599c40afe586ef473a4e63fc10b7302a6e8d7aa7c6853a6e87c3  research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
2731
```

## R1. Values, exponents and methods (fixed strings, from the repair-4 gate)

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gamma' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gamma' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'γ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'γ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_19` through its unit (`systemctl --user stop s11c-guard-4f8825ff30c7`), with its launcher alive; the journal record
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a +' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'a +' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'S_k' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'S_k' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ṙ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ṙ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Rdot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Rdot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Ẋ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Ẋ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Log' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Laurent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Laurent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'expansion' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'expansion' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Taylor' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Taylor' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Gauss' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Gauss' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'quadrature' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'quadrature' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pochhammer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Pochhammer' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'hypergeom' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'hypergeom' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'stratum' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'stratum' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'strata' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'strata' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'co-moving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'co-moving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'comoving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'comoving' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'continuity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'continuity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'flux' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'flux' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Reduce' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Reduce' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Solve' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'Solve' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'MemoryConstrained' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'MemoryConstrained' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'matrix' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'matrix' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pivot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'pivot' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'eliminate' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -F -- 'eliminate' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

## R3. The builder stays in its lane

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
25:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_19` through its unit (`systemctl --user stop s11c-guard-4f8825ff30c7`), with its launcher alive; the journal record
17:- every item-4 object that carries a label carries one of the two kinds, the 6-hour stop or running out of memory, and the run record says which.
23:- Finish as the brief says: the final guarded run of each engine from the final code, the updated run record, and a short final report.
count: 3
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md)"
1:Note from the orchestrator. I stopped `repair3_90441_wolfram_final_19`, and I have read your report of that stop (`repair4_external_stop_final19.md`). The user has made
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_19` through its unit (`systemctl --user stop s11c-guard-4f8825ff30c7`), with its launcher alive; the journal record
11:So `region/component_2/rate_labels_{drain1, drain2}` carries no label from `repair3_90441_wolfram_final_19`; attempt it again. The two labels for running out of memory
count: 3
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
2
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4f.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
3
```

```
$ /usr/bin/grep -c -F -- 'brief_s15_reading_repair_3.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- 'silent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
5
```

```
$ /usr/bin/grep -c -F -- 'attempt it again' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- 'line 221' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
2
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- 'one at a time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- '6 hours' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
1
```

```
$ /usr/bin/grep -c -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
0
```

```
$ sed -n '221p' _scratch/s15/reading_build/brief_s15_reading_repair_3.md
Every engine run, probe or final, prints and flushes a progress line before each integration or solve it starts, naming the object and the step. A run that stops printing for a long stretch has failed: find from its last lines what it was doing and change the construction; do not wait on a silent run. A run you stop yourself is part of the repair: change the construction and continue in this session.
```

```
$ /usr/bin/grep -n -F 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:18:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:7:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:9:Everything under "What must be true" in `resume_s15_reading_repair_4.md` stands. For an object whose attempt ran out of memory, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md` and in `resume_s15_reading_repair_4d.md` stands. For an object whose attempt was stopped at 6 hours, these are also true:
```

```
$ /usr/bin/grep -n -F -e 'silent' -e '6 hours' research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md | cut -c1-260
9:The user chose to label and finish. An exact item-4 attempt on an object that is still running 6 hours after it started counts as an attempt that failed, as one that runs out of memory does under `resume_s15_reading_repair_4d.md`. That object then carries a 
15:Everything under "What must be true" in `resume_s15_reading_repair_4.md` and in `resume_s15_reading_repair_4d.md` stands. For an object whose attempt was stopped at 6 hours, these are also true:
17:- its label says what was attempted, in which run, that the attempt was stopped after 6 hours or more without returning, the last step it reported, how long it ran, and the peak memory of the run;
24:- A stop at 6 hours is not a reason to stop and report; it is the case above.
25:- For an item-4 object, an attempt that goes silent for a long stretch (`brief_s15_reading_repair_3.md` line 221) is handled the same way: label it, with the last step it reported, instead of changing the construction.
```

## R5. Facts the note states

S0, captured by the orchestrator at 2026-11-03T04:40:15Z, before the stop (commands and literal output):

````
snapshot_utc: 2026-11-03T04:40:15Z
unit: s11c-guard-4f8825ff30c7 active=active
cgroup: /user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-4f8825ff30c7.service
anon 11732848640
inactive_file 2281508864
memory.peak: 14142205952
stdout: 12651226 2026-11-02 21:40:15.651771130 -0700
WORKING_ON count: 85628
last: WORKING_ON: reading/region/component_3/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / branch_683 Solve exponent stratum equalities
execution.json sha256: 5dfa41da5f1c51e574b14e1949f6a4133df22384d1107963a9695d7a976991e8
kernels: 7 order: 48
  base  2026-11-03T03:32:52.229425+00:00 2026-11-03T03:36:10.335880+00:00 
  item4 region/component_1/momentum_labels_{drain1, drain2} 2026-11-03T03:36:10.499828+00:00 2026-11-03T03:36:25.504233+00:00 
  item4 region/component_1/rate_labels_{drain1, drain2} 2026-11-03T03:36:25.509000+00:00 2026-11-03T03:37:35.520500+00:00 
  item4 region/component_2/momentum_labels_{drain1, drain2} 2026-11-03T03:37:35.525139+00:00 2026-11-03T03:37:51.529061+00:00 
  item4 region/component_2/rate_labels_{drain1, drain2} 2026-11-03T03:37:51.532980+00:00 2026-11-03T03:55:13.847045+00:00 silent_stop
  item4 region/component_2/rate_labels_{drain1, drain2} 2026-11-03T03:55:13.848374+00:00 2026-11-03T03:56:20.856206+00:00 
  item4 region/component_3/momentum_labels_{drain1, drain2} 2026-11-03T03:56:20.859861+00:00 None 
failed_objects kinds: {'region/component_1/momentum_labels_{drain1, drain2}': 'memory_exhausted', 'region/component_1/rate_labels_{drain1, drain2}': 'memory_exhausted', 'region/component_2/momentum_labels_{drain1, drain2}': 'six_hour_stop', 'region/component_2/rate_labels_{drain1, drain2}': 'silent_stop'}
run.py sha256: 4d4fd2cb0c5c73d2f7485456d298b50d4e06bee4a7bd7bffb7ccd94c69dd9279
66:ATTEMPT_SECONDS = 6 * 60 * 60
67:SILENT_SECONDS = 5 * 60
codex: 2953592  0.2    01:09:48 codex
log: 403723 2026-11-02 21:39:31.281510405 -0700
````

```
$ journalctl --user -u s11c-guard-4f8825ff30c7.service --no-pager -o short-iso-precise --since '2026-11-03 04:40:00 UTC' --until '2026-11-03 04:41:00 UTC' --utc | cut -c1-200
2026-11-03T04:40:21.706420+00:00 pop-os systemd[1646]: Stopping s11c-guard-4f8825ff30c7.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-03T04:40:26.808055+00:00 pop-os systemd[1646]: s11c-guard-4f8825ff30c7.service: State 'final-sigterm' timed out. Killing.
2026-11-03T04:40:26.808131+00:00 pop-os systemd[1646]: s11c-guard-4f8825ff30c7.service: Killing process 2962415 (WolframKernel) with signal SIGKILL.
2026-11-03T04:40:26.808226+00:00 pop-os systemd[1646]: s11c-guard-4f8825ff30c7.service: Failed with result 'timeout'.
2026-11-03T04:40:26.808254+00:00 pop-os systemd[1646]: s11c-guard-4f8825ff30c7.service: Unit process 2962415 (WolframKernel) remains running after unit stopped.
2026-11-03T04:40:26.808496+00:00 pop-os systemd[1646]: Stopped s11c-guard-4f8825ff30c7.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-03T04:40:26.808654+00:00 pop-os systemd[1646]: s11c-guard-4f8825ff30c7.service: Consumed 1h 7min 19.456s CPU time, 13.1G memory peak, 0B memory swap peak.
```

```
$ cat _scratch/s11c/s15/repair3_90441_wolfram_final_19/outcome.json 2>&1; echo; /usr/bin/grep -c -F 'No more memory available' _scratch/s11c/s15/repair3_90441_wolfram_final_19/stdout
{
  "exitCode": 1,
  "wallSeconds": 4054.6722250099992,
  "unit": "s11c-guard-4f8825ff30c7",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}

0
```

```
$ /usr/bin/grep -n -E '^(ATTEMPT|SILENT)_SECONDS' research/pde_ledger_v3/mathematica/S15_separation_departure_run.py; sha256sum research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
66:ATTEMPT_SECONDS = 6 * 60 * 60
67:SILENT_SECONDS = 5 * 60
4d4fd2cb0c5c73d2f7485456d298b50d4e06bee4a7bd7bffb7ccd94c69dd9279  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ python3 -c "import json;s=json.load(open('_scratch/s15/reading_build_work/repair3/round90441/wolfram_final19_execution/execution.json'));o=s['item4_worker_order'];print(len(o), o.index('region/component_3/momentum_labels_{drain1, drain2}')+1);k=[k for k in s['kernels'] if k.get('author_stop_kind')];print([(x['object'],x['author_stop_kind'],x['attempt_seconds_at_stop'],x['guarded_run_peak_at_stop']) for x in k]);print({n:v['kind'] for n,v in s['failed_objects'].items()})"
48 5
[('region/component_2/rate_labels_{drain1, drain2}', 'silent_stop', 975.9356204340002, 3583213568)]
{'region/component_1/momentum_labels_{drain1, drain2}': 'memory_exhausted', 'region/component_1/rate_labels_{drain1, drain2}': 'memory_exhausted', 'region/component_2/momentum_labels_{drain1, drain2}': 'six_hour_stop', 'region/component_2/rate_labels_{drain1, drain2}': 'silent_stop'}
```

```
$ ls -l _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final19.md | awk '{print $5, $9}'; sha256sum _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final19.md; /usr/bin/grep -n -F -e 'recorded peak' -e 'WORKING_ON' -e 'silent-stop' -e 'do not attempt' _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final19.md | cut -c1-220
2560 _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final19.md
805f0add4b2c4a66719aec79c61f668c668fd55e222145585c11407da84cd15a  _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final19.md
9:Guard exit **1**, child outcome `null`, wall time **4054.6722250099992 s**, recorded peak **14142205952 bytes**. The last sample reports `max 0`, `oom 0`, `oom_kill 0`, `oom_group_kill 0`. The journal records an extern
14:WORKING_ON: reading/region/component_3/momentum_labels_{drain1, drain2} / item4 exponent-sensitive / branch_690 PiecewiseExpand Q and departure strata
17:The preceding component-2 mixed rate attempt received a **silent-stop** label: exact-attempt elapsed `975.9356204340002` s; last step `rational direction reconstruction residual`; run peak at failure `3583213568` byte
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c codex; cat _scratch/s15/reading_build/repair4f.DONE
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
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4g_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```
