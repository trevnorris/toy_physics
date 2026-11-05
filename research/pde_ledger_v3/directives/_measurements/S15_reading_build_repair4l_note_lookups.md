# Measurements: S15 reading build repair-4l note (user: prototype the expansion route), leak gate, fact check and review

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes, process and
unit state. Lines cut to a byte width can end inside a character; the record drops such partial characters. Every hit is read by the orchestrator. The run, work and review directories named in these commands are local
scratch and disposable; every output this record relies on is quoted here in full. The note supplies an exact expansion of
the distance from drain 1 over the region (supplied, unfalsifiable in the trial) and asks for a one-object trial of an
item-4 construction that uses it, in each engine, under the 1-hour trial stop of `resume_s15_reading_repair_4i.md`
(gated at `4d2c29fc`). It names no value, exponent condition or expected outcome. Because both engines read the supplied
equations, the note is physics-bearing and was reviewed by two legs (orchestrator-written: Codex and Grok); see R6.

## The user's question and decision (2026-11-04 / 2026-11-05)

After the guard refused the 55 GiB Wolfram rerun of `resume_s15_reading_repair_4k.md` (gated at `04fa23f6`), the user
asked, in their own words: "I do genuinely wonder how an actual computational physicist would do this work. Feels like
we're doing something inefficient for this to take this long. Are we using the wrong tech to answer this question? I'd
like your thoughts." The orchestrator answered that the cost comes mainly from how item 4 is constructed, and that the
region lies wholly inside `r ≤ R < d`, so the distance from drain 1 has an exact series in `r/d` there.

AskUserQuestion. Question: "The 55 GiB Wolfram rerun was refused (only about 52 GiB fits beside SymPy). How should we
proceed with item 4?" Options offered:
- "Prototype expansion route (Recommended)": "SymPy keeps running in the background. A note supplying the series
  expansion of drain 1's flow in r/d gets two review legs, then Codex tries that route on one object in both engines,
  timed. If it's fast, the decision list changes (two legs) and all 48 objects run that way."
- "Wolfram at 52 GiB now": "Continue the current measurement: Wolfram reruns now beside SymPy at the most the guard
  allows (about 52 GiB), and a stop there counts as Wolfram's memory limit for this measurement."
- "Wait, then 55 GiB alone": "The plan as written: Wolfram reruns at 55 GiB only after both SymPy runs finish, which
  could take days. Nothing else changes meanwhile."
- "Narrow item 4 first": "Decide which of the 48 objects the record needs, for example only the totals, then change the
  decision list with two review legs before any more runs. SymPy keeps running meanwhile."

The user selected **"Prototype expansion route (Recommended)"**. The orchestrator carried the 1-hour trial stop of
`resume_s15_reading_repair_4i.md` into the trial and told the user so.

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
4d46c703ac54bc5d1f434c3ef2f822d40e70c9feffc9c8ae3d9d45ff82dc1ddd  research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
4202 research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
16:With the region's surfaces taken as spheres about `X_2` of radii `r_in` and `R` (`S15_decision_list.md:89–91`), write a point of the region as `x = X_2 + r e`, with `e` a unit vector and `r_in 
21:  `g_n(s, μ) = Σ_{k=0}^{⌊n/2⌋} (−1)^k (−s/2)_{n−k} (2μ)^{n−2k} / (k! (n−2k)!)`,
23:  with `(c)_0 = 1` and `(c)_m = c(c+1)⋯(c+m−1)`. The series converges absolutely and uniformly for `0 ≤ z ≤ z_0` and `−1 ≤ μ ≤ 1`, for every `z_0 < 1`.
24:- For `−s/2 > −1/2` and `s ≠ 0`, `g_n(s, μ)` is the Gegenbauer polynomial `C_n^{(−s/2)}(μ)` (DLMF 18.5.10 and 18.12.4). Here `g_n` is defined by the formula above for every real `s`.
count: 4
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
21:  `g_n(s, μ) = Σ_{k=0}^{⌊n/2⌋} (−1)^k (−s/2)_{n−k} (2μ)^{n−2k} / (k! (n−2k)!)`,
count: 1
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
16:With the region's surfaces taken as spheres about `X_2` of radii `r_in` and `R` (`S15_decision_list.md:89–91`), write a point of the region as `x = X_2 + r e`, with `e` a unit vector and `r_in 
count: 1
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
38:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md)"
11:- The Wolfram attempts that `resume_s15_reading_repair_4j.md` and `resume_s15_reading_repair_4k.md` ask for are not made.
35:- Write the construction in each engine, in your work directory, and attempt `region/component_1/momentum_labels_{drain1, drain2}` with it, Wolfram first. Run through the guard, one job at a time, 
count: 2
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4f.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4g.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4h.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4i.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4j.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4k.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
2
```

```
$ /usr/bin/grep -c -F -- 'S15_reading_decision_list.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'S15_decision_list.md:89–91' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'S15_decision_list.md:120' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'Supplied (unfalsifiable in this trial)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
2
```

```
$ /usr/bin/grep -c -F -- 'No limit or expansion is taken in place of an exact object' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'Each attempt stops 1 hour after it started' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not attempt any other object' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not start the final runs' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
0
```

```
$ /usr/bin/grep -c -F -- 'No production source changes' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
1
```

```
$ ls _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4{f,g,h,i,j,k}.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4j.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4k.md
```

```
$ /usr/bin/grep -n -F 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4{f,g,h,i}.md
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
$ /usr/bin/grep -n -F -e '1 hour after it started' -e 'No trial outcome is printed as an item-4 label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4i.md | cut -c1-220
9:The trial measures how fast your construction is. The user decided that in this trial, an attempt is stopped 1 hour after it started and is reported as not returned within 1 hour. The Wolfram trial of `region/component
16:- No trial outcome is printed as an item-4 label. A trial attempt ends when it returns its conditions, when it runs out of memory, or 1 hour after it started.
23:- Run through the guard, one job at a time, with `--memory-gib 55`, and poll each job until it ends. Each attempt stops 1 hour after it started. If an attempt ends on an error in your code, you may fix it and run that
```

```
$ sed -n '32,33p' _scratch/s15/reading_build/resume_s15_reading_repair_4.md
Stop and report only if:
- the guard refuses admission or ends a run you did not stop; or
```

## R5. Facts the note states

The refusal, and the S15 region and domain the note cites.

```
$ tail -1 _scratch/s15/reading_build_work/repair3/size_measurement/repair4_size_wl_c1_momentum_02_launcher.log; test -e _scratch/s11c/s15/repair4_size_wl_c1_momentum_02/outcome.json && echo 'outcome.json present' || echo 'outcome.json absent'; cat _scratch/s15/reading_build/repair4k.DONE
RuntimeError: available memory does not cover outstanding commitments plus this job
outcome.json absent
exit 0
```

```
$ sha256sum _scratch/s15/reading_build_work/repair3/size_measurement/measurement_guard_refusal_02.md; /usr/bin/grep -n -F -e 'refused at admission' -e 'No Wolfram kernel started' -e 'neither stopped nor restarted' _scratch/s15/reading_build_work/repair3/size_measurement/measurement_guard_refusal_02.md | cut -c1-220
2d98dac3432e81f34bf1e3321eeaaf7732e7cb12dd898bf6aca803d4b431a876  _scratch/s15/reading_build_work/repair3/size_measurement/measurement_guard_refusal_02.md
3:The requested 55 GiB Wolfram momentum retry was refused at admission. No Wolfram kernel started; the rate retry was not launched. The two original SymPy attempts remain running and were neither stopped nor restarted. T
```

```
$ sed -n '89,91p;120p' research/pde_ledger_v3/directives/S15_decision_list.md
- **The region** lies around drain 2, between an inner surface around its core and an outer surface, both moving
  with it. The build may take them as spheres about `X_2(t)` with live radii `r_in(t)` and `R(t)`. That is a choice
  of control surface, labelled where printed. Drain 1 lies outside the outer surface.
Every object is stated on its domain: `|X_1 − X_2| > R > r_in > 0`, with both cores excluded, and each pressure
```

```
$ sed -n '36p;78,82p;86,88p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
- **The separation:** `D = X_1 − X_2`, `d = |D|`. The domain requires `d > R` (L2).
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
   - `Q` on the family under each condition.
Every object is exact on S15's domain under the conditions it carries: S15's pressure conditions, and for each
departure object a condition under which it exists, as above. No limit or expansion is taken in place of an exact
object. Where an engine does not obtain an object exactly, it prints a label saying so, beside what it did print.
```

```
$ sed -n '445,449p' scripts/s11c_guarded_run.py; sed -n '284,286p' scripts/s11c_guarded_run.py
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-directory', type=Path)
    parser.add_argument('--seconds', type=int, default=0,
                        help='legacy compatibility argument; recorded but never imposes a deadline')
        raise RuntimeError('no unreserved physical CPU available')
    if available < AVAILABLE_MIN + headroom + requested:
        raise RuntimeError('available memory does not cover outstanding commitments plus this job')
```

The two SymPy measurement attempts.

```
$ systemctl --user is-active s11c-guard-c4214c0da91a.service s11c-guard-fd3dd4540a3c.service; python3 -c "import json;r=json.load(open('/home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json'));[print(k,v['memoryMax'],v['owner']['pid'],v['logDirectory'].rsplit('/',1)[1]) for k,v in r.items()]"
inactive
active
s11c-guard-fd3dd4540a3c.service 2147483648 237635 repair4_size_py_c1_rate_01
```

```
$ echo "WolframKernel processes: $(pgrep -c WolframKernel)"; echo "processes named codex: $(pgrep -c -x codex)"
WolframKernel processes: 0
processes named codex: 0
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
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4l_note_lookups.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ echo "final run dirs after final_20: $(ls -d _scratch/s11c/s15/*final_2[1-9]* 2>/dev/null | wc -l)"
final run dirs after final_20: 0
```

## R6. Review: two legs, one identical prompt, both reports

The note was reviewed uncommitted at the baseline below. Both legs received the same rendered prompt (quoted in full at
the end of this section). The report hashes were recorded before either report was read.

```
$ cat _scratch/s15/reading_note_4l/note_4l_baseline.sha256; sha256sum -c _scratch/s15/reading_note_4l/note_4l_baseline.sha256
4d46c703ac54bc5d1f434c3ef2f822d40e70c9feffc9c8ae3d9d45ff82dc1ddd  research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md
8da34e5e672465319c399e1a016382801fbfacdaea2192b13db76e8e27b02acb  research/pde_ledger_v3/directives/S15_reading_decision_list.md
ccabfa672e91c771e9f6a07ecbf6649fd33e2139c519cc7b2fadc0412d0fd8aa  research/pde_ledger_v3/directives/S15_decision_list.md
305e21e0f544f7933b9651963693dd106f237e89d1e2ca66438c7776b3c64352  research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md: OK
research/pde_ledger_v3/directives/S15_reading_decision_list.md: OK
research/pde_ledger_v3/directives/S15_decision_list.md: OK
research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md: OK
```

```
$ sha256sum _scratch/s15/reading_note_4l/note_4l_review_prompt_p1.md; cat _scratch/s15/reading_note_4l/codex_review_p1.sha256 _scratch/s15/reading_note_4l/grok_review_p1.sha256; sha256sum -c _scratch/s15/reading_note_4l/codex_review_p1.sha256 _scratch/s15/reading_note_4l/grok_review_p1.sha256
19125607b44a8bb542539db63a187929b2f27417781a72c6b077e2ccb460514e  _scratch/s15/reading_note_4l/note_4l_review_prompt_p1.md
41b37635f2b4e92ee2dacfd2901a917b8df66b9d4c93458bcf95df412baf9601  /var/projects/_review_scratch/s15-4l-codex/REVIEW.md
030c6137c1285c80faaca76b44d3535c7da0cb2ab12c8f82ffe02b652b6c7249  /var/projects/_review_scratch/s15-4l-grok/REVIEW.md
/var/projects/_review_scratch/s15-4l-codex/REVIEW.md: OK
/var/projects/_review_scratch/s15-4l-grok/REVIEW.md: OK
```

```
$ git status --porcelain | /usr/bin/grep -v -F S15_reading_build_repair4l_note_lookups.md | diff - _scratch/s15/reading_note_4l/status_before_note_4l.txt && echo 'working tree status unchanged since before the legs launched, apart from this record'
working tree status unchanged since before the legs launched, apart from this record
```

```
$ head -14 /var/projects/_review_scratch/s15-4l-codex/codex_review.log | /usr/bin/grep -E '^(model|sandbox|reasoning effort|session id):'; cat /var/projects/_review_scratch/s15-4l-codex/DONE /var/projects/_review_scratch/s15-4l-grok/DONE
model: gpt-6.1-sol
sandbox: danger-full-access
reasoning effort: xhigh
session id: 01a1aa65-1023-7c60-bfbf-e5439d36aff9
exit 0
exit 0
```

Verdicts, as the last line of each report.

```
$ tail -1 /var/projects/_review_scratch/s15-4l-codex/REVIEW.md; tail -1 /var/projects/_review_scratch/s15-4l-grok/REVIEW.md
**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
**CLEAR**: nothing outstanding changes what is computed or what may be claimed.
```

Codex: evidence lines quoted in its report, found in its own stdout files.

```
$ /usr/bin/grep -c -F -- 'squared-distance residual: 0' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'recurrence residuals n=0..11: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'Gegenbauer(n,-s/2,mu) residuals n=0..12: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'supplied minus independent coefficients n=0..20: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'reindex m=n-k sign residual: 0' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'moving-centre squared-distance time-derivative residual: 0' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'FORM CONTROL: reversed distance cross term residual: 4*d*mu*r' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'COEFFICIENT CONTROL: wrong rising-factorial parameter at n=1 residual: 3*mu*s' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'maximum normalized coefficient error: 6.4655278e-82' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout:1
```

```
$ /usr/bin/grep -c -F -- 'PASS: 432 series values, including q=R/d=0.999 and mu=+/-1' /var/projects/_review_scratch/s15-4l-codex/derive_before_note.stdout /var/projects/_review_scratch/s15-4l-codex/check_note_equations.stdout /var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-codex/numeric_equations.stdout:1
```

Grok: evidence lines quoted in its report, found in its own output files.

```
$ /usr/bin/grep -c -F -- 'squared_residual_after_two_constraints 0' /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion.out /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out /var/projects/_review_scratch/s15-4l-grok/check_edge_convergence.out /var/projects/_review_scratch/s15-4l-grok/check_power_identity.out | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out:1
```

```
$ /usr/bin/grep -c -F -- 'flipped_mu_residual_after_two_constraints -4*r*(D1*e1 + D2*e2 + D3*e3)' /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion.out /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out /var/projects/_review_scratch/s15-4l-grok/check_edge_convergence.out /var/projects/_review_scratch/s15-4l-grok/check_power_identity.out | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out:1
```

```
$ /usr/bin/grep -c -F -- 'residual_note_term_minus_binomial_term 0' /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion.out /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out /var/projects/_review_scratch/s15-4l-grok/check_edge_convergence.out /var/projects/_review_scratch/s15-4l-grok/check_power_identity.out | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out:1
```

```
$ /usr/bin/grep -c -F -- 'log_residual 0' /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion.out /var/projects/_review_scratch/s15-4l-grok/check_supplied_expansion_2.out /var/projects/_review_scratch/s15-4l-grok/check_edge_convergence.out /var/projects/_review_scratch/s15-4l-grok/check_power_identity.out | /usr/bin/grep -v ':0$'
/var/projects/_review_scratch/s15-4l-grok/check_power_identity.out:1
```

The lines of the engines and the guard that the legs cite.

```
$ sed -n '309,310p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200; sed -n '1797,1798p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
distanceChart={Sqrt[1-mu^2],0,mu};
squaredDistance=trackedFullSimplify[(radial distanceChart-distance {0,0,1}).(radial distanceChart-distance {0,0,1}),Assumptions->Element[{radial,distance,mu},Reals]];
    # The squared distance and Jacobian come from the rotated chart itself.
    squared=S.trigsimp(S.expand((rr*rotated-dd*av).dot(rr*rotated-dd*av)))
```

### The review prompt, quoted in full

````markdown
# Independent physics review: S15 reading build, builder note 4l (a supplied expansion and a one-object trial)

## Artifact
`/var/projects/toy_physics/research/pde_ledger_v3/directives/resume_s15_reading_repair_4l.md`, uncommitted, sha256
`4d46c703ac54bc5d1f434c3ef2f822d40e70c9feffc9c8ae3d9d45ff82dc1ddd` (baseline in
`/var/projects/toy_physics/_scratch/s15/reading_note_4l/note_4l_baseline.sha256`). The orchestrator wrote it. Ignore
`.claude/skills/prep-for-compact/SKILL.md` and the other uncommitted engine files, which are not part of the artifact.

## What to check
The note goes to the Codex builder of the S15 reading build, which is in its fourth repair. Item 4 of the reading
decision list asks, on a power-law profile family with symbolic exponents, for the exact and complete conditions under
which each closed-form geometric departure object vanishes for every separation `D` on the domain. Exact constructions
of item 4 for the region objects that carry both drains have not returned in hours of runs. The user chose to supply an
expansion of the distance from drain 1 over the region, and to have the builder try a construction that uses it on one
object in each engine, before any further run. Both engines will read the supplied equations, so the note is
physics-bearing: an error in it can make both engines agree on the same wrong thing. Check:

1. **The supplied equations.** Are they true, on the region and the domain the S15 lists define, for every real
   exponent `s` (zero, negative, positive, integer and non-integer), with the convergence the note states? Is the
   stated relation to the Gegenbauer polynomials, and its stated range, correct? Is anything the supplied equations need,
   for a construction of item 4 on the named object, missing or misstated, for example about the region, the domain,
   the variables `z` and `μ`, or the moving centres?
2. **Object, not answer.** Does the note state or hint at what any condition contains, what any object equals, or what
   outcome is expected? Does it say how to derive anything beyond supplying the equations?
3. **The contract.** Does the note keep the object item 4 names, its exactness and its completeness, together with the
   decision list's sentence "No limit or expansion is taken in place of an exact object"? Could its wording let a
   finite number of series coefficients, or a truncated series, stand as the complete conditions? Does it conflict with
   the earlier notes it says stand (`resume_s15_reading_repair_4.md` and `4d` in
   `/var/projects/toy_physics/_scratch/s15/reading_build/`; `4f`–`4k` in `research/pde_ledger_v3/directives/`)?
4. **Operation.** Is the trial well defined: one object, both engines, the guard, the cap beside the two running SymPy
   attempts (the guard's admission rule is in `scripts/s11c_guarded_run.py:236–293`), the 1-hour trial stop, and a
   report that is enough to decide whether the construction is cheaper and complete?

## What you are handed
The working tree and history of `/var/projects/toy_physics`, read only. In particular:
- the note, the reading decision list `research/pde_ledger_v3/directives/S15_reading_decision_list.md`, and the S15
  decision list `research/pde_ledger_v3/directives/S15_decision_list.md`;
- the earlier repair-4 notes named in item 3 above;
- the read-only method assessment `research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md`
  (advisory and unreviewed; it is the source of the idea, not a premise);
- the engines `research/pde_ledger_v3/scripts/S15_separation_departure.py` and
  `research/pde_ledger_v3/mathematica/S15_separation_departure.wl`, and the accepted S15 build
  `research/pde_ledger_v3/scripts/S15_two_drain_region.py` and `research/pde_ledger_v3/mathematica/S15_two_drain_region.wl`.

## Required method
The artifact is a document. Read the sources first and form your own view of what they establish, then read the note.
Quote both sides for every finding.

A prose derivation is worth nothing: show the script and its literal stdout, or the claim is discarded. Derive the
supplied equations yourself before you rely on the note's statement of them, and test them by computation: symbolic
coefficient checks, and numeric checks at sample exponents, including zero, negative, positive even integer and
non-integer exponents, and at `z` close to `R/d`. Any physics claim you make must carry a script you wrote and its
literal output. Save every script and its literal stdout to named absolute paths under your scratch directory, and
report those paths. Do not write the construction the note asks the builder for.

## Physics filter
Report a finding only if it catches a way the physics could be wrong, or a way the note could make the builder compute
the wrong object or the record claim something it may not. Do not report style, and do not report "the note would be
wrong on a different input".

## Ablation sandbox
Work only under your scratch directory: `/var/projects/_review_scratch/s15-4l-codex/` for Codex,
`/var/projects/_review_scratch/s15-4l-grok/` for Grok. Change into it before you run any script, so that every file a
script writes lands there. Copy anything you want to alter there. Never modify the working tree of
`/var/projects/toy_physics`. Run no `datalad get`, no `git` command that writes, and nothing that changes the
repository. Do not stop, signal or inspect the memory of the running guarded units. If you run a CAS kernel, wrap it in
`timeout 600` and run one at a time; another review leg may be using a kernel on this machine.

## Bounds
Write your report to `REVIEW.md` in your scratch directory, and also give it as your final message. End it with one
of:
- **CLEAR**: nothing outstanding changes what is computed or what may be claimed;
- **FINDINGS**: a numbered list, each with severity, the quoted note text, the quoted source text, and your command and
  output.

Do not spawn agents or build `run_all`, `watch` or supervisor orchestration. Do not read the orchestrator's memory
(~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
````
