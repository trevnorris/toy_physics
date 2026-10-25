# Measurements: S15 comparator build brief, leak gate and contract check

Generator: `_scratch/s15/gen/s15_cmp_build_gate.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit is read by
the orchestrator.

## G0. The brief and the cleared list

```
$ sha256sum _scratch/s15/cmp_build/brief_s15_comparator_build.md; git show e305f6d4:research/pde_ledger_v3/directives/S15_comparator_decision_list.md | sha256sum
55ff33c50424aa7ebf1b5fdbb369a4ed43f626ce2b1a019b19816cec052c54ed  _scratch/s15/cmp_build/brief_s15_comparator_build.md
929fec5370a2a09a10732feb643a2463636b49e00b7b4b83158d87e09a006f21  -
```

```
$ git log -1 --format='%h %s'
e305f6d4 Clear the S15 comparator decision list (pass 1, both legs clear)
```

## G1. The contract is the cleared list's builder-side sections, byte for byte

```
$ diff <(sed -n '/^-\{100\}$/,/^-\{100\}$/p' _scratch/s15/cmp_build/brief_s15_comparator_build.md | sed '1d;$d' | sed '$d') <(git show e305f6d4:research/pde_ledger_v3/directives/S15_comparator_decision_list.md | sed -n '/^## Object/,$p') && echo 'contract identical to e305f6d4 sections Object .. Builder report'
contract identical to e305f6d4 sections Object .. Builder report
```

```
$ /usr/bin/grep -n -E '^## ' _scratch/s15/cmp_build/brief_s15_comparator_build.md
6:## Object
11:## Inputs (read-only)
25:## What must be true
250:## Exclusions
254:## Builder report (at most 40 lines)
```

```
$ /usr/bin/grep -c -E '^-{100}$' _scratch/s15/cmp_build/brief_s15_comparator_build.md
2
```

## G2. Outcome words and leg results (fixed strings)

```
$ /usr/bin/grep -n -i -F -- 'agree' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agree' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
52:2. **Join through an explicit binding table; a binding is not agreement.**
66:   - **A binding is never agreement.** Agreement comes only from the content compared beneath it.
75:     Totals are compared as well. Agreement of a total is never reported as agreement of its pieces, and agreement of
76:     a piece never as agreement of its groups. The engines' own partition residuals show completeness only: a term
173:   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
count: 5
```

```
$ /usr/bin/grep -n -i -F -- 'match' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'match' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
97:     - A structural mismatch is printed as a comparison that could not be formed. Examples: different component counts,
173:   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
count: 2
```

```
$ /usr/bin/grep -n -i -F -- 'identical' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'identical' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
64:       similarity of names, by position, or because two payloads are identical. Whether a binding is made never
137:     forms, applied identically to both operands. That covers each field of item 3, both centres and both radii, with
143:     - In every form, no test function is constant in any of its arguments, and none vanishes identically on a surface
count: 3
```

```
$ /usr/bin/grep -n -i -F -- 'equal' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'equal' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
40:     - relations and equations (`Eq(…)` against `==`, a conjunction of inequalities against a chained inequality);
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'zero' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'zero' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
70:       including a group an engine prints as zero;
113:   - A mapped name is a substitution, never by itself a zero. A repoint control shows that pointing a mapped name at a
174:   and nonzero only on operational failure (a missing or ungrammatical input). Interpreting the output belongs to the
count: 3
```

```
$ /usr/bin/grep -n -i -F -- 'e-1' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'e-1' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'e-0' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'e-0' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'CLEAR' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'CLEAR' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
1:Status: the S15 comparator decision list has passed its review (commit e305f6d4), so you are cleared to PROCEED.
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'sound' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'sound' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'close' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'close' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'expected' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'expected' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
173:   expected to be. It emits no `PASS`/`FAIL`/`AGREE`/`MATCH`/`VERDICT`/`STATUS` token. It exits 0 whatever it finds,
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'should' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'should' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '309' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '309' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '315' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '315' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '318' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '318' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '4.1' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- '4.1' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'derivation' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'derivation' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'H1' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'H1' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'reviewed' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'reviewed' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

## G3. The builder stays in its lane (no review, agent or downstream instruction)

```
$ /usr/bin/grep -n -i -F -- 'review' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'review' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
1:Status: the S15 comparator decision list has passed its review (commit e305f6d4), so you are cleared to PROCEED.
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'agent' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
3:You are the author of the S15 cross-engine comparator, working in /var/projects/toy_physics. Read AGENTS.md first. Your contract is the text between the two rules below
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Claude' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
278:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' _scratch/s15/cmp_build/brief_s15_comparator_build.md | cut -c1-170; echo "count: $(/usr/bin/grep -c -i -F -- 'record' _scratch/s15/cmp_build/brief_s15_comparator_build.md)"
31:   - The SymPy stream repeats some tags as records (L2, L3): `dimension_argument/name`, `/operand_a`, `/operand_b` and
32:     `/residuals` (208 records), and `partition/unlabelled_coefficient` (27 lines). Each record, and each repeated
61:     - Each binding joins one end in each stream. An end is a tag, a record of item 1, or a component of a payload
79:     - Each tag occurs once per stream, apart from progress markers and item 1's records; any other repeat is
81:     - Every tag and record is bound, or is listed as unjoined with its engine and the reason. Unjoined tags may be
175:   S15 record.
206:   - in a repeated record of item 1, one record's payload changed, and one record removed;
216:   - split-against-grouped packaging: one side prints under several tags or records what the other prints as
235:   - item 1's repeated records printed in another order.
276:When the comparator and its tests are written, run the tests, then run the full comparison on the two accepted transcripts with the final code, each through the guard
count: 10
```

## G4. The guard path resolves under _scratch/s11c/

```
$ /usr/bin/grep -n -o -E -- '--log-directory [^ ]+' _scratch/s15/cmp_build/brief_s15_comparator_build.md; ls -ld _scratch/s11c/s15_cmp_guard; readlink -f _scratch/s11c/s15_cmp_guard
271:--log-directory /var/projects/toy_physics/_scratch/s11c/s15_cmp_guard/<a
drwxrwxr-x 2 trevnorris trevnorris 4096 Oct 24 20:39 _scratch/s11c/s15_cmp_guard
/var/projects/toy_physics/_scratch/s11c/s15_cmp_guard
```
