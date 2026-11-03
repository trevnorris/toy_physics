# Measurements: S15 reading build repair-4h note (user: try a faster method first), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes. Every hit
is read by the orchestrator. The run and work directories named in these commands are local scratch and disposable;
every output this note relies on is quoted here in full. The note reports that the orchestrator stopped the guarded run
`repair3_90441_wolfram_final_20` on 2026-11-03 at the user's decision, and carries that decision: before any final run,
look for an exact construction of item 4 that costs much less, try it on the six region objects that carry both drains, in
each engine, and report. It points at an advisory assessment committed beside it
(`S15_reading_item4_method_assessment.md`), and back at `resume_s15_reading_repair_4.md` (gated at `5878692b`), `_4d`
(`6f46583b`), `_4f` (`74fb5700`) and `_4g` (`1a703b67`). The note names no value, exponent or condition. The assessment
describes methods (it is about how item 4 is computed) and names no condition; R1/R2 below probe both files for values.

## The user's words and decisions (2026-11-03)

1. AskUserQuestion, ~11:05Z. Answer chosen: **"Keep going (Recommended)"**: "No change: stop only at 6 h or out of memory,
   one engine at a time. Fewest gaps. About 2–4 days before the review round. I keep counting and ask again once SymPy's pace
   is known."
2. The user, ~17:10Z, after asking for the average time per object: "If we hit another 6 hour stop, feels like that'll be
   50% failure and lost 18 hours. Are we sure we're computing these correctly? Maybe there's a faster way to do this. Have
   we asked that question?"
3. The orchestrator then obtained the read-only assessment. AskUserQuestion, ~17:55Z. Answer chosen: **"Try faster method
   first (Recommended)"**: "Stop the Wolfram run now (object 6, ~1 h in). The builder tries a faster exact method on the 6 hard
   objects and fixes the 2 problems, then reports times to you before any full rerun. If it works, both engines rerun with it,
   then review. The current run's output would likely be redone anyway because of the 2 problems."

## R0. The note and the assessment

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
924ae2216861b6a2d44fa3fee92ec5a64a495ad358e86faea42149bedc800813  research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
305e21e0f544f7933b9651963693dd106f237e89d1e2ca66438c7776b3c64352  research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
 4380 research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
17950 research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
22330 total
```

```
$ sha256sum /var/projects/_review_scratch/s15-item4-method-consult/prompt.md /var/projects/_review_scratch/s15-item4-method-consult/ANSWER.md; cat /var/projects/_review_scratch/s15-item4-method-consult/DONE; /usr/bin/grep -c -F 'session id: 01a1a280-56e3-7da0-8d0b-f9924ffe67de' /var/projects/_review_scratch/s15-item4-method-consult/codex.log; /usr/bin/grep -n -F -e 'model: gpt-6.1-sol' -e 'sandbox: read-only' -e 'reasoning effort: xhigh' /var/projects/_review_scratch/s15-item4-method-consult/codex.log | head -3
69f980d6ce9ae9086e311241cc9f7760a8c915071af84705a5fc2fb6434cb69b  /var/projects/_review_scratch/s15-item4-method-consult/prompt.md
e65338ec2be035204ff862d88100616fd4c6e41f0e0ec980bd3cb2f96843949b  /var/projects/_review_scratch/s15-item4-method-consult/ANSWER.md
exit 0
1
4:model: gpt-6.1-sol
7:sandbox: read-only
8:reasoning effort: xhigh
```

```
$ diff <(sed -n '17,$p' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md) /var/projects/_review_scratch/s15-item4-method-consult/ANSWER.md
0a1,4
> codex
> The new constructions **can return exponent conditions**, so they have moved beyond the incapable attempts described in Y2. They do not yet establish delivery of the complete item-4 classification. The dominant concern is excessive symbolic representation and parameter branching; exactness itself does not require this approach.
> 
> I could not save `ANSWER.md` because the filesystem sandbox rejected the write. The assessment follows.
7c11
< - **DRAFT**: `sympy_six_hour_candidate.py` (the builder's SymPy working draft)
---
> - **DRAFT**: [sympy_six_hour_candidate.py](/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/round90441/sympy_six_hour_candidate.py)
129c133
< The attempted write of the answer file was rejected by the read-only sandbox and disabled approval settings; the file was not created.
---
> The attempted write to `/var/projects/_review_scratch/s15-item4-method-consult/ANSWER.md` was rejected by the read-only sandbox and disabled approval settings; the file was not created.
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note and in the assessment

### research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_20` through its unit (`systemctl --user stop s11c-guard-1af66ae632bc`), with its launcher alive; the journal records the stop at 17:57:00Z. The r
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
9:None of the six region objects that carry both drains has returned its conditions in any run of this repair: two ran out of memory within minutes, three were stopped at 6 hours or more, and I stoppe
count: 1
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

### research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
1:# S15 reading item 4: an outside assessment of how it is computed (2026-11-03)
3:A read-only assessment the orchestrator obtained on 2026-11-03, after the user asked whether item 4 is computed correctly
7:- Source: a fresh `codex exec` session (model gpt-6.1-sol, reasoning xhigh, sandbox read-only), session id `01a1a280-56e3-7da0-8d0b-f9924ffe67de`, 2026-11-03 17:40–17:52Z.
count: 3
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
91:   G(z,\mu)=(1-2\mu z+z^2)^\gamma
95:   (1-2\mu z+z^2)G'=\gamma(-2\mu+2z)G.
count: 2
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
95:   (1-2\mu z+z^2)G'=\gamma(-2\mu+2z)G.
count: 1
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
24:- **DL**: [S15_reading_decision_list.md](/var/projects/toy_physics/research/pde_ledger_v3/directives/S15_reading_decision_list.md)
25:- **R2/R3**: the corresponding build-review disposition files in `research/pde_ledger_v3/directives/_measurements/`.
133:An all-orders convergent representation with a proved tail identity also does not impose a small-radius regime: throughout the domain, \(r/d<1\). If the intended rule instead forbids any internal 
count: 3
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
63:- **Increasing series order and matrix size:** Wolfram generates rows until there are as many pivots as columns (`WL:880–891`). SymPy repeatedly differentiates and computes rank toward the same f
count: 1
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
7:- Source: a fresh `codex exec` session (model gpt-6.1-sol, reasoning xhigh, sandbox read-only), session id `01a1a280-56e3-7da0-8d0b-f9924ffe67de`, 2026-11-03 17:40–17:52Z.
count: 1
```

```
$ /usr/bin/grep -n -F -- 'b1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b1' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '2a' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '2a' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q=' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q=' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'q =' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'q =' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'w=' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'w=' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
61:- **Larger parameter problems than the physical problem:** Wolfram introduces coordinates for parameter powers, logarithms, amplitude monomials and then products of these coordinates (`WL:1081–11
count: 1
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
137:**5. Item 4 is meaningful; its wording does not force the current cost.**
count: 1
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
27:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'record' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_20` through its unit (`systemctl --user stop s11c-guard-1af66ae632bc`), with its launcher alive; the journal records the stop at 17:57:00Z. The r
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
9:None of the six region objects that carry both drains has returned its conditions in any run of this repair: two ran out of memory within minutes, three were stopped at 6 hours or more, and I stoppe
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
9:None of the six region objects that carry both drains has returned its conditions in any run of this repair: two ran out of memory within minutes, three were stopped at 6 hours or more, and I stoppe
count: 1
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'numeric' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'RAM' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md)"
1:Note from the orchestrator. I stopped `repair3_90441_wolfram_final_20`, and I have read your report of that stop (`repair4_external_stop_final20.md`). The user has made a further decision about item
5:On 2026-11-03 I stopped `repair3_90441_wolfram_final_20` through its unit (`systemctl --user stop s11c-guard-1af66ae632bc`), with its launcher alive; the journal records the stop at 17:57:00Z. The r
24:- Try the construction on the six objects `region/component_1/momentum_labels_{drain1, drain2}`, `region/component_1/rate_labels_{drain1, drain2}`, `region/component_2/momentum_labels_{drain1, drai
count: 3
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4d.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4f.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4g.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'S15_reading_item4_method_assessment.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'repair4_external_stop_final20.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
2
```

```
$ /usr/bin/grep -c -F -- 'memory-gib 55' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'one job at a time' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- '6 hours' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
3
```

```
$ /usr/bin/grep -c -F -- 'Wolfram first' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not start the final runs' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'poll each job until it ends' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
1
```

```
$ ls _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md
research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
```

```
$ /usr/bin/grep -n -F 'What must be true' _scratch/s15/reading_build/resume_s15_reading_repair_4.md _scratch/s15/reading_build/resume_s15_reading_repair_4d.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md
_scratch/s15/reading_build/resume_s15_reading_repair_4.md:18:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:7:## What must be true
_scratch/s15/reading_build/resume_s15_reading_repair_4d.md:9:Everything under "What must be true" in `resume_s15_reading_repair_4.md` stands. For an object whose attempt ran out of memory, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4f.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md` and in `resume_s15_reading_repair_4d.md` stands. For an object whose attempt was stopped at 6 hours, these are also true:
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:13:## What must be true
research/pde_ledger_v3/directives/resume_s15_reading_repair_4g.md:15:Everything under "What must be true" in `resume_s15_reading_repair_4.md`, `resume_s15_reading_repair_4d.md` and `resume_s15_reading_repair_4f.md` stands. This decision replaces the bullet in `resume_s15_reading_repair_4f.md` about an attempt that goes silent. Also:
```

The source lines the note cites, retrieved verbatim:

```
$ sed -n '1390,1394p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200
 solo[prefix<>"/unresolved_parameter_branches",nativeAnswer[[2]]/.inverse];
 If[nativeAnswer[[2]]=!={},
  emit[prefix<>"/status","NOT OBTAINED EXACTLY on the displayed branches after native Euler radial lifting, primitive-support partition, differential elimination, series row reduction, Euler-spectrum 
  AppendTo[exactGaps,{name,"item4",nativeAnswer[[2]]/.inverse}],
  emit[prefix<>"/status","The exponent-sensitive traversal returned the printed conditions, restricted Q objects, and unresolved-branch list."]];
```

```
$ sed -n '2143,2148p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
            solo(tag+'/unresolved_parameter_branches',unresolved)
            solo(tag+'/condition_count',len(answers))
            if unresolved:
                reason='NOT OBTAINED EXACTLY on the displayed unresolved parameter branches: Euler radial lifting of the printed closed Q, primitive-support partition, exact differential elimination, 
                emit(tag+'/status',reason)
                gaps.append((name,'item4',reason))
```

```
$ sed -n '444,462p;831,847p;1094,1114p;1235,1240p;1277,1287p;1307,1309p' research/pde_ledger_v3/scripts/S15_separation_departure.py | /usr/bin/grep -n -F -e 'self.hi' -e 'self.lo' -e 'subs(rules' -e 'radial_moment' | cut -c1-200
30:                result+=self.call(scope+'/integrate_radial_remainder_'+str(index),S.integrate,term,(self.r,self.lo,self.hi),conds='piecewise')
32:            if self.radial_moment is None:
33:                self.radial_moment=self.call(scope+'/integrate_symbolic_radial_moment',S.integrate,self.r**self.moment_exponent,(self.r,self.lo,self.hi),conds='piecewise')
34:            moment=self.radial_moment.subs(self.moment_exponent,exponent)
35:            self.show(scope+'/radial_moment_substitution',(self.radial_moment,self.moment_exponent,exponent,moment))
37:                if set(variables) <= {self.lo,self.hi}:
39:                        radial_polynomial=S.Poly(base,self.lo,self.hi,domain='QQ')
45:                        rescaled=S.cancel(base.subs(self.hi,self.lo*ratio)/self.lo**degree)
47:                        self.show(scope+'/radial_ratio_operands',(base,{self.hi:self.lo*ratio},degree,ratio_polynomial.as_expr()))
48:                        self.show(scope+'/radial_ratio_reconstruction_residual',S.cancel(base-self.lo**degree*ratio_polynomial.as_expr().subs(ratio,self.hi/self.lo)))
56:                                if admitted is S.true:candidates.append({self.hi:self.lo*root})
64:        support=S.simplify(support.subs(rules,simultaneous=True))
67:        q=q.subs(rules,simultaneous=True)
73:            h,b=(self.normalize(value.subs(rules,simultaneous=True),scope+'/radial_restriction') for value in radial)
77:            evaluated=S.expand_log(evaluated,force=True).subs(rules,simultaneous=True)
```

```
$ sha256sum research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1  research/pde_ledger_v3/mathematica/S15_separation_departure.wl
48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f  research/pde_ledger_v3/scripts/S15_separation_departure.py
b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

## R5. Facts the note states

S0, captured by the orchestrator at 2026-11-03T17:56Z, before the stop (commands and literal output):

````
2026-11-03T17:56:54.544696521Z
● s11c-guard-1af66ae632bc.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/toy_physics/_scratch/s11c/s15/repair3_90441_wolfram_final_20/invocation.json
     Loaded: loaded (/run/user/1000/systemd/transient/s11c-guard-1af66ae632bc.service; transient)
  Transient: yes
     Active: active (running) since Mon 2026-11-02 21:43:27 MST; 13h ago
   Main PID: 3086480 (python3)
      Tasks: 20 (limit: 32)
     Memory: 8.7G (max: 55.0G swap max: 0B available: 31.4G peak: 35.7G)
        CPU: 13h 13min 6.114s
     CGroup: /user.slice/user-1000.slice/user@1000.service/app.slice/s11c-guard-1af66ae632bc.service
             ├─3086480 /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/toy_physics/_scratch/s11c/s15/repair3_90441_wolfram_final_20/invocation.json
             ├─3086482 python3 /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure_run.py --source /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl --transcript /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out --work-directory /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/round90441/wolfram_final20_execution
             └─3490680 /var/projects/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel -script /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departure.wl
anon 6082359296
inactive_file 3254177792
38382227456
max 0
oom 0
oom_kill 0
3085929  0.1 237880   13:13:48 codex
3086407  0.0 33720    13:13:38 codex-code-mode
3490680 99.8 6017688  01:05:46 WolframKernel
17669303 2026-11-03 10:53:21.061172911 -0700 _scratch/s11c/s15/repair3_90441_wolfram_final_20/stdout
24021 2026-11-03 09:51:08.541188332 -0700 _scratch/s15/reading_build_work/repair3/round90441/wolfram_final20_execution/execution.json
4170217 2026-11-03 10:56:12.866361017 -0700 _scratch/s15/reading_build/repair4g.log
WORKING_ON: reading/region/component_3/rate_labels_{drain1, drain2} / item4 exponent-sensitive / PolynomialReduce normalization direction monomial
b7f32997595560494d2bf9951aa145f65d3f973408459cb438cf985d9e612302  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
4f26580821bb7333bca27718c29f8cff4d32056880359a04e9eb9db2754e51e1  research/pde_ledger_v3/mathematica/S15_separation_departure.wl
48768927d30b0e1458cf7c9758ff9bbf8378a6f1473816e4a72a29d9602a142f  research/pde_ledger_v3/scripts/S15_separation_departure.py
````

```
$ journalctl --user -u s11c-guard-1af66ae632bc.service --no-pager -o short-iso-precise --since '2026-11-03 17:56:30 UTC' --until '2026-11-03 17:58:00 UTC' --utc | cut -c1-200
2026-11-03T17:57:00.061900+00:00 pop-os systemd[1646]: Stopping s11c-guard-1af66ae632bc.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-03T17:57:05.308061+00:00 pop-os systemd[1646]: s11c-guard-1af66ae632bc.service: State 'final-sigterm' timed out. Killing.
2026-11-03T17:57:05.308148+00:00 pop-os systemd[1646]: s11c-guard-1af66ae632bc.service: Killing process 3490680 (WolframKernel) with signal SIGKILL.
2026-11-03T17:57:05.308247+00:00 pop-os systemd[1646]: s11c-guard-1af66ae632bc.service: Failed with result 'timeout'.
2026-11-03T17:57:05.308275+00:00 pop-os systemd[1646]: s11c-guard-1af66ae632bc.service: Unit process 3490680 (WolframKernel) remains running after unit stopped.
2026-11-03T17:57:05.308519+00:00 pop-os systemd[1646]: Stopped s11c-guard-1af66ae632bc.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-03T17:57:05.308671+00:00 pop-os systemd[1646]: s11c-guard-1af66ae632bc.service: Consumed 13h 13min 16.875s CPU time, 35.7G memory peak, 0B memory swap peak.
```

```
$ cat _scratch/s11c/s15/repair3_90441_wolfram_final_20/outcome.json 2>&1; echo; /usr/bin/grep -c -F 'No more memory available' _scratch/s11c/s15/repair3_90441_wolfram_final_20/stdout
{
  "exitCode": 1,
  "wallSeconds": 47618.05529435002,
  "unit": "s11c-guard-1af66ae632bc",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}

0
```

```
$ python3 -c "import json;s=json.load(open('/tmp/claude-1000/-var-projects-toy-physics/a418167f-6971-4b74-ba3d-d85e8401ca3f/scratchpad/s15_4h/execution_final20_at_stop.json'));k=[x for x in s['kernels'] if isinstance(x,dict) and x.get('mode')=='item4'];print([(x['object'],x['started'],x.get('ended'),x.get('author_stop_kind')) for x in k]);print({n:v['kind'] for n,v in s['failed_objects'].items()})"
[('region/component_1/momentum_labels_{drain1, drain2}', '2026-11-03T04:46:40.529865+00:00', '2026-11-03T04:46:55.534773+00:00', None), ('region/component_1/rate_labels_{drain1, drain2}', '2026-11-03T04:46:55.539523+00:00', '2026-11-03T04:48:04.550821+00:00', None), ('region/component_2/momentum_labels_{drain1, drain2}', '2026-11-03T04:48:04.555622+00:00', '2026-11-03T04:48:19.560430+00:00', None), ('region/component_2/rate_labels_{drain1, drain2}', '2026-11-03T04:48:19.564251+00:00', '2026-11-03T10:49:27.850703+00:00', 'six_hour_stop'), ('region/component_2/rate_labels_{drain1, drain2}', '2026-11-03T10:49:28.606502+00:00', '2026-11-03T10:50:35.618996+00:00', None), ('region/component_3/momentum_labels_{drain1, drain2}', '2026-11-03T10:50:35.623158+00:00', '2026-11-03T16:50:52.527968+00:00', 'six_hour_stop'), ('region/component_3/momentum_labels_{drain1, drain2}', '2026-11-03T16:50:52.529322+00:00', '2026-11-03T16:51:08.533484+00:00', None), ('region/component_3/rate_labels_{drain1, drain2}', '2026-11-03T16:51:08.538656+00:00', None, None)]
{'region/component_1/momentum_labels_{drain1, drain2}': 'memory_exhausted', 'region/component_1/rate_labels_{drain1, drain2}': 'memory_exhausted', 'region/component_2/momentum_labels_{drain1, drain2}': 'six_hour_stop', 'region/component_2/rate_labels_{drain1, drain2}': 'six_hour_stop', 'region/component_3/momentum_labels_{drain1, drain2}': 'six_hour_stop'}
```

```
$ ls -l _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final20.md | awk '{print $5, $9}'; sha256sum _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final20.md; /usr/bin/grep -n -F -e 'recorded peak' -e 'external service stop' -e 'Unfinished object' -e 'WORKING_ON' -e 'do not attempt' _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final20.md | cut -c1-220
2590 _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final20.md
605ef512c59ea54ca1b6e0b518a124378310103314f415652a1e9ba69bc58920  _scratch/s15/reading_build_work/repair3/round90441/repair4_external_stop_final20.md
9:Guard exit **1**, child outcome `null`, wall time **47618.05529435002 s**, recorded peak **38382227456 bytes**. The last sample reports `max 0`, `oom 0`, `oom_kill 0`, `oom_group_kill 0`. The journal records an externa
11:Unfinished object: `region/component_3/rate_labels_{drain1, drain2}`, kernel started at `2026-11-03T16:51:08.538656Z`. It had not reached six hours and had not run out of memory. Last reported step:
14:WORKING_ON: reading/region/component_3/rate_labels_{drain1, drain2} / item4 exponent-sensitive / PolynomialReduce normalization direction monomial
24:Preserve these two failures on continuation; do not attempt those objects again. The earlier two memory labels and final18 six-hour label were retained.
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json; echo; pgrep -c WolframKernel; pgrep -c codex; ls _scratch/s15/reading_build/repair4g.DONE 2>&1
{}

0
2
ls: cannot access '_scratch/s15/reading_build/repair4g.DONE': No such file or directory
```

```
$ git diff --quiet 90441b89 -- research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md; echo "run record diff vs 90441b89 exit $?"; git status --short -- research/pde_ledger_v3
run record diff vs 90441b89 exit 0
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4h_note_lookups.md
?? research/pde_ledger_v3/directives/_measurements/S15_reading_item4_method_assessment.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4h.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```
