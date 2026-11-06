# Measurements: S15 reading build repair-4o note (user: SymPy attempts no item-4 object), leak gate and fact check

Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts, hashes, process and
unit state. Every hit is read by the orchestrator. The run and work directories named in these commands are local scratch
and disposable; every output this note relies on is quoted here in full. The note carries a user decision that changes
what SymPy computes for item 4 and what the S15 run record may claim (a fourth label kind), so it is physics-bearing and
gets two review legs before commit.

## The user's decision (2026-11-06)

The user's words, verbatim: "Feels like we need to stop using sympy for this work. It was great for symbolic work, but
without being compiled it's just too slow"

AskUserQuestion: "What should happen to SymPy in the S15 build that's running now?" The user selected **"Drop it from item 4
(Recommended)"**, whose text read: "Stop SymPy's final run. SymPy still prints the other items (1–3, 5, 6), so those keep
their two-engine check. Each item-4 object gets a 'not attempted in SymPy, user decision' label, and Wolfram's run
continues untouched. Item 4 stands on Wolfram alone until a replacement engine checks it. I'd write a short note for the
builder and put it through two reviewers before using it."

AskUserQuestion: "What should replace SymPy as the second engine for heavy algebra?" The user selected **"Wait for the
Singular test (Recommended)"**. (Not in the note: it changes nothing the builder computes.)

## The stops, as run (literal output of the two commands, copied from the session at the time)

```
$ cd /var/projects/toy_physics && echo "unit check:"; systemctl --user show -p ExecStart s11c-guard-06d785790b54.service | /usr/bin/grep -a -o 'repair4_decision_final_[a-z_0-9]*' | head -1; TZ=UTC date '+stop requested %Y-%m-%dT%H:%M:%SZ'; systemctl --user stop s11c-guard-06d785790b54.service; TZ=UTC date '+stop returned %Y-%m-%dT%H:%M:%SZ'; for i in $(seq 1 30); do [ -e _scratch/s11c/s15/repair4_decision_final_sympy_02/outcome.json ] && break; sleep 1; done; cat _scratch/s11c/s15/repair4_decision_final_sympy_02/outcome.json; echo; pgrep -fa run_sympy02.py | /usr/bin/grep -a -v pgrep | cut -c1-80
unit check:
repair4_decision_final_sympy_02
stop requested 2026-11-06T17:42:28Z
stop returned 2026-11-06T17:42:28Z
{
  "exitCode": 1,
  "wallSeconds": 7309.415392439114,
  "unit": "s11c-guard-06d785790b54",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}

2038153 python3 _scratch/s15/reading_build_work/repair3/final_decision/run_sympy

$ cd /var/projects/toy_physics && TZ=UTC date '+kill %Y-%m-%dT%H:%M:%SZ'; kill -TERM 2035710; for i in $(seq 1 15); do pgrep -f 'codex exec resume' >/dev/null || break; sleep 1; done; pgrep -fa 'codex exec resume' | /usr/bin/grep -a -v pgrep | cut -c1-60 || true; echo "--- after:"; ps -eo pid,ppid,etime,args | /usr/bin/grep -a -E 'run_parallel|run_sympy02|wait_event|s11c_guarded_run|wolfram_01_runner|WolframKernel' | /usr/bin/grep -a -v grep | cut -c1-150; systemctl --user list-units --no-legend 's11c*' | awk '{print $1,$3,$4}'; cat ~/.local/state/s11c-guard/resource-pool/reservations.json | python3 -c 'import json,sys; d=json.load(sys.stdin); print(list(d))'
kill 2026-11-06T17:42:36Z
--- after:
2037429    1646    02:02:51 python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 8 --log-directory /var/projects/toy_
2037432 2037429    02:02:51 systemd-run --user --quiet --wait --pipe --collect --unit=s11c-guard-b4d073f2bd93 --property=MemoryMax=8589934592 --proper
2037433    1646    02:02:51 /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/toy_physics/_scratch
2037435 2037433    02:02:51 python3 /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01
2043520 2037435    01:55:03 /var/projects/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel -script /var/projects/toy_physic
s11c-guard-b4d073f2bd93.service active running
['s11c-guard-b4d073f2bd93.service']
(PID 2035710 was the `codex exec resume … 01a18cb6…` process, as listed by `ps` at 17:3xZ in the session: "2035710 2035709 2035710 2035709 01:59:12 codex exec resume -m gpt-6-astra …".)
```

## R0. The note

```
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md; wc -c research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1196390b199b67509b3cd65096a5b86d35a15f27a1d967c70a38071ea46d7a22  research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
2652 research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
```

## R1. Values and exponents (fixed strings, from the repair-4 gate), in the note

```
$ /usr/bin/grep -n -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'alpha' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'α' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'λ' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
5:On 2026-11-06, on the user's decision:
count: 1
```

```
$ /usr/bin/grep -n -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '−3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '-3' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '+2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '= 0' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '1/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd²' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd^2' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd**' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'd³' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- '/d' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'inverse' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'square' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b_1' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'b₁' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'aa[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'bb[1]' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'a+b' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Pizzetti' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'Lagally' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'pull' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'mass' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'force' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'gravity' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'attract' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -F -- 'series' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

## R2. Symbols in proximity: every line naming an exponent or amplitude symbol

```
$ /usr/bin/grep -n -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-220; echo "count: $(/usr/bin/grep -c -E '(^|[^A-Za-z])(a_k|b_k|a_[12]|b_[12]|A_[12]|B_[12]|aa\[|bb\[|AA\[|BB\[)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

## R3. The builder stays in its lane (the note)

```
$ /usr/bin/grep -n -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Claude' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
31:Do not read the orchestrator's memory (~/.claude/projects/) or any other checkout of this repository (the /var/projects/toy_physics-* directories).
count: 1
```

```
$ /usr/bin/grep -n -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Grok' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Codex' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'agent' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'spawn' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'comparator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'review' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'REVIEW' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_review_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- '_scratch' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'disposition' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'lookups' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'leg' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Singular' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Singular' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'Maple' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'Maple' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

```
$ /usr/bin/grep -n -i -F -- 'engine' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md | cut -c1-200; echo "count: $(/usr/bin/grep -c -i -F -- 'engine' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md)"
count: 0
```

## R4. Pointers and operational sentences

```
$ /usr/bin/grep -c -F -- 'resume_s15_reading_repair_4n.md' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
4
```

```
$ /usr/bin/grep -c -F -- 'What must be true' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
2
```

```
$ /usr/bin/grep -c -F -- 'fourth kind of label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'the third kind of' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'No Wolfram object carries it' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'not attempted in SymPy' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'items 1, 2, 3, 5 and 6' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'as history' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'Neither gives any object a label' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not stop, signal or change it' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'near its own expected peak' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'make every attempt again' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'Stop there' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -c -F -- 'Do not read the orchestrator' research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
1
```

```
$ /usr/bin/grep -n -F -e '## What must be true' -e 'third kind of label' -e 'An item-4 attempt in a final run is labelled only' -e 'make every attempt again that the rule above requires' -e 'the updated run record and a short final report' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md | cut -c1-200
18:## What must be true
21:- In each engine's final run, each of the six objects is printed with a not-obtained label that says it was set by this decision and that the object was not attempted in that run. This is a third k
23:- An item-4 attempt in a final run is labelled only when it is still running 6 hours after it started, or when it runs out of memory. It runs out of memory when it ends for lack of memory, whether 
29:- Then make the final guarded run of each engine from the final code. Run the two at the same time, each with a `--memory-gib` near its own expected peak, the two together within what the guard adm
30:- Then finish as the brief says: the updated run record and a short final report. Stop there; the orchestrator will reply.
```

```
$ git log -1 --format='%h %s' -- research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
266f2acf Accept the S15 reading build repair-4n note after two CLEAR legs (user: stop the SymPy attempts, label the six, finish)
```

```
$ /usr/bin/grep -n -F -e '1. `Q` itself' -e '4. For each `Q` whose item-3 object' -e '5. For each drain' -e '6. Dimensional homogeneity' research/pde_ledger_v3/directives/S15_reading_decision_list.md | cut -c1-120
72:1. `Q` itself.
78:4. For each `Q` whose item-3 object is in closed form:
83:5. For each drain `k`, the part of `Q` that depends on that drain's velocity: `Q` minus `Q` at `Ẋ_k = 0`.
84:6. Dimensional homogeneity of items 2 and 5, each change `F′` carrying the dimension of `F`.
```

## R5. Facts the note states

The stopped SymPy run: its guard files and outcome, and the unit's journal around the stop.

```
$ ls _scratch/s11c/s15/repair4_decision_final_sympy_02; cat _scratch/s11c/s15/repair4_decision_final_sympy_02/outcome.json
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
outcome.json
resource-samples.jsonl
stderr
stdout
{
  "exitCode": 1,
  "wallSeconds": 7309.415392439114,
  "unit": "s11c-guard-06d785790b54",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}
```

```
$ python3 -c "import json;i=json.load(open('_scratch/s11c/s15/repair4_decision_final_sympy_02/invocation.json'));print(i.get('unit'), i.get('command'))" | cut -c1-400
s11c-guard-06d785790b54 ['env', 'S15_PY_TRANSCRIPT=/var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_separation_departure.out', 'S15_PY_WORK_DIRECTORY=/var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_sympy_02_work', 'S15_PY_RUN_NAME=repair4_decision_final_sympy_02', 'python3', '-u', '/var/projects/toy_physics/_scratch/s15/readin
```

```
$ journalctl --user -u s11c-guard-06d785790b54.service --since '2026-11-06 17:42:00 UTC' --until '2026-11-06 17:43:00 UTC' -o short-iso-precise --no-pager | cut -c1-200
2026-11-06T17:42:28.019213+00:00 pop-os systemd[1646]: Stopping s11c-guard-06d785790b54.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/
2026-11-06T17:42:28.041332+00:00 pop-os systemd[1646]: Stopped s11c-guard-06d785790b54.service - /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/t
2026-11-06T17:42:28.041589+00:00 pop-os systemd[1646]: s11c-guard-06d785790b54.service: Consumed 2h 1min 48.041s CPU time, 607.0M memory peak, 0B memory swap peak.
```

```
$ ls _scratch/s11c/s15/repair4_decision_final_sympy_01; python3 -c "import json;o=json.load(open('_scratch/s11c/s15/repair4_decision_final_sympy_01/outcome.json'));print({k:o.get(k) for k in ('exitCode','wallSeconds')})"
child-outcome.json
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
outcome.json
resource-samples.jsonl
stderr
stdout
{'exitCode': 1, 'wallSeconds': 0.5348064810968935}
```

The session end: its marker, and what is still running.

```
$ cat _scratch/s15/reading_build/repair4n.DONE
exit 143
```

```
$ for s in run_parallel.py run_sympy02.py wait_event.py 'codex exec resume'; do printf '%s: %s\n' "$s" "$(pgrep -c -f "$s")"; done
run_parallel.py: 0
run_sympy02.py: 0
wait_event.py: 0
codex exec resume: 0
```

```
$ ps -eo pid,ppid,lstart,args | /usr/bin/grep -a -E 's11c_guarded_run|wolfram_01_runner|WolframKernel' | /usr/bin/grep -a -v grep | cut -c1-200
2037429    1646 Fri Nov  6 15:39:59 2026 python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15 --memory-gib 8 --log-directory /var/projects/toy_physics/_scratch/s15/guard/repair4_de
2037432 2037429 Fri Nov  6 15:39:59 2026 systemd-run --user --quiet --wait --pipe --collect --unit=s11c-guard-b4d073f2bd93 --property=MemoryMax=8589934592 --property=MemorySwapMax=0 --property=TasksMa
2037433    1646 Fri Nov  6 15:39:59 2026 /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --child-manifest /var/projects/toy_physics/_scratch/s11c/s15/repair4_decision_final_wolf
2037435 2037433 Fri Nov  6 15:39:59 2026 python3 /var/projects/toy_physics/_scratch/s15/reading_build_work/repair3/final_decision/repair4_decision_final_wolfram_01_runner.py --source /var/projects/toy
2043520 2037435 Fri Nov  6 15:47:47 2026 /var/projects/Wolfram/Wolfram/15.0/SystemFiles/Kernel/Binaries/Linux-x86-64/WolframKernel -script /var/projects/toy_physics/_scratch/s15/reading_build_work/rep
```

```
$ systemctl --user list-units --no-legend 's11c-guard-*' | cut -c1-120
  s11c-guard-b4d073f2bd93.service loaded active running /usr/bin/python3 /var/projects/toy_physics/scripts/s11c_guarded_
```

```
$ python3 -c "import json;r=json.load(open('/home/trevnorris/.local/state/s11c-guard/resource-pool/reservations.json'));[print(k,v['memoryMax'],v['owner']['pid'],v['logDirectory'].rsplit('/',1)[1]) for k,v in r.items()]"
s11c-guard-b4d073f2bd93.service 8589934592 2037429 repair4_decision_final_wolfram_01
```

```
$ ls -l --time-style=+%Y-%m-%dT%H:%M:%SZ _scratch/s11c/s15/repair4_decision_final_wolfram_01/; test -e _scratch/s11c/s15/repair4_decision_final_wolfram_01/outcome.json && echo 'outcome.json present' || echo 'outcome.json absent'
total 6744
-rw-rw-r-- 1 trevnorris trevnorris     164 2026-11-06T15:40:00Z duration-validation.json
-rw-rw-r-- 1 trevnorris trevnorris     594 2026-11-06T15:40:00Z effective-limits.json
-rw-rw-r-- 1 trevnorris trevnorris    1446 2026-11-06T15:40:00Z invocation.json
-rw-rw-r-- 1 trevnorris trevnorris     676 2026-11-06T15:40:00Z limit-validation.json
-rw-rw-r-- 1 trevnorris trevnorris 1510640 2026-11-06T17:44:59Z resource-samples.jsonl
-rw-rw-r-- 1 trevnorris trevnorris       0 2026-11-06T15:40:00Z stderr
-rw-rw-r-- 1 trevnorris trevnorris 5367900 2026-11-06T17:44:59Z stdout
outcome.json absent
```

Run names the note gives, as files in the work directory.

```
$ ls _scratch/s15/reading_build_work/repair3/final_decision | /usr/bin/grep -E '^(run_parallel|run_sympy02|wait_event)\.py$|repair4_decision_final_(sympy_0[12]|wolfram_01)\.(py|wl)$'
repair4_decision_final_sympy_01.py
repair4_decision_final_sympy_02.py
repair4_decision_final_wolfram_01.wl
run_parallel.py
run_sympy02.py
wait_event.py
```

Production sources (unchanged by the orchestrator).

```
$ sha256sum research/pde_ledger_v3/scripts/S15_separation_departure.py research/pde_ledger_v3/mathematica/S15_separation_departure.wl research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
fcbb23bcd27aa3e3d899318c3e1f8333ccf895b93468006ec0935d8baea7f204  research/pde_ledger_v3/scripts/S15_separation_departure.py
366298d82e62d1f36a4ec715feeaa2a6a37903b69a9853d83a4fde50089f19d2  research/pde_ledger_v3/mathematica/S15_separation_departure.wl
ca79558aedf0fcfbb3c3c89e79fa48e5cb40767433947c121eed4d427bcce6b4  research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```

```
$ git status --short -- research/pde_ledger_v3
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
 T research/pde_ledger_v3/scripts/out/S15_separation_departure.out
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair4o_note_lookups.md
?? research/pde_ledger_v3/directives/_measurements/S15_separation_departure_attempt_history.json
?? research/pde_ledger_v3/directives/_measurements/S15_separation_departure_attempt_history.md
?? research/pde_ledger_v3/directives/resume_s15_reading_repair_4o.md
?? research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
```
