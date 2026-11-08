# S15 reading build, note 4p: fact check and leak gate

Mechanical lookups only. Run from `/var/projects/toy_physics`. The finding and its confirmation are in
`S15_reading_build_review_r4_disposition.md` and `S15_reading_build_review_r4_lookups.md` (L4–L9).

## F1. The note

```text
$ sha256sum research/pde_ledger_v3/directives/resume_s15_reading_repair_4p.md; wc -w research/pde_ledger_v3/directives/resume_s15_reading_repair_4p.md
df9d4bb35f0d5f7f2799795e4450618e653d068cf91d33f48a7f88af7dd136d1  research/pde_ledger_v3/directives/resume_s15_reading_repair_4p.md
448 research/pde_ledger_v3/directives/resume_s15_reading_repair_4p.md
```

## F2. The reviewed version's commit, and the preserved transcripts

```text
$ git log --oneline -1 0d9675eb
0d9675eb Preserve the S15 reading build repair 4 as reviewed in round 4 (not accepted)
```

```text
$ sha256sum research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r4.out research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r4.out
10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6  research/pde_ledger_v3/scripts/out/preserved/S15_separation_departure.r4.out
5779d92eb7820b11fa51068526145c3f61e04b69d2f10f86d6b1af923531abce  research/pde_ledger_v3/mathematica/out/preserved/S15_separation_departure.r4.out
```

```text
$ grep '\.out$' research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r4_baseline.sha256
10160f89704d910139b0401fbbcc1419078e45e915a2f1e663f24807647182d6  /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_separation_departure.out
5779d92eb7820b11fa51068526145c3f61e04b69d2f10f86d6b1af923531abce  /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
```

## F3. Each quotation in the note, found in its source (count of matching lines)

```text
$ grep -c 'The retained sufficient conditions are printed above.' research/pde_ledger_v3/mathematica/S15_separation_departure_run.py
1
```

```text
$ grep -c 'were not attempted in these final segments; each retained its previously printed sufficient conditions and corresponding Q objects' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md
1
```

```text
$ grep -c 'Every condition printed now that is true stays printed.' research/pde_ledger_v3/directives/resume_s15_reading_repair_4.md
1
```

```text
$ grep -c 'are not attempted again in this build' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

```text
$ grep -c 'the object was not attempted in that run' research/pde_ledger_v3/directives/resume_s15_reading_repair_4n.md
1
```

## F4. Leak gate: symbols of the item-4 objects, relations and values in the note (every hit is read below)

```text
$ grep -n -E 'aa\[|bb\[|AA\[|BB\[|==|<|>|=[^=]|exponent [a-z]+ =|equals|vanish|zero|nonzero|-[0-9]|[0-9]/[0-9]|expected|should (be|equal)' research/pde_ledger_v3/directives/resume_s15_reading_repair_4p.md || echo 'no hits'
7:For each of them, your final Wolfram run started an item-4 kernel. That kernel computed the sufficient conditions and `Q` under each, and printed them before the third-kind label. Review found these conditions true.
```

Reading of the F4 hit: line 7 matches only on `-4` in "item-4", the item's name. No relation, value, symbol of an
item-4 object or expected outcome appears in the note.
