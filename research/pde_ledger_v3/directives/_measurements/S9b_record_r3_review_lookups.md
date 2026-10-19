# Measurements: S9b record review, round 3 (generated 2026-10-19 06:57)

Generator: `_scratch/s9b_build/gen/s9b_record_r3_lookups.sh` (sha256sum/grep/sed/cut/wc only). The legs'
evidence is hard-linked at `_scratch/s9b_build/record/s9b_record_review_r3_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_baseline_r3.sha256 | grep -c ': OK$')
4
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_prompt_r3.md | cut -c1-64
42dc76e2a9dd5052c56e14edd3c2e0d18fcefcd09aea7afb441e6d37e46d7d6c
```

## Verdicts
```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r3_claude.md
Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict: [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r3_grok.txt
Verdict: NEEDS REVISION
```

## Both legs: the generator reproduces the filed measurements
```
$ sha256sum /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md | cut -c1-64
b1c4393a4cb872f6c1b487d2fb9f827b446cc23c206a6e49c37e5222f525cddd
```

```
$ grep -h 'WHOLE_FILE_IDENTICAL' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r3_claude_evidence/regen/compare_sections.stdout | head -1
WHOLE_FILE_IDENTICAL True
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r3_grok_evidence/rerun_measurements.md | cut -c1-64
b1c4393a4cb872f6c1b487d2fb9f827b446cc23c206a6e49c37e5222f525cddd
```

## Both legs (Claude 1, Grok 1): the river/acoustic-flow check's owner and source
```
$ sed -n '766,772p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
The supplied GR/PPN comparison references do not constitute an applied oracle check.

**User decision, 2026-10-19:** a check of the flow-only condition against the river/acoustic-flow
prior art is queued for **after this record is accepted**. It is separate new computation, outside
this record; **no owner named** in that decision. This queues the check without selecting an
oracle's applicability, importing a flow profile, or answering its second-route question. The
identified river/inflow/acoustic candidates in the table remain open with that handoff.
```

```
$ sed -n '1073,1079p' /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
The supplied GR/PPN comparison references do not constitute an applied oracle check.

**User decision, 2026-10-19:** a check of the flow-only condition against the river/acoustic-flow
prior art is queued for **after this record is accepted**. It is separate new computation, outside
this record; **no owner named** in that decision. This queues the check without selecting an
oracle's applicability, importing a flow profile, or answering its second-route question. The
identified river/inflow/acoustic candidates in the table remain open with that handoff.
```

```
$ git -C /var/projects/toy_physics show 4760e6e2:research/pde_ledger_v3/directives/_measurements/S9b_record_r2_review_disposition.md | sed -n '19,21p'
- **A user decision this round.** Asked whether a check of the flow-only condition against the river/acoustic-flow
  prior art should be queued, the user chose "Queue it after record". The check is new computation, done after the
  record is accepted and outside it: Codex writes the check and it is reviewed under G1.
```

```
$ grep -c 'S9b_record_r2_review_disposition' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
/var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md:0
/var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md:0
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:0
```

```
$ grep -c 'Queue it after record' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md /var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md /var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md
/var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md:0
/var/projects/toy_physics/research/pde_ledger_v3/steps/_measurements/S9b_record_measurements.md:0
/var/projects/toy_physics/research/pde_ledger_v3/SUBSTRATE_REQUIREMENTS.md:0
```

```
$ sed -n '3,5p' /var/projects/toy_physics/_scratch/s9b_build/record/resume_record_r3.md
A user decision to record, dated 2026-10-19. A check of the flow-only condition against the river/acoustic-flow prior art is queued for after this record is accepted. It is separate new computation, outside this record.

After this round, all of these must be true:
```

## Claude 2: the rows that carry the unresolved brane-zero value samples
```
$ sed -n '559p;561p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
| A_DEFLECTION, A_ROUND_TRIP, A_ONE_WAY_ER/RE, A_NONRECIPROCAL, A_RADAR_LOG; conditional F_FLOW/F_LIVE copies | Formed value samples and markers; all twelve applicable grades retained, four radar grades | No full along-ray domain comparison; no all-profile equality; unresolved brane-zero samples and printed numeric differences remain. Shared nonreciprocal objects retain §2/M14: K11 supports local
| B_RESIDUAL_DEFLECTION/RADAR and R_*/conditional F_* copies | Formed value samples, with printed differences retained | No selected profiles or GR match, and no domain equality |
```

```
$ sed -n '628,631p;641,643p' /var/projects/toy_physics/research/pde_ledger_v3/steps/S9b_light_bending.md
- Vanishing-brane-density samples retain interior-denominator-zero/integrability or pole-set reasons.
  The r4 acceptance counts 183 and 27 occurrences of its two pole-reason families and 105 unresolved
  value samples in that stratum. The broader printed value-unresolved census remains 319; these
  differently scoped counts are not substituted for one another.
  changed their scope. In particular the 36 mixed-grade rows now form value samples on generic
  inverse/shifted/exponential profiles, with the remaining brane-zero gap, rather than remaining
  wholly untested. The old `approx_30` caveat is not a claim of 30 reliable production digits.
```

```
$ grep '^{"kind":"comparison","name":"B_RESIDUAL_DEFLECTION"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
0
```

```
$ grep '^{"kind":"comparison","name":"F_FLOW_B_RESIDUAL_DEFLECTION"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
3
```

```
$ grep '^{"kind":"comparison","name":"F_LIVE_B_RESIDUAL_RADAR"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
3
```

```
$ grep '^{"kind":"comparison","name":"R_FLOW_ONLY_B_RESIDUAL_RADAR"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
0
```

```
$ grep '^{"kind":"comparison","name":"A_DEFLECTION_G020"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
0
```

```
$ grep '^{"kind":"comparison","name":"F_FLOW_A_DEFLECTION_G020"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
3
```

```
$ grep '^{"kind":"comparison","name":"R_FLOW_ONLY_A_DEFLECTION_G020"' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S9b_cross_engine_comparator.out | grep -o 'operand evaluations unavailable for subtraction' | wc -l
0
```

```
$ grep -h '^F_\|^R_\|^general_' /var/projects/toy_physics/_scratch/s9b_build/record/s9b_record_review_r3_claude_evidence/brane_zero_by_class.stdout
F_A {'residual formed': 165, 'unresolved coverage': 93}
F_B {'unresolved coverage': 12}
F_other {'unresolved coverage': 12}
R_A {'residual formed': 90}
R_B {'residual formed': 18}
general_A {'residual formed': 129, 'boolean_difference': 3}
general_B {'residual formed': 6}
general_other {'unresolved coverage': 6}
```

## The B_RESIDUAL row-class cell at each reviewed version
```
$ git -C /var/projects/toy_physics show 36f59ff1:research/pde_ledger_v3/steps/S9b_light_bending.md | grep -c 'Formed value samples, with printed differences retained | No selected profiles or GR match, and no domain equality'
1
```

```
$ git -C /var/projects/toy_physics show d1511b20:research/pde_ledger_v3/steps/S9b_light_bending.md | grep -c 'Formed value samples, with printed differences retained | No selected profiles or GR match, and no domain equality'
1
```

```
$ git -C /var/projects/toy_physics show 4760e6e2:research/pde_ledger_v3/steps/S9b_light_bending.md | grep -c 'Formed value samples, with printed differences retained | No selected profiles or GR match, and no domain equality'
1
```

