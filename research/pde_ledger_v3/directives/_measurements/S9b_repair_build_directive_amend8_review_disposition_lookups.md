# Measurements — S9b directive amendment 8, review round 0 (generated 2026-10-11 20:55)

Generator: `_scratch/s9b_build/gen/s9b_amend8_r0_lookups.sh` (sha256sum/grep/sed/cat/ls only).

## The round-0 reviewed version
```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend8_review_baseline_r0.sha256
cf258688a4f193d61a67ac7ff2ab72cfc13d2fa07168e753d20b2b55a3f7032d  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md
cf258688a4f193d61a67ac7ff2ab72cfc13d2fa07168e753d20b2b55a3f7032d  _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r0.md
4fba4aca1e0ab26677952195d7cad72d7a5651b2458a03e45476e2ff9d947132  _scratch/s9b_build/s9b_repair_build_directive_amend8_review_prompt_r0.md
```

```
$ grep -n 'that positivity is' _scratch/s9b_build/S9b_repair_build_directive_amend8_reviewed_r0.md
165:   `ρ = ρ₀(1 + f)` is positive, and that positivity is a supplied input of each row whose response reads `ρ`, so it
```

## Both legs ran to completion, and both verdicts
```
$ cat _scratch/s9b_build/amend8_r0_codex.DONE _scratch/s9b_build/amend8_r0_grok.DONE
EXIT=0
EXIT=0
```

```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
_scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt:1:Verdict: NEEDS REVISION
_scratch/s9b_build/s9b_amend8_review_r0_grok.txt:1:Verdict: NEEDS REVISION
```

## Codex's finding, and what it says of item 14
```
$ sed -n 3p _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt
1. **Item 9, lines 164–166:** strict positivity adds a premise the spec does not supply. The independent script constructs a nonnegative density with a local zero and positive asymptotic density. At `n = 1`, the fixed-ratio response remains exactly defined and positive there, but the amendment excludes that stratum. **Fix:** p
```

```
$ grep -n -o 'No additional findings on item 14, coverage, or M2' _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt
5:No additional findings on item 14, coverage, or M2
```

## Grok's findings, and what it says of item 14
```
$ sed -n 9p _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
**1. Item 9 states a sign the spec does not write, and conjoining it changes the Part C domains.**
```

```
$ sed -n 43p _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
**2. As written, "whose response reads `ρ`" does not fix the K9 copies.**
```

```
$ grep -n -o 'This stands only while finding 1.s sign remains' _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
45:This stands only while finding 1's sign remains
```

```
$ grep -n -o '\*\*Item 14\.\*\* The new sentences match C1 and C2' _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
74:**Item 14.** The new sentences match C1 and C2
```

## Grok's script output for finding 1
```
$ grep -n 'DENSITY_GT_0\|RESPONSE_AT_NEG_DENSITY\|ZERO_DENSITY' _scratch/s9b_build/s9b_amend8_review_r0_grok_evidence/amend8_domain_predicates.stdout
6:ZERO_DENSITY_IDENTITY_RESIDUAL: 0.0
7:ZERO_DENSITY_PRODUCT: 0
8:ZERO_DENSITY_PRODUCT_GT_0: False
9:ZERO_DENSITY_RHO0_GT_0: True
10:ZERO_DENSITY_RHO_GT_0: False
21:POWER_RESPONSE_AT_NEG_DENSITY: 1
24:POWER_DENSITY_GT_0: False
27:FIXED_RESPONSE_AT_NEG_DENSITY: 1
31:FIXED_DENSITY_GT_0: False
34:FIXED_LOCAL_CS2_AT_ZERO_DENSITY: 1
35:FIXED_FAR_CS2_AT_ZERO_DENSITY: 1
36:FIXED_RESPONSE_AT_ZERO_DENSITY: 1
38:FIXED_ZERO_DENSITY_GT_0: False
39:FIXED_ZERO_DENSITY_GE_0: True
```

## The two leftover Wolfram units: reported by the legs, stopped and reconciled
```
$ grep -n -o 'supplementary Wolfram run hit the task ceiling' _scratch/s9b_build/s9b_amend8_review_r0_codex_final.txt
9:supplementary Wolfram run hit the task ceiling
```

```
$ grep -n -o 'It was not stopped' _scratch/s9b_build/s9b_amend8_review_r0_grok.txt
5:It was not stopped
```

```
$ cat _scratch/s11c/s9b-amend8-review-r0/independent-wolfram-FqQ9xRVa/RECONCILIATION.md
# Reconciliation (orchestrator, 2026-10-11 20:54)

Unit `s11c-guard-462fd9341920` was a supplementary Wolfram run by the Codex leg of the amendment 8 review, round 0. The leg exited (its
`.DONE` marker exists) and reported the run as stuck at the 32-task ceiling, not used as evidence. The unit's owner,
pid 2943850, was dead. Nothing was queued. I stopped the unit with `systemctl --user stop`, confirmed it inactive and no
WolframKernel left, backed the registry up to `resource-pool/reservations.json.stale-amend8-r0-legs.bak`, and wrote
`{}`.
```

```
$ cat _scratch/s11c/s9b-amend8-review-r0/amend8-k9-forms/RECONCILIATION.md
# Reconciliation (orchestrator, 2026-10-11 20:54)

Unit `s11c-guard-4e9c9677fa88` was a supplementary Wolfram run by the Grok leg of the amendment 8 review, round 0. The leg exited (its
`.DONE` marker exists) and reported the run as stuck at the 32-task ceiling, not used as evidence. The unit's owner,
pid 2952936, was dead. Nothing was queued. I stopped the unit with `systemctl --user stop`, confirmed it inactive and no
WolframKernel left, backed the registry up to `resource-pool/reservations.json.stale-amend8-r0-legs.bak`, and wrote
`{}`.
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

```
$ systemctl --user list-units 's11c-guard-*' --no-legend | wc -l
0
```

