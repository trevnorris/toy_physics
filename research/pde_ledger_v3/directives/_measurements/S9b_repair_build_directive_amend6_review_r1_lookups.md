# Measurements — S9b directive amendment 6, review round 1 (generated 2026-10-10 11:58)

Generator: `_scratch/s9b_build/gen/s9b_amend6_r1_lookups.sh` (sha256sum/grep/sed/cat only).

## The reviewed version is the version committed
```
$ sha256sum research/pde_ledger_v3/directives/S9b_repair_build_directive.md _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r1.md
3072d5116210d0bfc24489133c41390b662905591ab8ec61de999801964ca261  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
3072d5116210d0bfc24489133c41390b662905591ab8ec61de999801964ca261  _scratch/s9b_build/S9b_repair_build_directive_amend6_reviewed_r1.md
```

```
$ head -1 _scratch/s9b_build/s9b_repair_build_directive_amend6_review_baseline_r1.sha256
3072d5116210d0bfc24489133c41390b662905591ab8ec61de999801964ca261  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend6_review_prompt_r1.md
a0394cbd8783d861906780773df9a03816179954800b83ac366932430c42c7af  _scratch/s9b_build/s9b_repair_build_directive_amend6_review_prompt_r1.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend6_review_r1_codex_final.txt _scratch/s9b_build/s9b_amend6_review_r1_grok.txt
_scratch/s9b_build/s9b_amend6_review_r1_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_amend6_review_r1_grok.txt:1:Verdict: CLEAR
```

```
$ grep -n -o 'Findings: none[^.]*' _scratch/s9b_build/s9b_amend6_review_r1_codex_final.txt
1:Findings: none within the requested scope
```

```
$ grep -n -o 'CLEAR.\*\* No findings' _scratch/s9b_build/s9b_amend6_review_r1_grok.txt
1:CLEAR.** No findings
```

## Codex's stalled Wolfram probe: what the leg claims from it
```
$ grep -n -o 'I do not claim a completed Wolfram probe[^.]*' _scratch/s9b_build/s9b_amend6_review_r1_codex_evidence/report.md
17:I do not claim a completed Wolfram probe
```

```
$ cat _scratch/s11c/s9b-amend6-review-r1/codex-independent-wl-4seadt-01/outcome.json
{
  "exitCode": 1,
  "wallSeconds": 3972.3076907179784,
  "unit": "s11c-guard-079018730038",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}
```

## The guard registry after the stop (owner alive, released by the guard)
```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

