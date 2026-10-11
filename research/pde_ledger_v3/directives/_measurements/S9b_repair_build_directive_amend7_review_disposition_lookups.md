# Measurements — S9b directive amendment 7, review round 0 (generated 2026-10-10 19:44)

Generator: `_scratch/s9b_build/gen/s9b_amend7_r0_lookups.sh` (sha256sum/grep/sed/cat only).

## The reviewed version is the version committed
```
$ sha256sum research/pde_ledger_v3/directives/S9b_repair_build_directive.md _scratch/s9b_build/S9b_repair_build_directive_amend7_reviewed_r0.md
fd0d8931e9b8468e3a35357c23a78407a6d8a8d844e31ece3ac4ae9ec6d48516  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
fd0d8931e9b8468e3a35357c23a78407a6d8a8d844e31ece3ac4ae9ec6d48516  _scratch/s9b_build/S9b_repair_build_directive_amend7_reviewed_r0.md
```

```
$ cat _scratch/s9b_build/s9b_repair_build_directive_amend7_review_baseline_r0.sha256
fd0d8931e9b8468e3a35357c23a78407a6d8a8d844e31ece3ac4ae9ec6d48516  research/pde_ledger_v3/directives/S9b_repair_build_directive.md
f916eb0459dfb82d6a5087cc4d5a12a8c00f7bbe73ae311a6054a35c55e79cff  research/pde_ledger_v3/directives/_measurements/S9b_repair_build_r3_review_disposition.md
f4c9d3516a50c091a56ef4dd09366b9213b0c3a38e463daf0f43dbbaf550daaa  research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend7_review_prompt_r0.md
f090caad9131e41336d4273c357cc5c77bc11b0776060d22d2cd5ce07e0c33fc  _scratch/s9b_build/s9b_repair_build_directive_amend7_review_prompt_r0.md
```

## Both verdicts
```
$ grep -n -o 'Verdict: [A-Z][A-Z ]*' _scratch/s9b_build/s9b_amend7_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend7_review_r0_grok.txt
_scratch/s9b_build/s9b_amend7_review_r0_codex_final.txt:1:Verdict: CLEAR
_scratch/s9b_build/s9b_amend7_review_r0_grok.txt:1:Verdict: CLEAR
```

```
$ grep -n -o 'Findings: none' _scratch/s9b_build/s9b_amend7_review_r0_codex_final.txt _scratch/s9b_build/s9b_amend7_review_r0_codex_evidence/report.md
_scratch/s9b_build/s9b_amend7_review_r0_codex_final.txt:1:Findings: none
_scratch/s9b_build/s9b_amend7_review_r0_codex_evidence/report.md:3:Findings: none
```

```
$ grep -n -o 'Findings:\*\* none' _scratch/s9b_build/s9b_amend7_review_r0_grok.txt
3:Findings:** none
```

## Codex's note on the repair route
```
$ grep -n -o 'C1 nevertheless already violated[^.]*' _scratch/s9b_build/s9b_amend7_review_r0_codex_evidence/report.md
13:C1 nevertheless already violated the governing spec’s “everywhere along the ray” requirement; the route should not be read as removing the builder’s existing spec obligation
```

## Codex's stuck Wolfram diagnostic: what the leg claims from it, and the stop
```
$ grep -n -o 'Run exception: the additional native Wolfram diagnostic is incomplete and contributes no physics evidence to this verdict' _scratch/s9b_build/s9b_amend7_review_r0_codex_evidence/report.md
21:Run exception: the additional native Wolfram diagnostic is incomplete and contributes no physics evidence to this verdict
```

```
$ tail -3 _scratch/s11c/s9b-amend7-review-r0/independent-wolfram-mylKeM/stdout
ClearAll::ssym: Global<>FromCharacterCode[96]<>* is not a symbol or a valid string pattern.
radial_dispersion: -((chi*kr^2)/a) + (omega - kr*v)^2
radial_discriminant_over_4kr2: chi/a
```

```
$ cat _scratch/s11c/s9b-amend7-review-r0/independent-wolfram-mylKeM/RECONCILIATION.md
# Reconciliation (orchestrator, 2026-10-10)

The amendment-7 Codex review leg finished and reported this Wolfram diagnostic as stuck at its task cap
(pids.current 32 = pids.max 32), excluded from its verdict, and left running. I stopped the unit
`s11c-guard-db76ab8a59db` with `systemctl --user stop`. Afterwards: the unit was inactive, its owner pid 1003050
was dead, no Wolfram kernel was running and no guard unit was queued, but the registry still held the reservation.
I backed it up to `resource-pool/reservations.json.stale-amend7-codex-wl-diag.bak` and wrote `{}`.
```

```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

## Grok's guarded run
```
$ cat _scratch/s11c/s9b-amend7-review-r0/domain-predicate-2/outcome.json
{
  "exitCode": 0,
  "wallSeconds": 8.547820543171838,
  "unit": "s11c-guard-9f6754c92d47",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": {
    "exitCode": 0,
    "guardReason": null,
    "wallSeconds": 8.470627164002508,
    "childPid": 1009034
  }
}
```

