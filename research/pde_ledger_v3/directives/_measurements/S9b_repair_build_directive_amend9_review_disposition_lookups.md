# Measurements — S9b directive amendment 9 review, round 0 (generated 2026-10-12 08:35)

Generator: `_scratch/s9b_build/gen/s9b_amend9_r0_lookups.sh` (sha256sum/grep/cat/ls only). The legs' scripts and
stdout are copied to `_scratch/s9b_build/s9b_amend9_review_r0_{codex,grok}_evidence/` (files under 2 MB).

## The reviewed text is the working-tree directive
```
$ sha256sum -c _scratch/s9b_build/s9b_repair_build_directive_amend9_review_baseline_r0.sha256
research/pde_ledger_v3/directives/S9b_repair_build_directive.md: OK
```

```
$ sha256sum _scratch/s9b_build/s9b_repair_build_directive_amend9_review_prompt_r0.md
99af0fac9080039df61fce3991a4cbad411f9a78dea66de187037f9cdd020ea2  _scratch/s9b_build/s9b_repair_build_directive_amend9_review_prompt_r0.md
```

## Verdicts
```
$ cat _scratch/s9b_build/s9b_amend9_review_r0_codex.DONE _scratch/s9b_build/s9b_amend9_review_r0_grok.DONE _scratch/s9b_build/s9b_amend9_review_r0b_codex.DONE
EXIT=0
EXIT=0
EXIT=0
```

```
$ cat _scratch/s9b_build/s9b_amend9_review_r0_codex_final.txt | head -3
Verdict withheld: the guard refused both permitted attempts with `aggregate 16 GiB reservation budget exceeded`. The independent physics script never started, so neither CLEAR nor NEEDS REVISION is substantiated.

No physics findings were established. I stopped after the required wait and single retry.
```

```
$ cat _scratch/s9b_build/s9b_amend9_review_r0b_codex_final.txt | head -3
**Verdict: CLEAR.**

**Findings:** None.
```

```
$ grep -o 'Verdict: [A-Z ]*' _scratch/s9b_build/s9b_amend9_review_r0_grok.txt | head -1
Verdict: CLEAR
```

## The Codex leg's first pass was blocked by the orchestrator's production runs
```
$ grep -c 'aggregate 16 GiB reservation budget exceeded' _scratch/s9b_build/s9b_amend9_review_r0_codex_log.txt
8
```

```
$ grep -h 'startedUtc' _scratch/s11c/s9b-production/py-audit/invocation.json _scratch/s11c/s9b-production/wl-audit/invocation.json
  "startedUtc": "2026-10-12T13:52:48.854838+00:00",
  "startedUtc": "2026-10-12T13:59:52.001870+00:00",
```

```
$ cat _scratch/s9b_build/s9b_amend9_review_r0_codex_resume.md
Status: the guard refusals you hit came from an orchestrator job (the S9b production runs) that held the aggregate 16 GiB reservation budget at the time. That job has ended, the registry is now `{}`, and no guard unit is active. This is not a retry of a failed physics run; your physics script never 
```

```
$ ls _scratch/s9b_build/s9b_amend9_review_r0_codex_evidence | grep -E 'physics_contract|wolfram_contract'
physics_contract.command
physics_contract.guard.stderr
physics_contract.guard.stdout
physics_contract.py
physics_contract.retry.exit
physics_contract.retry.guard.stderr
physics_contract.retry.guard.stdout
physics_contract.run03.exit
physics_contract.run03.guard.stderr
physics_contract.run03.guard.stdout
physics_contract.stderr
physics_contract.stdout
wolfram_contract.command
wolfram_contract.exit
wolfram_contract.guard.stderr
wolfram_contract.guard.stdout
wolfram_contract.stderr
wolfram_contract.stdout
wolfram_contract_v2.command
wolfram_contract_v2.exit
wolfram_contract_v2.guard.stderr
wolfram_contract_v2.guard.stdout
wolfram_contract_v2.stderr
wolfram_contract_v2.stdout
wolfram_contract_v2.wl
wolfram_contract.wl
```

## Guard state after both legs
```
$ cat _scratch/s11c/resource-pool/reservations.json
{}
```

