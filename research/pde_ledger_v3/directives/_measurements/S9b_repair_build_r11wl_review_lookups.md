# Measurements: S9b build review, round 11, Wolfram (generated 2026-10-19 00:48)

Generator: `_scratch/s9b_build/gen/s9b_build_r11wl_lookups.sh` (sha256sum/grep/sed/head/wc/ls only). The legs'
evidence is hard-linked at `_scratch/s9b_build/s9b_repair_build_review_r11wl_{claude,grok}_evidence/`.

## The reviewed files and the prompt
```
$ (cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r11wl.sha256 | grep -c ': OK$')
30
```

```
$ wc -l < /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_baseline_r11wl.sha256
30
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_prompt_r11wl.md | cut -c1-64
95bd53609afe43b459fdd8f811a56d68582480426f6434fc3f1195b100b93435
```

```
$ git -C /var/projects/toy_physics diff --quiet 86097a44 -- /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl && echo engine-unchanged-since-86097a44
engine-unchanged-since-86097a44
```

## Verdicts
```
$ head -1 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_claude.md
**Verdict: NEEDS REVISION.** The only finding is about the record, and fixing it needs no rerun. The K10 repair, the harness, the checks script and every printed output are confirmed by my own computation. Nothing else moved.
```

```
$ grep -o 'Verdict:\*\* [A-Z][A-Z ]*' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_grok.txt
Verdict:** CLEAR
```

## K10's mutation: both legs
```
$ grep -n 'K10' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | head -3
29: "K10" -> {"rayGate = propagationGate && traversalGate;", "rayGate = With[{domain = traversedRadii}, propagationGate /.\n  HoldPattern[Inactive[ForAll][r, domain, predicate_]] :>\n    Not[Inactive[ForAll][r, domain, predicate]]] && traversalGate;"},
195:(* These three searches evaluate only parts reached by K6, K7 or K10.
283:(* K10: decide the observable's post-negation domain. The flat profiles are
```

```
$ sed -n '262,266p' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl
propagationGate = rayConstruction && Inactive[ForAll][r, traversedRadii,
  physicalPropagation && speedChart && rhoBr[r] != 0 &&
  Element[{velocityProfile[r],xiProfile[r],muPerp[r],rhoBr[r]},Reals]];
traversalGate = rayConstruction && Inactive[ForAll][r, traversedRadii, physicalTraversal];
rayGate = propagationGate && traversalGate;
```

```
$ grep -h 'DOMAIN_DECISION\|WHOLE_OBJECTS_SAMEQ\|PATH_TRAVERSAL_AT_INPUT\|CORRUPTED_NOT_SUBTERMS\|BRANCH_EXISTENCE_SAMEQ' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_claude_evidence/out/k10_witness_check_v3.stdout | head -8
CORRUPTED_NOT_SUBTERMS: {fa15 >= rTurnFlyby || rTurnRadar <= fa15 <= Max[Sqrt[b^2 + ZE^2], Sqrt[b^2 + ZR^2]]}
DOMAIN_DECISION corrupted: False
DOMAIN_DECISION baseline: True
WHOLE_OBJECTS_SAMEQ: False
PATH_TRAVERSAL_AT_INPUT corrupted: True
PATH_TRAVERSAL_AT_INPUT baseline: True
BRANCH_EXISTENCE_SAMEQ_ACROSS_COPIES: True
```

```
$ grep -h 'CORRUPTED_WHOLE\|BASELINE_WHOLE\|MATCH_PRINTED\|OBJECTS_DIFFER\|EXISTENCE_FORALL' /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_grok_evidence/k10_independent_eval_stdout.txt
CORRUPTED_WHOLE: False
BASELINE_WHOLE: True
MATCH_PRINTED: True
OBJECTS_DIFFER: True
TRAVERSAL_MATCH_PRINTED: True
EXISTENCE_FORALL: True
NEGATED_EXISTENCE_FORALL: False
```

## Form ablations
```
$ head -12 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_claude_evidence/out/form_vs_control_attempt.diff
2c2
< "Mutation" -> "rayGate = With[{domain = traversedRadii}, propagationGate /.\n  HoldPattern[Inactive[ForAll][r, domain, predicate_]] :>\n    Not[Inactive[ForAll][r, domain, predicate]]] && traversalGate;",
---
> "Mutation" -> "rayGate = With[{domain = traversedRadii}, propagationGate /.\n  HoldPattern[Inactive[ForAll][r, domain, predicate_]] :>\n    (Inactive[ForAll][r, domain, predicate] && Inactive[ForAll][r, domain, r > 0])] && traversalGate;",
5c5
< WL_LOCAL_S9B_ABLATION_K10_ATTEMPT: <|"Result" -> "witness",
---
> WL_LOCAL_S9B_ABLATION_K10_ATTEMPT: <|"Result" -> "undecided",
12c12
< "Source" -> "/var/projects/_review_scratch/s9b-r11wl-claude/rootC/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl",
---
> "Source" -> "/var/projects/_review_scratch/s9b-r11wl-claude/rootF/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl",
```

```
$ head -12 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_grok_evidence/form_ablation_diff_stdout.txt
SAME_AFTER_PATH_NORM 5336
DIFFERING 347 EXTRA_R11 False EXTRA_FORM False
BY_KNIFE Counter({'K10': 320, 'K9a': 3, 'K9b': 3, 'K1': 2, 'K2': 2, 'K3': 2, 'K4a': 2, 'K4b': 2, 'K4c': 2, 'K5': 2, 'K6': 2, 'K8': 2, 'K11': 2, 'MANIFEST': 1})
HUNK WL_LOCAL_S9B_ABLATION_MANIFEST line 1 at 1351 lens 1488 1488
  R11  'te /.\\n  HoldPattern[Inactive[ForAll][r, domain, predicate_]] :>\\n    Not[Inactive[ForAll][r, domain, predicate]]] && traversalGate;"}, "K11" -> {"vVector = {vL'
  FORM 'te /.\\n  HoldPattern[Inactive[ForAll][r, domain, predicate_]] :>\\n    Inactive[ForAll][r, domain, Not[predicate]]] && traversalGate;"}, "K11" -> {"vVector = {vL'
HUNK WL_LOCAL_S9B_ABLATION_K1_ATTEMPT line 379 at 4016 lens 12659 12659
  R11  'ness`q}, -(10*S9bHarness`q^2)^(-1)]}, "DerivativeVariable" -> radius391, "DerivativeRange" -> radius391 > 0, "DerivativesThroughOrder4" -> <|deltaProfile -> {1/'
  FORM 'ness`q}, -(10*S9bHarness`q^2)^(-1)]}, "DerivativeVariable" -> radius392, "DerivativeRange" -> radius392 > 0, "DerivativesThroughOrder4" -> <|deltaProfile -> {1/'
HUNK WL_LOCAL_S9B_ABLATION_K1_WITNESS_SEARCH line 380 at 4023 lens 12842 12842
  R11  'ness`q}, -(10*S9bHarness`q^2)^(-1)]}, "DerivativeVariable" -> radius391, "DerivativeRange" -> radius391 > 0, "DerivativesThroughOrder4" -> <|deltaProfile -> {1/'
  FORM 'ness`q}, -(10*S9bHarness`q^2)^(-1)]}, "DerivativeVariable" -> radius392, "DerivativeRange" -> radius392 > 0, "DerivativesThroughOrder4" -> <|deltaProfile -> {1/'
```

## Non-K10 knives against round 10
```
$ head -8 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_claude_evidence/out/nonk10_diff_spans.stdout
1 K1_A 
   r10: '}, "SitePart" -> "value", "Source" -> "/var/projects/toy_physics-s9b-r8-wl/r' 
   r11: '}, "SitePart" -> "value", "Source" -> "/var/projects/toy_physics-s9b-r11-wl/r'
1 K1_W 
   r10: '}, "SitePart" -> "value", "Source" -> "/var/projects/toy_physics-s9b-r8-wl/r' 
   r11: '}, "SitePart" -> "value", "Source" -> "/var/projects/toy_physics-s9b-r11-wl/r'
1 K2_A 
   r10: '}, "SitePart" -> "value", "Source" -> "/var/projects/toy_physics-s9b-r8-wl/r' 
```

```
$ head -8 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_claude_evidence/out/k10_r10_r11.stdout
K10 per-tag lines: r10 374 r11 374 Baseline operand identical across rounds: 374 differing: []
r10 marks {'"identical"': 54, '"not compared"': 320}
r11 marks {'"identical"': 57, '"not compared"': 317}
not-compared only in r10: ['WL_S9B_BRANCH_EXISTENCE', 'WL_S9B_F_FLOW_BRANCH_EXISTENCE', 'WL_S9B_F_LIVE_BRANCH_EXISTENCE']
not-compared only in r11: []
```

```
$ head -12 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_grok_evidence/compare_rounds_stdout.txt
LINES_R11 5683
KNIFE_LINE_COUNTS_R11 [('HARNESS', 374), ('K1', 379), ('K10', 379), ('K11', 379), ('K2', 379), ('K3', 379), ('K4a', 379), ('K4b', 379), ('K4c', 379), ('K5', 379), ('K6', 379), ('K7', 379), ('K8', 379), ('K9a', 380), ('K9b', 380), ('MANIFEST', 1)]
KNIFE_LINE_COUNTS_R10 [('HARNESS', 374), ('K1', 379), ('K10', 379), ('K11', 379), ('K2', 379), ('K3', 379), ('K4a', 379), ('K4b', 379), ('K4c', 379), ('K5', 379), ('K6', 379), ('K7', 379), ('K8', 379), ('K9a', 380), ('K9b', 380), ('MANIFEST', 1)]
TAG_MISMATCHES 0
NON_K10_DIFFERING_LINES 29
  (1, 'WL_LOCAL_S9B_ABLATION_MANIFEST', 1488, 1667)
  (379, 'WL_LOCAL_S9B_ABLATION_K1_ATTEMPT', 12687, 12686)
  (380, 'WL_LOCAL_S9B_ABLATION_K1_WITNESS_SEARCH', 12870, 12869)
  (758, 'WL_LOCAL_S9B_ABLATION_K2_ATTEMPT', 12668, 12667)
  (759, 'WL_LOCAL_S9B_ABLATION_K2_WITNESS_SEARCH', 12851, 12850)
  (1137, 'WL_LOCAL_S9B_ABLATION_K3_ATTEMPT', 12719, 12718)
  (1138, 'WL_LOCAL_S9B_ABLATION_K3_WITNESS_SEARCH', 12902, 12901)
```

## Claude finding 1: the superseded session's runs
```
$ sed -n '3p' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_measurements.md
This packet uses only the final round-11 sources and their current executions. The engine is unchanged. The only harness repair is K10's observable gate, site inventory and domain witness; the scratch output destination is isolated for this round. The checks exercise synthetic gate factors and conditional-object readings. Input snapshots are scope evidence, not a source of measured results. This i
```

```
$ grep -n 'Guard admission refusals' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_measurements.md /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_measurements.md:128:Guard admission refusals: **10**, all solely because available memory did not cover outstanding commitments plus this job. Admission was retried after waits of about three minutes at the unchanged 8 GiB limit; no CAS workload started in those refused attempts. Exact refusal trac
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md:3:The final-code baseline/K1–K11 harness, checks and both engine portability runs completed with exit code 0 through `--pool s9b --memory-gib 8 --tasks-max 64`, one Wolfram kernel at a time, without execution deadlines. No guard kill, memory exhaustion, task-limit or item-12 event occurre
```

```
$ grep -c -iE 'prior-code|supersed|4dc92515' /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_measurements.md /var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_round11_measurements.md:0
/var/projects/toy_physics/research/pde_ledger_v3/directives/_measurements/S9b_repair_wolfram_build_report.md:0
```

```
$ sed -n '423,433p' /var/projects/toy_physics/research/pde_ledger_v3/directives/S9b_repair_build_directive.md
17. **Report** (in the final message):
    - the tags emitted;
    - each deferred item (item 13) and how it was handled;
    - each binding and its spec line;
    - each family restriction and what it leaves out;
    - each knife, and where it acts;
    - every `NOT_ESTABLISHED`;
    - every guard refusal or kill;
    - every stop-and-report event (item 12);
    - for K11, the compact method, its coverage, and every tag printed as not evaluated;
    - every shared-vocabulary name (item 4) the engine does not emit, and why.
```

```
$ grep -n 'superseded' /var/projects/toy_physics/_scratch/s9b_build/repair_build/brief_wl_r11.md
16:3. Regenerate your measurement files from those runs, and rewrite the report as item 17 asks, so that nothing in them describes a superseded run. Then stop. Commit nothing.
```

```
$ grep -o '"event": "guard_[a-z]*", "time": [0-9.]*, "name": "[a-z0-9-]*"' /var/projects/toy_physics/_scratch/s9b_build/r11_harness_output/wl/builder_scratch/driver-events.jsonl
"event": "guard_start", "time": 1792369544.6603575, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792369545.1610677, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792369725.1837928, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792369725.6844957, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792369905.705272, "name": "s9b-r11-checks-final"
"event": "guard_complete", "time": 1792369907.7202291, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792369907.733723, "name": "s9b-r11-full"
"event": "guard_complete", "time": 1792370214.1520822, "name": "s9b-r11-full"
"event": "guard_start", "time": 1792370214.1658719, "name": "s9b-r11-engine-repo"
"event": "guard_complete", "time": 1792370230.4138472, "name": "s9b-r11-engine-repo"
"event": "guard_start", "time": 1792370230.427149, "name": "s9b-r11-engine-alone"
"event": "guard_complete", "time": 1792370246.6887434, "name": "s9b-r11-engine-alone"
"event": "guard_start", "time": 1792370405.2887216, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792370405.7894845, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792370585.8116033, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792370586.3123384, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792370766.3346782, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792370766.8354545, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792370946.8574672, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792370947.3581972, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792371127.3803427, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792371127.8810875, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792371307.9033103, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792371308.4040585, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792371488.4261498, "name": "s9b-r11-checks-final"
"event": "guard_refusal", "time": 1792371488.9275527, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792371668.949198, "name": "s9b-r11-checks-final"
"event": "guard_complete", "time": 1792371670.9505432, "name": "s9b-r11-checks-final"
"event": "guard_start", "time": 1792371680.459241, "name": "s9b-r11-full"
"event": "guard_refusal", "time": 1792371680.960049, "name": "s9b-r11-full"
"event": "guard_start", "time": 1792371860.9825902, "name": "s9b-r11-full"
"event": "guard_complete", "time": 1792372177.9816213, "name": "s9b-r11-full"
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s9b_build/r11_harness_output/wl/builder_scratch/prior-code/harness-pinned-final.wl /var/projects/toy_physics/_scratch/s9b_build/r11_harness_output/wl/builder_scratch/harness-pinned-final.wl /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl | cut -c1-16,66-
4dc92515b3066a39 /var/projects/toy_physics/_scratch/s9b_build/r11_harness_output/wl/builder_scratch/prior-code/harness-pinned-final.wl
76659e389e9e17de /var/projects/toy_physics/_scratch/s9b_build/r11_harness_output/wl/builder_scratch/harness-pinned-final.wl
76659e389e9e17de /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_ablation.wl
```

```
$ cat /var/projects/toy_physics/_scratch/s9b_build/r11_harness_output/wl/builder_scratch/launcher.json
{
  "jobPid": 195465,
  "watcherPid": 195466
}
```

```
$ ls /var/projects/toy_physics/../toy_physics-s9b-r11-wl/_scratch/s11c/
resource-guard.lock
s9b-r11-checks-final
s9b-r11-checks-final-memory-refusal-0
s9b-r11-checks-final-memory-refusal-1
s9b-r11-checks-final-memory-refusal-2
s9b-r11-checks-final-memory-refusal-3
s9b-r11-checks-final-memory-refusal-4
s9b-r11-checks-final-memory-refusal-5
s9b-r11-checks-final-memory-refusal-6
s9b-r11-checks-final-prior-code
s9b-r11-engine-alone
s9b-r11-engine-repo
s9b-r11-full
s9b-r11-full-memory-refusal-0
s9b-r11-full-prior-code
```

```
$ head -12 /var/projects/toy_physics/_scratch/s9b_build/s9b_repair_build_review_r11wl_claude_evidence/out/superseded_run_evidence.stdout
SHA256 prior-code/harness-pinned-final.wl 4dc92515b3066a3926b524e419260693826849b12ecc2622eb0aec5db9098914
SHA256 prior-code/checks-pinned-final.wl 5688ac1cd41c902ef3d5c4757cebf87aa192491e7ca647e564087019c39e99b6
SHA256 harness-pinned-final.wl 76659e389e9e17ded980c1b704279072bcb7d57d90e4d1db2d46780dfb40f008
SHA256 checks-pinned-final.wl 1fd77df88c1e8065fc0ec4ec225db02d3ef254883ab82719c3d11c8e875f9516
SHA256 engine-pinned-final.wl 2028badd3772d68d43659e06d6ebc6644843dd0fb77563fc6fff37255f241561
LAUNCHER {   "jobPid": 195465,   "watcherPid": 195466 }
EVENT 17:31:28 waiting_for_kernel    0
EVENT 18:31:47 guard_complete s9b-r11-checks-final s11c-guard-9bfe4f206cb8 0 2
EVENT 18:36:54 guard_complete s9b-r11-full s11c-guard-d420f0fc3eba 0 306
EVENT 18:37:10 guard_complete s9b-r11-engine-repo s11c-guard-4d39d728ca3c 0 16
EVENT 18:37:26 guard_complete s9b-r11-engine-alone s11c-guard-489323f085ab 0 16
EVENT 18:37:26 complete    0
```

