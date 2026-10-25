# Measurements: S15 build, review round 1 lookups

Generator: `_scratch/s15/gen/s15_build_r1_lookups.sh` (GNU grep, sed, diff, cmp, git, sha256sum, `LC_ALL=C`). Each
block is a command and its literal output: existence, verbatim retrieval by line, literal comparison, match counts,
hashes. The legs' scripts were not re-run here; their filed stdout is retrieved.

## M0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/build_review/s15_build_baseline_r1.sha256
research/pde_ledger_v3/scripts/S15_two_drain_region.py: OK
research/pde_ledger_v3/mathematica/S15_two_drain_region.wl: OK
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out: OK
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out: OK
research/pde_ledger_v3/directives/_measurements/S15_build_run_record.md: OK
```

```
$ sha256sum _scratch/s15/build_review/s15_build_review_prompt_r1.md
b8bb2e031e532eff456f6a5bf7059400c18e2d0c130c5e947d16f865f06323be  _scratch/s15/build_review/s15_build_review_prompt_r1.md
```

```
$ cat /var/projects/_review_scratch/s15-build-r1-grok/DONE; ls -la /var/projects/_review_scratch/s15-build-r1-claude/REVIEW.md /var/projects/_review_scratch/s15-build-r1-grok/REVIEW.md | cut -c1-110
exit 0
-rw-rw-r-- 1 trevnorris trevnorris 7840 Oct 24 19:58 /var/projects/_review_scratch/s15-build-r1-claude/REVIEW.
-rw-rw-r-- 1 trevnorris trevnorris 8871 Oct 24 18:06 /var/projects/_review_scratch/s15-build-r1-grok/REVIEW.md
```

```
$ git status --porcelain | /usr/bin/grep -v -F -e 'S15_build_review_r1_lookups.md' -e 'S15_build_review_r1_disposition.md' | diff - <(sed -n '1,6p' _scratch/s15/build_review/status_before_r1.txt) && echo 'status unchanged apart from the r1 disposition and lookups files'
status unchanged apart from the r1 disposition and lookups files
```

```
$ sed -n '7p' _scratch/s15/build_review/status_before_r1.txt; git rev-parse HEAD
efe5ce73b8534b4ce6adf489ecad16bee6714987
efe5ce73b8534b4ce6adf489ecad16bee6714987
```

## M1. The verdicts

```
$ /usr/bin/grep -n -E '^\*\*Verdict' /var/projects/_review_scratch/s15-build-r1-claude/REVIEW.md /var/projects/_review_scratch/s15-build-r1-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-160
s15-build-r1-claude/REVIEW.md:3:**Verdict: CLEAR.** There are no findings. H1's "What must be true" now holds in both engines, and the repair changed
s15-build-r1-grok/REVIEW.md:3:**Verdict: CLEAR**
```

## M2. The repair runs that wrote the transcripts

```
$ for d in sympy_h1_01 wolfram_h1_01; do python3 -c "import json,sys;o=json.load(open('_scratch/s11c/s15/'+sys.argv[1]+'/outcome.json'));c=o.get('childOutcome') or {};print(sys.argv[1],'exitCode',o.get('exitCode'),'limitsVerified',o.get('limitsVerified'),'child',c.get('exitCode'),'guardReason',c.get('guardReason'),'stderrBytes',o.get('stderrBytes'),'wallSeconds',o.get('wallSeconds'))" $d; done
sympy_h1_01 exitCode 0 limitsVerified True child 0 guardReason None stderrBytes 0 wallSeconds 387.798950522003
wolfram_h1_01 exitCode 0 limitsVerified True child 0 guardReason None stderrBytes 0 wallSeconds 10.292321781016653
```

```
$ sha256sum research/pde_ledger_v3/scripts/out/S15_two_drain_region.out research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out
4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c  research/pde_ledger_v3/scripts/out/S15_two_drain_region.out
811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6  research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out
```

## M3. H1: what each repaired engine computes and prints for the region balance

```
$ sed -n '294,299p;311p' research/pde_ledger_v3/scripts/S15_two_drain_region.py | cut -c1-200
        parts=[surface_parts[0][i],surface_parts[1][i],forcing]
        rhs=sum(obj for surface in surface_parts for obj in surface[i].values())+forcing
        emit(tag+'/balance_parts',parts)
        pair(tag+'/balance',rate,rhs)
        emit(tag+'/balance/reduction', 'No route in this engine reduces this balance residual with the supplied local balance.')
        emit(tag+'/balance_equation',S.Eq(rate,rhs,evaluate=False))
    emit('balance_construction','Printed objects: the moving-domain momentum integral, its time derivative, the surface parts, the integrated force density, and the balance operands, equation and resi
```

```
$ sed -n '169,172p;177p' research/pde_ledger_v3/mathematica/S15_two_drain_region.wl | cut -c1-200
 rhs=Total[Flatten[Table[Values[surfaceResults[[j,i]]],{j,2}]]]+forcing;
 emit[tag<>"/balance_parts",{surfaceResults[[1,i]],surfaceResults[[2,i]],forcing}];
 pair[tag<>"/balance",lhs,rhs];emit[tag<>"/balance_equation",Equal[lhs,rhs]];
 emit[tag<>"/balance/reduction", "This balance residual has no reduction route in this engine using the supplied local balance."];
emit["balance_construction", "Printed objects: the chart momentum integral, its Leibniz time derivative, each surface contribution, the undivided force integral, and the balance operands, equation and
```

```
$ /usr/bin/grep -n -E '^(balance_construction|region/component_[123]/balance/reduction)' research/pde_ledger_v3/scripts/out/S15_two_drain_region.out research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out | cut -c1-230
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out:1294:region/component_1/balance/reduction: No route in this engine reduces this balance residual with the supplied local balance.
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out:1349:region/component_2/balance/reduction: No route in this engine reduces this balance residual with the supplied local balance.
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out:1404:region/component_3/balance/reduction: No route in this engine reduces this balance residual with the supplied local balance.
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out:2244:balance_construction: Printed objects: the moving-domain momentum integral, its time derivative, the surface parts, the integrated force density, and the balance ope
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out:1254:region/component_1/balance/reduction: "This balance residual has no reduction route in this engine using the supplied local balance."
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out:1308:region/component_2/balance/reduction: "This balance residual has no reduction route in this engine using the supplied local balance."
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out:1362:region/component_3/balance/reduction: "This balance residual has no reduction route in this engine using the supplied local balance."
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out:1368:balance_construction: "Printed objects: the chart momentum integral, its Leibniz time derivative, each surface contribution, the undivided force integral, and th
```

```
$ for f in _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.py.out research/pde_ledger_v3/scripts/out/S15_two_drain_region.out; do echo "$f: $(/usr/bin/grep -c -F 'supplied local balance produces' $f)"; done; for f in _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.wl.out research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out; do echo "$f: $(/usr/bin/grep -c -F 'Integrating the supplied local balance' $f)"; done
_scratch/s15/build_review/r0_transcripts/S15_two_drain_region.py.out: 1
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out: 0
_scratch/s15/build_review/r0_transcripts/S15_two_drain_region.wl.out: 1
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out: 0
```

## M4. Round 0 against round 1: the lines that differ (literal diff)

```
$ sha256sum _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.py.out _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.wl.out
85d99aed185464f9a9e68b2d0cc7b64193c073fdc5c321c0c60825069b2e8259  _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.py.out
8774462d38a926329ce15f524001972b9af2a5b0225ec19a2f48b7c3b2982fd1  _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.wl.out
```

```
$ for p in 'py _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.py.out research/pde_ledger_v3/scripts/out/S15_two_drain_region.out' 'wl _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.wl.out research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out'; do set -- $p; echo "$1 removed: $(diff $2 $3 | /usr/bin/grep -c '^<')  added: $(diff $2 $3 | /usr/bin/grep -c '^>')"; diff $2 $3 | /usr/bin/grep -E '^[<>]' | cut -c1-120; done
py removed: 1  added: 4
> region/component_1/balance/reduction: No route in this engine reduces this balance residual with the supplied local ba
> region/component_2/balance/reduction: No route in this engine reduces this balance residual with the supplied local ba
> region/component_3/balance/reduction: No route in this engine reduces this balance residual with the supplied local ba
< balance_construction: The rate differentiates the moving-domain momentum integral. Applying divergence and moving-doma
> balance_construction: Printed objects: the moving-domain momentum integral, its time derivative, the surface parts, th
wl removed: 1  added: 4
> region/component_1/balance/reduction: "This balance residual has no reduction route in this engine using the supplied 
> region/component_2/balance/reduction: "This balance residual has no reduction route in this engine using the supplied 
> region/component_3/balance/reduction: "This balance residual has no reduction route in this engine using the supplied 
< balance_construction: "Leibniz differentiation of the moving spherical-chart momentum integral supplies the rate. Inte
> balance_construction: "Printed objects: the chart momentum integral, its Leibniz time derivative, each surface contrib
```

```
$ cmp -s <(/usr/bin/grep -v -E '^(balance_construction|region/component_[123]/balance/reduction)' _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.py.out) <(/usr/bin/grep -v -E '^(balance_construction|region/component_[123]/balance/reduction)' research/pde_ledger_v3/scripts/out/S15_two_drain_region.out) && echo 'py: all other lines identical' || echo 'py: other lines differ'; cmp -s <(/usr/bin/grep -v -E '^(balance_construction|region/component_[123]/balance/reduction)' _scratch/s15/build_review/r0_transcripts/S15_two_drain_region.wl.out) <(/usr/bin/grep -v -E '^(balance_construction|region/component_[123]/balance/reduction)' research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out) && echo 'wl: all other lines identical' || echo 'wl: other lines differ'
py: all other lines identical
wl: all other lines identical
```

## M5. The legs' H1 checks and repair-undone controls (their stdout)

```
$ /usr/bin/grep -n -E 'CHECK|claim_lines|reduction_lines' /var/projects/_review_scratch/s15-build-r1-grok/logs/check_h1.out | cut -c1-200
2:claim_lines: 0
3:reduction_lines: 3
7:CHECK: PASS
9:claim_lines: 0
10:reduction_lines: 3
14:CHECK: PASS
16:claim_lines: 1
18:reduction_lines: 0
19:CHECK: FAIL
21:claim_lines: 1
23:reduction_lines: 0
24:CHECK: FAIL
```

```
$ /usr/bin/grep -n -E 'UNDONE' /var/projects/_review_scratch/s15-build-r1-grok/logs/ablation_summary.out | cut -c1-200
293:UNDONE_CLAIM_EQUAL_R0 True
294:UNDONE_CLAIM balance_construction: "Leibniz differentiation of the moving spherical-chart momentum integral supplies the rate. Integrating the supplied local balance and applying divergence and mo
296:UNDONE_RESIDUAL_EQUAL_R0 False
```

```
$ /usr/bin/grep -n -E 'H1_CHECK|closure residual' /var/projects/_review_scratch/s15-build-r1-claude/item0/h1_check.stdout | cut -c1-200
2:closure residual C_1: baseline=-4.809e-13 ablated(X2 motion dropped)=3.799e-01
3:closure residual C_2: baseline=1.283e-13 ablated(X2 motion dropped)=-1.444e-01
4:closure residual C_3: baseline=-2.864e-13 ablated(X2 motion dropped)=1.579e-01
7:H1_CHECK FAILS (sympy_r0_repair_undone)
9:closure residual C_1: baseline=-4.809e-13 ablated(X2 motion dropped)=3.799e-01
10:closure residual C_2: baseline=1.283e-13 ablated(X2 motion dropped)=-1.444e-01
11:closure residual C_3: baseline=-2.864e-13 ablated(X2 motion dropped)=1.579e-01
17:H1_CHECK passes (sympy_r1_repaired)
19:closure residual C_1: baseline=-4.663e-13 ablated(X2 motion dropped)=3.799e-01
20:closure residual C_2: baseline=1.317e-13 ablated(X2 motion dropped)=-1.444e-01
21:closure residual C_3: baseline=-2.849e-13 ablated(X2 motion dropped)=1.579e-01
24:H1_CHECK FAILS (wolfram_r0_repair_undone)
26:closure residual C_1: baseline=-4.663e-13 ablated(X2 motion dropped)=3.799e-01
27:closure residual C_2: baseline=1.317e-13 ablated(X2 motion dropped)=-1.444e-01
28:closure residual C_3: baseline=-2.849e-13 ablated(X2 motion dropped)=1.579e-01
34:H1_CHECK passes (wolfram_r1_repaired)
36:closure residual C_1: baseline=-4.809e-13 ablated(supplied balance corrupted: Pi = rho V V + T + Pi_o)=-4.051e+00
37:closure residual C_2: baseline=1.283e-13 ablated(supplied balance corrupted: Pi = rho V V + T + Pi_o)=-7.457e-01
38:closure residual C_3: baseline=-2.864e-13 ablated(supplied balance corrupted: Pi = rho V V + T + Pi_o)=-2.065e+00
44:H1_CHECK passes (sympy_r1_corrupt)
46:closure residual C_1: baseline=-4.663e-13 ablated(supplied balance corrupted: Pi = rho V V + T + Pi_o)=-4.051e+00
47:closure residual C_2: baseline=1.317e-13 ablated(supplied balance corrupted: Pi = rho V V + T + Pi_o)=-7.457e-01
48:closure residual C_3: baseline=-2.849e-13 ablated(supplied balance corrupted: Pi = rho V V + T + Pi_o)=-2.065e+00
54:H1_CHECK passes (wolfram_r1_corrupt)
```

## M6. The legs' object comparisons and reproduction (their stdout)

```
$ tail -3 /var/projects/_review_scratch/s15-build-r1-claude/compare_sympy.stdout | cut -c1-200; tail -3 /var/projects/_review_scratch/s15-build-r1-claude/compare_wolfram.stdout | cut -c1-200
  region/component_3/balance/operand_a: diff=3.844e-13 (engine=6.290291e-02, reviewer=6.290291e-02)
  region/component_3/rate_labels_drain1+remainder: diff=2.975e-13 (engine=5.797523e-02, reviewer=5.797523e-02)
  region/component_1/momentum_rate: diff=-2.829e-13 (engine=2.072459e+00, reviewer=2.072459e+00)
  region/component_2/balance/residual: diff=-6.809e-09 (engine=1.753466e+00, reviewer=1.753466e+00)
  region/component_2/balance/operand_b: diff=4.086e-13 (engine=-1.321984e+00, reviewer=-1.321984e+00)
  region/component_3/momentum_rate: diff=3.845e-13 (engine=6.290291e-02, reviewer=6.290291e-02)
```

```
$ /usr/bin/grep -n -E '4839b124|811bea8c' /var/projects/_review_scratch/s15-build-r1-claude/REVIEW.md /var/projects/_review_scratch/s15-build-r1-grok/REVIEW.md /var/projects/_review_scratch/s15-build-r1-claude/item8/sympy_rerun.stdout /var/projects/_review_scratch/s15-build-r1-claude/item8/wolfram_rerun.stdout | sed 's|/var/projects/_review_scratch/||' | cut -c1-200
s15-build-r1-claude/REVIEW.md:99:  (`4839b124…` and `811bea8c…`).
s15-build-r1-grok/REVIEW.md:56:The run record names, for each engine, the guarded command, exit 0, `limitsVerified: true`, stderr bytes 0, the transcript path and its sha256. `outcome.json` in `_scrat
s15-build-r1-claude/item8/sympy_rerun.stdout:2:4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c  abl/py_r1_rerun/out/S15_two_drain_region.out
s15-build-r1-claude/item8/sympy_rerun.stdout:3:4839b124811e88b6c0cd99dda9a9518776418df93e16ac8073063e142f34fc8c  /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_two_drain_region.out
s15-build-r1-claude/item8/wolfram_rerun.stdout:2:811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6  abl/wl_r1_rerun/out/S15_two_drain_region.out
s15-build-r1-claude/item8/wolfram_rerun.stdout:3:811bea8c58f18b52566f3de9661e501994aa5d8cd9fcde5e3610ae2aefc215a6  /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_two_drain_region
```

## M7. The Claude leg's own stopped evaluator run

```
$ cat _scratch/s11c/s15-build-review-r1-guard/claude_evalwl_full01/outcome.json; ls _scratch/s11c/s15-build-review-r1-guard/claude_evalwl_full01/
{
  "exitCode": 1,
  "wallSeconds": 854.2176596270001,
  "unit": "s11c-guard-293a6f9936b7",
  "stderrBytes": 0,
  "limitsVerified": true,
  "childOutcome": null
}
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
outcome.json
resource-samples.jsonl
stderr
```

```
$ systemctl --user is-active s11c-guard-293a6f9936b7.service; ps -eo comm | /usr/bin/grep -c -E '^(WolframKernel|MathKernel)$'
inactive
0
```

## M8. The legs' scripts and stdout files exist (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-build-r1-claude/derive.py /var/projects/_review_scratch/s15-build-r1-claude/derive.stdout /var/projects/_review_scratch/s15-build-r1-claude/h1_check.py /var/projects/_review_scratch/s15-build-r1-claude/item0/h1_check.stdout /var/projects/_review_scratch/s15-build-r1-claude/item0/transcript_diff_r0_r1.stdout /var/projects/_review_scratch/s15-build-r1-claude/item0/prose_sweep.stdout /var/projects/_review_scratch/s15-build-r1-claude/compare.py /var/projects/_review_scratch/s15-build-r1-claude/compare_sympy.stdout /var/projects/_review_scratch/s15-build-r1-claude/compare_wolfram.stdout /var/projects/_review_scratch/s15-build-r1-claude/item9/sympy_ablation_diffs.stdout /var/projects/_review_scratch/s15-build-r1-claude/item9/wolfram_ablation_diffs.stdout /var/projects/_review_scratch/s15-build-r1-grok/derive_s15.py /var/projects/_review_scratch/s15-build-r1-grok/derive_s15.out /var/projects/_review_scratch/s15-build-r1-grok/check_h1.py /var/projects/_review_scratch/s15-build-r1-grok/logs/check_h1.out /var/projects/_review_scratch/s15-build-r1-grok/logs/r0_r1_line_multiset.out /var/projects/_review_scratch/s15-build-r1-grok/logs/ablation_summary.out | sed 's|/var/projects/_review_scratch/||'
7e2d79d31f18231e5b28228ee0d25e2f253b1654e41960498114517188b76877  s15-build-r1-claude/derive.py
e6cae55b6d02230f835922a011166fa2e8a2a72086cedce15dab208695eeb740  s15-build-r1-claude/derive.stdout
0c74c70280ce99bef17d77f27d8f60808d555365225341cab473635bd63701ed  s15-build-r1-claude/h1_check.py
074abe4d7b95d2d8a07734972a0dfdb2ccd35fa576be1ffd979dc7ff3ee7b37f  s15-build-r1-claude/item0/h1_check.stdout
05955b74f6c66067e449eec06d46d149250b01316970d24a3e91552d8916b6f9  s15-build-r1-claude/item0/transcript_diff_r0_r1.stdout
bf00771ce50c5dde1d3899922e297d018db9e25487345b38a60ee6cf9b86daa3  s15-build-r1-claude/item0/prose_sweep.stdout
1deaa59e83fe381bdaf57676da485753b79f55debe192963a1a3ac3ee1712f49  s15-build-r1-claude/compare.py
c5bc7a1fa2102c9ea902259e16cc3ff8d3194533da568ccb837710d40f54ad0e  s15-build-r1-claude/compare_sympy.stdout
e990f5c60acc82852305b15299d76efea534a77b1d603174df379417a3da52b2  s15-build-r1-claude/compare_wolfram.stdout
b1203daa0476691fee3784b9b0f960ba09d7dd725b1e917b0f1994c11905c6e4  s15-build-r1-claude/item9/sympy_ablation_diffs.stdout
762f82d6061cd287fad46bd39419ea72dd48cce7b69b14326e53e40c37058e35  s15-build-r1-claude/item9/wolfram_ablation_diffs.stdout
516e5c2d52ee3ed4e6e0f6c846d6e6ba9aa006decbe34bf9ee11b2d5a28739ac  s15-build-r1-grok/derive_s15.py
abd8a3a1f7d5b551f184441dbd9f948544af31532f8a5d7cc2a8e56190ac0044  s15-build-r1-grok/derive_s15.out
56e66aac86ad6f04f49f494d2e6af15ba889db91658c744744a5768d8dd5914f  s15-build-r1-grok/check_h1.py
4508d33bd75bad00347a8fc0020ec4708bb3d7d941cb9a748ccd223c0845226c  s15-build-r1-grok/logs/check_h1.out
3db302954f6341752e795221d808505ddeb140baa2f2f53a8e32aa593e7478d2  s15-build-r1-grok/logs/r0_r1_line_multiset.out
175717a3afbbbbafaa714bfd2594da1f50ac1f9815b148c4b2e35d61309cf1d8  s15-build-r1-grok/logs/ablation_summary.out
```
