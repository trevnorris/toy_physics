# Measurements: S15 comparator build review, round 0 lookups

Generator: `_scratch/s15/gen/s15_cmp_r0_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts and hashes. The legs' scripts
were not re-run here; their filed stdout is retrieved.

## M0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/cmp_review/s15_cmp_build_baseline_r0.sha256
research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py: OK
research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md: OK
research/pde_ledger_v3/directives/_measurements/S15_comparator_builder_report.md: OK
scripts/s11c_guarded_run.py: OK
_scratch/s15/cmp_work/full_07.jsonl: OK
_scratch/s15/cmp_work/full_07.audit.json: OK
_scratch/s15/cmp_work/full_07.coverage-audit.json: OK
_scratch/s15/cmp_work/tests_12.json: OK
research/pde_ledger_v3/directives/S15_comparator_decision_list.md: OK
```

```
$ sha256sum _scratch/s15/cmp_review/s15_cmp_build_review_prompt_r0.md _scratch/s15/cmp_review/s11c_guard_diff_r0.patch
28ef471f77a2b18896b095117b985cef7573daa6bc65ddd55b4935913464cee4  _scratch/s15/cmp_review/s15_cmp_build_review_prompt_r0.md
2f2667a09fa32773938c4824dafcde535483074ccc632901643d9bc381487f6e  _scratch/s15/cmp_review/s11c_guard_diff_r0.patch
```

```
$ cat /var/projects/_review_scratch/s15-cmp-r0-grok/DONE; ls -la /var/projects/_review_scratch/s15-cmp-r0-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md | cut -c24-
exit 0
 trevnorris 11608 Oct 25 16:18 /var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md
 trevnorris  6342 Oct 25 11:47 /var/projects/_review_scratch/s15-cmp-r0-grok/REVIEW.md
```

```
$ diff <(git status --porcelain | grep -v -F -e 'S15_comparator_build_review_r0_lookups.md' -e 'S15_comparator_build_review_r0_disposition.md'; git rev-parse HEAD) _scratch/s15/cmp_review/status_before_cmp_r0.txt && echo 'status unchanged apart from the r0 disposition and lookups files'
status unchanged apart from the r0 disposition and lookups files
```

## M1. The verdicts

```
$ grep -n -E '^\*\*Verdict' /var/projects/_review_scratch/s15-cmp-r0-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md | sed 's|/var/projects/_review_scratch/||'
s15-cmp-r0-grok/REVIEW.md:3:**Verdict: CLEAR**
s15-cmp-r0-claude/REVIEW.md:3:**Verdict: NEEDS REVISION.** Two findings, both about what may be claimed for the **controls** (list item 7).
```

## M2. Finding C1: the comparison-disable switch

```
$ grep -n -F 'disable_comparison' research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200
research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py:1180:    if options.get('disable_comparison'): return {'not_compared':'comparison disabled'}
research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py:773:            OVERRIDE['disable_comparison']=True
```

```
$ grep -c -F 'disable' research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py
1
```

```
$ sed -n '96,99p' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200

def detects(base,changed,evaluation=False,engine='py',exact=False):
    assert comps(base) and all('not_compared' not in q and 'not_formed' not in q for q in comps(base)), 'baseline comparison inactive'
    assert computed(base)!=computed(changed),'mutation invisible to comparison'
```

```
$ sed -n '768,780p' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200
        # Pure test-profile construction and duplicate-table checks do not target
        # the comparison path, and have their direct defect controls above.
        for fn in TESTS:
            if fn in (test_functions,exclusive_guard_controls):continue
            print('WORKING_ON: comparison ablation '+fn.__name__,flush=True)
            OVERRIDE['disable_comparison']=True
            try:fn()
            except (AssertionError,KeyError):COUNTS['disable_'+fn.__name__]='control rejected reduction'
            else:raise AssertionError('control survived disabled comparison: '+fn.__name__)
            finally:OVERRIDE.clear()
    except Exception:
        args.output.write_text(json.dumps({'completed':COUNTS,'exception':traceback.format_exc(),'runtime_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_m
        raise
```

```
$ grep -c -F 'baseline comparison inactive' _scratch/s11c/s15-cmp-review-r0-claude/tests_base_disable_01/stdout; grep -F -e ': ASSERTION' -e ': ERROR' _scratch/s11c/s15-cmp-review-r0-claude/tests_base_disable_01/stdout | sed 's/ (.*//' | sort | uniq -c | sort -rn | head -12
15
      1 supplied_equation_and_alias: ASSERTION: baseline comparison inactive
      1 split_integral_evaluation: ASSERTION: linear integral evaluation changed value
      1 sparse_integral_oracle: ASSERTION: sparse moments differ from dense expansion
      1 repeated_records: ASSERTION: baseline comparison inactive
      1 quadrature: ASSERTION: baseline comparison inactive
      1 pieces_groups_and_components: ASSERTION: baseline comparison inactive
      1 partial_time_semantics: ASSERTION: partial time became chart time
      1 parallel_schedule: ASSERTION: parallel comparison inactive
      1 packaging_and_accounting: ASSERTION: baseline comparison inactive
      1 nested_radial_derivative: ASSERTION: nested radial derivative lost its exact polynomial value
      1 nested_integral_parameters: ASSERTION:            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      1 name_repoints: ASSERTION: baseline comparison inactive
```

```
$ grep -n -i -F 'comparison-disable' research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md | cut -c1-200
15:Tests: guard exit `0`, `limitsVerified: true`; child exit `0`, `guardReason: null`. Twenty-seven synthetic control groups completed and twenty-five comparison-disable groups rejected their reductio
```

## M3. Finding C2: the sparse-moment oracle

```
$ sed -n '695,708p' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200

def sparse_integral_oracle():
    # Independent dense expansion on synthetic polynomials only. No measured
    # transcript or source is loaded. Compare through both native readers.
    th,ph,r=s.symbols('theta phi s')
    sn,cs,sp,cp=s.sin(th),s.cos(th),s.sin(ph),s.cos(ph)
    def dense(integrand,radial=False):
        generators=(sn,cs,sp,cp,r) if radial else (sn,cs,sp,cp)
        polynomial=s.Poly(s.expand(integrand),*generators)
        terms=[]
        for powers,coefficient in polynomial.terms():
            value=coefficient*c.angular_moment(*powers[:2],False)*c.angular_moment(*powers[2:4],True)
            if radial:value*=(s.Rational(5,2)**(powers[4]+1)-s.Rational(2,3)**(powers[4]+1))/(powers[4]+1)
            terms.append(value)
```

```
$ grep -n -F 'def angular_moment' research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py; grep -c -F 'c.angular_moment' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py
853:def angular_moment(p,q,full):
1
```

```
$ sed -n '/^== A10_sparse_alter_angular/,/^== A11/p' /var/projects/_review_scratch/s15-cmp-r0-claude/ablate/ablation_diffs.txt | head -30
== A10_sparse_alter_angular/
856c856
<         return 2*s.gamma(s.Rational(p+1,2))*s.gamma(s.Rational(q+1,2))/s.gamma(s.Rational(p+q+2,2))
---
>         return s.gamma(s.Rational(p+1,2))*s.gamma(s.Rational(q+1,2))/s.gamma(s.Rational(p+q+2,2))
syntax-ok
== A11_linearity_drop_term/
```

```
$ grep -c -F ': PASSED' _scratch/s11c/s15-cmp-review-r0-claude/abl_A10_01/stdout; grep -c -E ': (ASSERTION|ERROR)' _scratch/s11c/s15-cmp-review-r0-claude/abl_A10_01/stdout; grep -F 'sparse_integral_oracle' _scratch/s11c/s15-cmp-review-r0-claude/abl_A10_01/stdout | cut -c1-120
26
0
sparse_integral_oracle: PASSED (13.8s)
```

```
$ head -12 /var/projects/_review_scratch/s15-cmp-r0-claude/ablate/out/diff_A10.txt | cut -c1-200
A out/A10_sparse_alter_angular.jsonl bindings 17
B /var/projects/toy_physics/_scratch/s15/cmp_work/full_07.jsonl
binding-table rows of A also present identically in B: 17 of 17
name_table identical: True | strata True | forms True
test_functions records identical: True (45 vs 45)
prose:0 operands=same(1/1) exact_comparison=same(0/0) comparison=same(1/1) binding_accounting=same(1/1) sample=same(0/0)
global_domain:9 operands=same(1/1) exact_comparison=same(1/1) comparison=same(1/1) binding_accounting=same(1/1) sample=same(45/45)
dimension_basis:11 operands=same(1/1) exact_comparison=same(0/0) comparison=same(1/1) binding_accounting=same(1/1) sample=same(0/0)
supplied_Pi:51 operands=same(1/1) exact_comparison=same(1/1) comparison=same(1/1) binding_accounting=same(1/1) sample=same(45/45)
supplied_Pi:53 operands=same(1/1) exact_comparison=same(1/1) comparison=same(1/1) binding_accounting=same(1/1) sample=same(45/45)
homogeneity_dimension_operands:145 operands=same(1/1) exact_comparison=same(0/0) comparison=same(1/1) binding_accounting=same(1/1) sample=same(0/0)
surface_pressure_first_order_choice:252 operands=same(1/1) exact_comparison=same(1/1) comparison=DIFF(1/1) binding_accounting=same(1/1) sample=DIFF(45/45)
```

## M4. The Claude leg's independent value checks (filed stdout)

```
$ tail -5 /var/projects/_review_scratch/s15-cmp-r0-claude/own/compare_own.stdout | cut -c1-200
   A/base[1]: py -2.517779287379287  wl -2.517779287379287  residual 0.000e+00  rel 0.000e+00
   A/base[2]: py -3.079134545454545  wl -3.079134545454545  residual 0.000e+00  rel 0.000e+00
   B/base[0]: py 0.1152221210170731  wl 0.1152221210170731  residual -1.388e-17  rel 1.204e-16
   B/base[1]: py 1.360192992350824  wl 1.360192992350824  residual 2.220e-16  rel 1.632e-16
   B/base[2]: py 0.8926054111670357  wl 0.8926054111670357  residual 0.000e+00  rel 0.000e+00
```

```
$ grep -E 'rel|max' /var/projects/_review_scratch/s15-cmp-r0-claude/own/evalcheck_02.stdout | tail -6 | cut -c1-200
region_rate_groups:1320 (region_c1_rate_d1r): comparator py -2.679837264632786e+17 wl -2.679837264632786e+17 | mine -2.679837264632802e+17 fine -2.679837264632782e+17 | mine-comparator -1.600e+03 rel 
region_momentum_groups:1373 (region_c2_momentum_d2): comparator py 0 wl 0 | mine 0.0001047253608703613 fine 0.0001479387283325195 | mine-comparator 1.047e-04 rel 1.05e+296
region_balance_parts:1340 (region_c1_balance_parts_outer_pressure): comparator py 9.759059902143628e+49 wl 9.759059902143628e+49 | mine 9.759059902143682e+49 fine 9.759059902143611e+49 | mine-comparat
region_balance_parts:1338 (region_c1_balance_parts_inner_transport_radius): comparator py 9.207341868692515e+23 wl 9.207341868692515e+23 | mine 9.207341868692562e+23 fine 9.207341868692497e+23 | mine-
supplied_Pi:51 (point_Pi_12): comparator py 6.653454613973127e+34 wl 6.653454613973127e+34 | mine 6.653454613973127e+34 | mine-comparator 0.000e+00 rel 0.00e+00
supplied_Pi:53 (point_Pi_21): comparator py 6.653454613973127e+34 wl 6.653454613973127e+34 | mine 6.653454613973127e+34 | mine-comparator 0.000e+00 rel 0.00e+00
```

```
$ grep -F 'nonzero exact residual AND all samples zero' /var/projects/_review_scratch/s15-cmp-r0-claude/scan/zero_scan.stdout | cut -c1-200
nonzero exact residual AND all samples zero: []
```

```
$ head -3 /var/projects/_review_scratch/s15-cmp-r0-claude/ablate/out/diff_corrupt_vs_full07.txt | cut -c1-200
A out/corrupt_wl68.jsonl bindings 5
B /var/projects/toy_physics/_scratch/s15/cmp_work/full_07.jsonl
binding-table rows of A also present identically in B: 5 of 5
```

## M5. The Grok leg's reproduction and splitting lines (filed stdout)

```
$ grep -c -E '^REPRO .*exact_ok=True accounting_ok=True sample_mismatches=0' /var/projects/_review_scratch/s15-cmp-r0-grok/02_ablate_03.stdout; grep -c -E '^REPRO ' /var/projects/_review_scratch/s15-cmp-r0-grok/02_ablate_03.stdout
19
19
```

```
$ grep -h -E 'split_1_vs_2|split_forward_vs_reverse|^reader ' /var/projects/_review_scratch/s15-cmp-r0-grok/*.stdout | cut -c1-200
reader py_lines=2246 py_records=2231 py_bad=0 wl_lines=1370 wl_records=1370 wl_bad=0
reader py_lines=2246 py_records=2231 py_bad=0 wl_lines=1370 wl_records=1370 wl_bad=0
split_1_vs_2 True
split_forward_vs_reverse True
reader py_lines=2246 py_records=2231 py_bad=0 wl_lines=1370 wl_records=1370 wl_bad=0
```

```
$ grep -l -F 'disable' /var/projects/_review_scratch/s15-cmp-r0-grok/*.stdout /var/projects/_review_scratch/s15-cmp-r0-grok/*.txt 2>/dev/null | sed 's|/var/projects/_review_scratch/||'
```

## M6. The guard change (item 9)

```
$ grep -n -i -F 'guard change' /var/projects/_review_scratch/s15-cmp-r0-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-200
s15-cmp-r0-grok/REVIEW.md:5:**Findings.** None. On the accepted transcripts, and on the synthetic controls rerun here, the comparator computes and prints the comparison the decision list asks for. The
s15-cmp-r0-grok/REVIEW.md:45:**Guard change.** `s11c_guard_diff_r0.patch` adds opt-in `--all-cpus-exclusive`. Without the flag, pooled admission, `MemoryMax`, `MemorySwapMax=0`, `TasksMax`, and single
s15-cmp-r0-claude/REVIEW.md:108:  `--all-cpus-exclusive --memory-gib 55`, and the report does not mention the guard change.
```

```
$ tail -8 /var/projects/_review_scratch/s15-cmp-r0-claude/guardcheck/compare_guard.stdout | cut -c1-200
memory.max='27917287424': HEAD False working-tree False
memory.swap.max='0': HEAD True working-tree True
memory.swap.max='1': HEAD False working-tree False
pids.max='32': HEAD True working-tree True
pids.max='64': HEAD False working-tree False
nativeAddressSpaceBytes=[26843545600, 26843545600]: HEAD True working-tree True
nativeAddressSpaceBytes=[27917287424, 27917287424]: HEAD False working-tree False
cases 20 differing 2 | new accepts where HEAD refuses: 2
```

```
$ git diff --stat scripts/s11c_guarded_run.py
 scripts/s11c_guarded_run.py | 60 +++++++++++++++++++++++++++++++++------------
 1 file changed, 45 insertions(+), 15 deletions(-)
```

```
$ grep -n -F -e '--memory-gib' -e '--all-cpus-exclusive' research/pde_ledger_v3/directives/S15_comparator_decision_list.md research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md | cut -c1-200
research/pde_ledger_v3/directives/S15_comparator_decision_list.md:249:   - Run every CAS job through `scripts/s11c_guarded_run.py` in pooled mode (`--pool s15-comparator --memory-gib 50`),
research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md:12:python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-comparator --all-cpus-exclusive --memory-gib 55
research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md:18:python3 /var/projects/toy_physics/scripts/s11c_guarded_run.py --pool s15-comparator --all-cpus-exclusive --memory-gib 55
```

```
$ grep -c -F 's11c_guarded_run' research/pde_ledger_v3/directives/_measurements/S15_comparator_builder_report.md
0
```

## M7. The stopped A08 unit

```
$ cat _scratch/s11c/s15-cmp-review-r0-claude/abl_A08_01/outcome.json | tr -d '\n' | cut -c1-200; echo; ls -la --time-style=+%H:%M _scratch/s11c/s15-cmp-review-r0-claude/abl_A08_01/stderr | cut -c24-
{  "exitCode": 1,  "wallSeconds": 12096.200187718001,  "unit": "s11c-guard-1b2527c930bf",  "stderrBytes": 322782,  "limitsVerified": true,  "childOutcome": null}

 trevnorris 322782 13:38 _scratch/s11c/s15-cmp-review-r0-claude/abl_A08_01/stderr
```

## M8. Binding region_momentum_groups:1373: the printed operands and samples, and the Claude leg's quadrature lines

```
$ python3 -c "import json,collections; R=[json.loads(l) for l in open('_scratch/s15/cmp_work/full_07.bindings/001373.jsonl') if l.strip()]; print(dict(collections.Counter(r.get('kind') for r in R))); E=[r for r in R if r.get('kind')=='exact_comparison'][0]; print('exact_residual', E['exact_residual']); print('exact_method', E['exact_method']); S=[r for r in R if r.get('kind')=='sample']; print('samples', len(S), dict(collections.Counter((str(r.get('py')), str(r.get('wl')), str(r.get('residual'))) for r in S)))"
{'operands': 1, 'exact_comparison': 1, 'sample': 45, 'comparison': 1, 'binding_accounting': 1}
exact_residual 0
exact_method symbolic subtraction; integral linearity; predicates use symmetric difference; no supplied equation
samples 45 {('0', '0', '0'): 45}
```

```
$ grep -n -F '1373' /var/projects/_review_scratch/s15-cmp-r0-claude/own/evalcheck_02.stdout | cut -c1-250
91:region_momentum_groups:1373 (region_c2_momentum_d2): comparator py 0 wl 0 | mine 1.923181116580963e-07 fine 3.252644091844559e-07 | mine-comparator 1.923e-07 rel 1.92e+293
186:region_momentum_groups:1373 (region_c2_momentum_d2): comparator py 0 wl 0 | mine 0.0001047253608703613 fine 0.0001479387283325195 | mine-comparator 1.047e-04 rel 1.05e+296
```

```
$ grep -n -F -e '22 integrals' -e '5.5e-15' /var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md | cut -c1-200
120:   - 22 integrals × 2 configurations agree with my quadrature at ≤5.5e-15 relative, e.g.
```

## M9. The legs' scripts and stdout (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-cmp-r0-claude/REVIEW.md /var/projects/_review_scratch/s15-cmp-r0-claude/ablate/ablation_diffs.txt /var/projects/_review_scratch/s15-cmp-r0-claude/ablate/out/diff_A10.txt /var/projects/_review_scratch/s15-cmp-r0-claude/own/compare_own.py /var/projects/_review_scratch/s15-cmp-r0-claude/own/compare_own.stdout /var/projects/_review_scratch/s15-cmp-r0-claude/own/check_evaluator_compare.py /var/projects/_review_scratch/s15-cmp-r0-claude/own/evalcheck_02.stdout /var/projects/_review_scratch/s15-cmp-r0-claude/scan/zero_scan.py /var/projects/_review_scratch/s15-cmp-r0-claude/scan/zero_scan.stdout /var/projects/_review_scratch/s15-cmp-r0-claude/guardcheck/compare_guard.py /var/projects/_review_scratch/s15-cmp-r0-claude/guardcheck/compare_guard.stdout /var/projects/_review_scratch/s15-cmp-r0-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r0-grok/01_independent_compare.py /var/projects/_review_scratch/s15-cmp-r0-grok/01_independent_compare.stdout /var/projects/_review_scratch/s15-cmp-r0-grok/02_ablate.py /var/projects/_review_scratch/s15-cmp-r0-grok/02_ablate_03.stdout /var/projects/_review_scratch/s15-cmp-r0-grok/03_nested_param.py /var/projects/_review_scratch/s15-cmp-r0-grok/03_nested_param.stdout | sed 's|/var/projects/_review_scratch/||'
988826cc82c9f4ec0b29eaf01a62209f6f893a6abf710b3dd4ca22c2798189b6  s15-cmp-r0-claude/REVIEW.md
248f97b29f7d5314cdb327d230b86614a6368fb715fc950cabec5f27fc11c0f8  s15-cmp-r0-claude/ablate/ablation_diffs.txt
203e3a8dbe19dd1e1ad00ddfe57e2777ab43d21665471d4d82426ac83f68ed2b  s15-cmp-r0-claude/ablate/out/diff_A10.txt
b1c06658160e84bc9abc05d55fd2970204cf330eb87a37e4db9e29430a1b973f  s15-cmp-r0-claude/own/compare_own.py
3f7d39e7175c8ab73a0903179b776fcc3e38e7e04b3c9f219d86976ec1da73c8  s15-cmp-r0-claude/own/compare_own.stdout
6360e8ffcb4c85eca351047b2ce47a3a330206dda339a491bba79382d273e3c3  s15-cmp-r0-claude/own/check_evaluator_compare.py
be53d285120e837c596211d4dbafa13bcdc7cafa6ba316acaf2d98d50601dc5c  s15-cmp-r0-claude/own/evalcheck_02.stdout
7c201222819fc9306c41dd0f4c60f6a8aaa60fd0d974a3abdadac33aac913933  s15-cmp-r0-claude/scan/zero_scan.py
d1c38666b7cfa8205e2a8ccbfd2e08852cb01b1123b9c4131e7be50b151c2c8a  s15-cmp-r0-claude/scan/zero_scan.stdout
89e757ddc3b142c58882071ca58de3fd6b78ff206869f1d0179a8463d49d0c7b  s15-cmp-r0-claude/guardcheck/compare_guard.py
bc3d854dada701dc9c12869668b5661683c1d29a3d42bb3f498f9b4c23e45f63  s15-cmp-r0-claude/guardcheck/compare_guard.stdout
a411f5294d5d84d8876d2acde25fb6888c9448c4b1c7f6a725e612878f96c6cc  s15-cmp-r0-grok/REVIEW.md
da1cb5506a907ee3a5cb9296e0642856bacc0a8070a3f7bc3d84e3a4cdb03f55  s15-cmp-r0-grok/01_independent_compare.py
50eaf4e7c30737059b5832921efe4e06e00e6ecc892042d17aeba1fe32244cf5  s15-cmp-r0-grok/01_independent_compare.stdout
54c9cde32698f0974886bd49f11fa1c54df2e3bfc389932cfcae7a891af7e197  s15-cmp-r0-grok/02_ablate.py
d00a03fe41f2117755f822ef66785ce93aade9f5b8b7695f7632302d79a6ae4a  s15-cmp-r0-grok/02_ablate_03.stdout
069d0f0c8ea8cee56a8d5506672882ef8e9959bdd18c6ab0c2efc55c898b1f0a  s15-cmp-r0-grok/03_nested_param.py
3f9da8116680338b0d1f6cf70937699863bd0f042ba2f9833d74a7f2f1758f5a  s15-cmp-r0-grok/03_nested_param.stdout
```
