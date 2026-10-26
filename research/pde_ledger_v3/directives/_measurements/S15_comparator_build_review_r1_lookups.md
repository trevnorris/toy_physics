# Measurements: S15 comparator build review, round 1 lookups

Generator: `_scratch/s15/gen/s15_cmp_r1_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval, literal-match counts and hashes. The legs' scripts
were not re-run here; their filed stdout is retrieved. `C/` is `/var/projects/_review_scratch/s15-cmp-r1-claude/`; `G/` is `/var/projects/_review_scratch/s15-cmp-r1-grok/`.

## M0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/cmp_review/s15_cmp_build_baseline_r1.sha256
research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py: OK
research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md: OK
research/pde_ledger_v3/directives/_measurements/S15_comparator_builder_report.md: OK
scripts/s11c_guarded_run.py: OK
_scratch/s15/cmp_work/full_07.jsonl: OK
_scratch/s15/cmp_work/full_07.audit.json: OK
_scratch/s15/cmp_work/full_07.coverage-audit.json: OK
_scratch/s15/cmp_work/tests_15.json: OK
research/pde_ledger_v3/directives/S15_comparator_decision_list.md: OK
```

```
$ sha256sum _scratch/s15/cmp_review/s15_cmp_build_review_prompt_r1.md
546ba3e76e7ed1d10c7585757b60367136f7e4da13ede86630d87ba56dac5c77  _scratch/s15/cmp_review/s15_cmp_build_review_prompt_r1.md
```

```
$ cat /var/projects/_review_scratch/s15-cmp-r1-grok/DONE; ls -la /var/projects/_review_scratch/s15-cmp-r1-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r1-claude/REVIEW.md | cut -c24-
exit 0
 trevnorris 11892 Oct 26 00:14 /var/projects/_review_scratch/s15-cmp-r1-claude/REVIEW.md
 trevnorris 10436 Oct 25 18:13 /var/projects/_review_scratch/s15-cmp-r1-grok/REVIEW.md
```

```
$ diff <(git status --porcelain | grep -v -F -e 'S15_comparator_build_review_r1_lookups.md' -e 'S15_comparator_build_review_r1_disposition.md'; git rev-parse HEAD) _scratch/s15/cmp_review/status_before_cmp_r1.txt && echo 'status unchanged apart from the r1 disposition and lookups files'
status unchanged apart from the r1 disposition and lookups files
```

```
$ systemctl --user list-units 's11c-guard-*' --state=active --no-legend | wc -l
0
```

## M1. The verdicts and findings

```
$ grep -n -E '^\*\*Verdict' /var/projects/_review_scratch/s15-cmp-r1-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r1-claude/REVIEW.md | sed 's|/var/projects/_review_scratch/||'
s15-cmp-r1-grok/REVIEW.md:3:**Verdict: CLEAR**
s15-cmp-r1-claude/REVIEW.md:3:**Verdict: CLEAR**
```

```
$ grep -n -A2 -F '## Findings' /var/projects/_review_scratch/s15-cmp-r1-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r1-claude/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-200
s15-cmp-r1-grok/REVIEW.md:7:## Findings
s15-cmp-r1-grok/REVIEW.md-8-
s15-cmp-r1-grok/REVIEW.md-9-None.
--
s15-cmp-r1-claude/REVIEW.md:12:## Findings
s15-cmp-r1-claude/REVIEW.md-13-
s15-cmp-r1-claude/REVIEW.md-14-None. No result below changes what the comparator computes on the accepted streams or what may be claimed from
```

## M2. What changed since the preserved round-0 baseline 0ad8773e

```
$ git diff --stat 0ad8773e -- research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py scripts/s11c_guarded_run.py research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md research/pde_ledger_v3/directives/_measurements/S15_comparator_builder_report.md
 .../_measurements/S15_comparator_builder_report.md |  32 +-
 .../_measurements/S15_comparator_run_record.md     |  80 ++++-
 .../scripts/test_S15_cross_engine_comparator.py    | 377 +++++++++++++++++----
 3 files changed, 408 insertions(+), 81 deletions(-)
```

```
$ grep -H -F -e 'research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py' -e 'scripts/s11c_guarded_run.py' _scratch/s15/cmp_review/s15_cmp_build_baseline_r0.sha256 _scratch/s15/cmp_review/s15_cmp_build_baseline_r1.sha256 | sed 's|_scratch/s15/cmp_review/||'
s15_cmp_build_baseline_r0.sha256:786ff802b0bcf54b1e8103f596bc7af3216eafc800a524106e23ee293fdd11aa  research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py
s15_cmp_build_baseline_r0.sha256:86be9b55a56f7a3538cdf3852b74562bdc8c5bab19c9ad894cbc2c3b6bbe257b  scripts/s11c_guarded_run.py
s15_cmp_build_baseline_r1.sha256:786ff802b0bcf54b1e8103f596bc7af3216eafc800a524106e23ee293fdd11aa  research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py
s15_cmp_build_baseline_r1.sha256:86be9b55a56f7a3538cdf3852b74562bdc8c5bab19c9ad894cbc2c3b6bbe257b  scripts/s11c_guarded_run.py
```

## M3. C1: each control altered at its own step, the rest of the comparison on

```
$ grep -n -F 'disable_comparison' research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200
research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py:1180:    if options.get('disable_comparison'): return {'not_compared':'comparison disabled'}
research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py:96:# In particular, no fault uses disable_comparison or removes all comparisons.
```

```
$ grep -n -o -F 'No global `disable_comparison` option is used.' research/pde_ledger_v3/directives/_measurements/S15_comparator_run_record.md
37:No global `disable_comparison` option is used.
```

```
$ for f in /var/projects/_review_scratch/s15-cmp-r1-claude/tests_out/c1_*.stdout.txt; do echo "$(basename $f): $(grep -h -E 'AssertionError: |^ERRORS' $f | sed 's/^ *| *//' | tail -2 | tr '\n' ' ')"; done | cut -c1-220
c1_dimension.stdout.txt: AssertionError: mutation invisible to exact residual ERRORS ['binding_repoints'] 
c1_domain.stdout.txt: AssertionError: mutation has no targeted comparison assertion ERRORS ['labels_and_conditions'] 
c1_exact.stdout.txt: AssertionError: mutation invisible to exact residual ERRORS ['basic_mutations'] 
c1_formation.stdout.txt: AssertionError: targeted control survived: mutation invisible to exact_residual ERRORS ['basic_mutations'] 
c1_leafcredit.stdout.txt: AssertionError: omitted leaf credit invisible ERRORS ['packaging_and_accounting'] 
c1_ordinal.stdout.txt: AssertionError: mutation invisible to comparison ERRORS ['repeated_records'] 
c1_r0_global_disable.stdout.txt: AssertionError: residual fault disabled evaluation ERRORS ['basic_mutations'] 
c1_sample.stdout.txt: AssertionError: mutation invisible to samples ERRORS ['basic_mutations'] 
c1_text.stdout.txt: AssertionError: mutation invisible to comparison ERRORS ['labels_and_conditions'] 
```

```
$ sed -n '77,79p' /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2.stdout | cut -c1-200
======== C1 ablation restoring global disable_comparison
disabled comparison record {'not_compared': 'comparison disabled', 'exact_residual': None, 'samples': None}
ROUND0_DISABLE_FAIL mutation invisible to samples
```

## M4. C2: the sparse-moment oracle

```
$ grep -n -F -e 'c.angular_moment' -e 'def fourier_moment' -e 'def sparse_integral_oracle' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200
908:def sparse_integral_oracle():
915:    def fourier_moment(p,q,full):
946:        original=c.angular_moment
950:            clear();c.angular_moment=lambda p,q,full:original(p,q,full)/(2 if full else 1)
952:        finally:c.angular_moment=original;clear()
```

```
$ sed -n '909,911p;928,930p' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200
    # Independent dense Fourier expansion on synthetic polynomials only.
    # The oracle integrates each exponential mode at its endpoints. It neither
    # calls angular_moment nor uses that function's beta/gamma moment formula.
        for powers,coefficient in polynomial.terms():
            value=coefficient*fourier_moment(*powers[:2],False)*fourier_moment(*powers[2:4],True)
            if radial:value*=(s.Rational(5,2)**(powers[4]+1)-s.Rational(2,3)**(powers[4]+1))/(powers[4]+1)
```

```
$ for f in /var/projects/_review_scratch/s15-cmp-r1-claude/tests_out/c2_*.stdout.txt; do echo "$(basename $f): $(grep -h -E 'AssertionError: |^ERRORS' $f | sed 's/^ *| *//' | tail -2 | tr '\n' ' ')"; done | cut -c1-220
c2_full_halved.stdout.txt: AssertionError: sparse moments differ from independent Fourier expansion ERRORS ['sparse_integral_oracle'] 
c2_full_halved_r0oracle.stdout.txt: ERRORS [] 
c2_half_halved.stdout.txt: AssertionError: sparse moments differ from independent Fourier expansion ERRORS ['sparse_integral_oracle'] 
c2_r0oracle_only.stdout.txt: ERRORS [] 
```

```
$ sed -n '8p' /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2.stdout; sed -n '1,5p' /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_02.stdout | cut -c1-200
pairs 98 mismatches 0
form_ablate start
======== clarify halved moment against an oracle that calls angular_moment
halved_integral 2*pi/3
oracle_from_halved_angular_moment 2*pi/3
round0_oracle_matches_halved_production True
```

```
$ sed -n '80,92p' /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2.stdout | cut -c1-200
======== C2 ablation: oracle built from angular_moment after it is halved
unhalved integrate_exact 4*pi/3
fourier product 4*pi/3
unhalved_equals_fourier True
halved integrate_exact 2*pi/3
halved minus itself 0
shared_halved_oracle_equal True
fourier_reference 4*pi/3
halved_production_moment 2*pi/3
fourier_minus_halved 2*pi/3
round0_shared_moment_would_match False
repaired_fourier_detects_halving True
c1c2 probe done
```

```
$ grep -n -F -e 'reference = comparator.angular_moment' -e 'production = s.simplify(comparator.angular_moment' -e 'round0_shared_moment_would_match' /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2_probe.py | cut -c1-200
132:        reference = comparator.angular_moment(1, 0, False) * comparator.angular_moment(0, 2, True)
143:        production = s.simplify(comparator.angular_moment(3, 0, False) * comparator.angular_moment(0, 2, True))
148:        print(f"round0_shared_moment_would_match {s.simplify(production - reference) == 0}", flush=True)
```

## M5. Item 10: the operands that evaluate to zero

```
$ head -10 /var/projects/_review_scratch/s15-cmp-r1-claude/mine/item10_join.stdout.txt | cut -c1-200
label | samples | max |py value|/int|f| | max |wl value|/int|f| | max |py-wl|/int|f| | comparator py/wl sample operands (distinct)
CONTROL(drain1 sibling of 1373) | 45 (wl 45; integrand identically 0: 0) | 1.00e+00 | 1.00e+00 | 1.23e-15 | n/a (control)
CONTROL(drain1 sibling of 151) | 45 (wl 45; integrand identically 0: 0) | 8.07e-01 | 8.07e-01 | 1.55e-15 | n/a (control)
region_momentum_groups:1373 | 45 (wl 45; integrand identically 0: 0) | 3.25e-16 | 3.76e-17 | 2.87e-16 | [('0', '0')]
region_rate_groups:1308 | 45 (wl 45; integrand identically 0: 0) | 4.80e-16 | 2.61e-16 | 2.66e-16 | [('0', '0')]
surface_groups:1187 | 45 (wl 45; integrand identically 0: 12) | 5.24e-16 | 2.92e-16 | 4.65e-16 | [('0', '0')]
surface_groups:151 | 45 (wl 45; integrand identically 0: 0) | 3.28e-16 | 2.25e-16 | 1.58e-16 | [('0', '0')]
surface_groups:601 | 45 (wl 45; integrand identically 0: 10) | 2.50e-16 | 1.76e-16 | 2.50e-16 | [('0', '0')]
surface_groups:737 | 45 (wl 45; integrand identically 0: 0) | 7.36e-16 | 4.11e-16 | 4.05e-16 | [('0', '0')]
```

```
$ grep -n -E '^(phi_integral_of_sin|independent_integral|value_over_scale|py_minus_wl|both_zero)' /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_03.stdout | cut -c1-200
18:phi_integral_of_sin 0
19:independent_integral 0
20:independent_integral_is_zero True
27:value_over_scale 0
43:phi_integral_of_sin 0
44:independent_integral 0
45:independent_integral_is_zero True
52:value_over_scale 0
53:py_minus_wl 0
59:both_zero True
62:value_over_scale undefined because the integrand is the zero expression
```

```
$ python3 -c "import json,collections; R=[json.loads(l) for l in open('_scratch/s15/cmp_work/full_07.bindings/001373.jsonl') if l.strip()]; S=[r for r in R if r.get('kind')=='sample']; print('samples', len(S), dict(collections.Counter((str(r.get('py')), str(r.get('wl')), str(r.get('residual'))) for r in S)))"
samples 45 {('0', '0', '0'): 45}
```

## M6. Reruns and splitting against full_07

```
$ grep -h '^SUMMARY' /var/projects/_review_scratch/s15-cmp-r1-claude/runs/compare_B_base_vs_full07.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-claude/runs/compare_A_w1_vs_full07.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-claude/runs/compare_A_w2_reverse_vs_full07.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-claude/runs/compare_A_w3_forward_vs_full07.stdout.txt
SUMMARY bindings=1470 identical=1470 differing=0
SUMMARY bindings=18 identical=18 differing=0
SUMMARY bindings=18 identical=18 differing=0
SUMMARY bindings=18 identical=18 differing=0
```

```
$ grep -h -E '^records|^differing' /var/projects/_review_scratch/s15-cmp-r1-claude/runs/whole_file_split.stdout.txt
records 4438 4438
differing record positions 0 []
records 4438 4438
differing record positions 0 []
```

```
$ grep -h -E '^RERUN|worker split' /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_02.stdout /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_03.stdout | cut -c1-200
RERUN global_domain:9 matches_full07 True samples 45
RERUN prose:0 matches_full07 True samples 0
RERUN surface_piece:476 matches_full07 True samples 45
RERUN surface_pressure_groups:281 matches_full07 True samples 45
RERUN region_momentum_groups:1373 matches_full07 True samples 45
RERUN region_rate_groups:1328 matches_full07 True samples 45
worker split mathematical records equal True count 3
```

```
$ sed -n '3,4p' /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2.stdout
bindings 1616 unique_py_ends 1616 unique_wl_ends 1616
names 45 unique_py 45 unique_wl 45
```

## M7. Ablations that move production and no test (f01, f02, f03, f06; Claude leg)

```
$ for f in f01 f02 f03 f06; do grep -h -E '^f0[0-9] +[0-9]+ ' /var/projects/_review_scratch/s15-cmp-r1-claude/runs/effect_$f.stdout.txt; done | cut -c1-200
f01    1288  identical
f01     170  exact_residual_changed+nonzero_sample_residuals
f01      12  exact_residual_changed
f02    1377  identical
f02      81  exact_residual_changed+nonzero_sample_residuals
f02      12  now_not_formed+dimension_difference_changed
f03    1434  identical
f03      36  exact_residual_changed+nonzero_sample_residuals
f06     904  identical
f06     566  sample_count_45_to_36+sample_operands_moved_residuals_still_zero
```

```
$ for f in f01_surfaces_by_index f02_drains_exchanged f03_pi_transposed f06_stratum_dropped; do echo "$f: $(grep -h '^ERRORS' /var/projects/_review_scratch/s15-cmp-r1-claude/tests_out/$f.stdout.txt)"; done
f01_surfaces_by_index: ERRORS []
f02_drains_exchanged: ERRORS []
f03_pi_transposed: ERRORS []
f06_stratum_dropped: ERRORS []
```

## M8. Test-function distinctness and the Grok leg's density collision

```
$ sed -n '716,724p' research/pde_ledger_v3/scripts/test_S15_cross_engine_comparator.py | cut -c1-200

def test_functions():
    # Distinctness, dependence and domain are checked on every declared profile
    # in every form/stratum, not inferred from the profile generator's loops.
    for form in c.FORMS:
        for stratum in c.STRATA:
            for index in range(len(c.POINTS)):
                polys={r.py:c.profile(r.py,r.arity,form,stratum,index) for r in c.NAME_TABLE if r.arity}
                assert len(set(polys.values()))==len(polys),'live test functions collided'
```

```
$ sed -n '825,826p' research/pde_ledger_v3/scripts/S15_cross_engine_comparator.py | cut -c1-200
def profile(name,arity,form,stratum,point,options=None):
    return profile_cached(name,arity,form,stratum,point,tuple(sorted((options or {}).items())))
```

```
$ sed -n '73,92p' /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_fix.py | cut -c1-200
    p1 = c.profile_cached("rho_1", 2, c.FORMS[0], "generic", 0, ())
    p2 = c.profile_cached("rho_2", 2, c.FORMS[0], "generic", 0, ())
    print(f"rho_profiles_differ {p1 != p2}", flush=True)
    original_profile = c.profile_cached

    def shared_density(name, arity, form, stratum, point, opts=()):
        if name == "rho_2":
            name = "rho_1"
        return original_profile(name, arity, form, stratum, point, opts)

    c.profile_cached = shared_density
    c.profile_jet.cache_clear()
    fa.clear(c)
    try:
        collided = fa.comparison(fa.run_rows(c, streams, names, conditions, [piece], options), piece.key)
    finally:
        c.profile_cached = original_profile
        c.profile_jet.cache_clear()
        fa.clear(c)
    fa.show_move("rho_2 evaluated as rho_1", base, collided)
```

```
$ grep -n -E 'rho_profiles_differ|rho_2 evaluated as rho_1' /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_04.stdout | cut -c1-200
14:rho_profiles_differ True
15:ABLATION rho_2 evaluated as rho_1 moved True
```

## M9. The legs' guarded jobs

```
$ grep -c -E ' 0 True 0 None$' /var/projects/_review_scratch/s15-cmp-r1-claude/runs/guard_outcomes.stdout.txt; grep -v -E ' 0 True 0 None$' /var/projects/_review_scratch/s15-cmp-r1-claude/runs/guard_outcomes.stdout.txt | grep -v '^$' | cut -c1-200
27
followup_01 (no outcome.json — running or not admitted)
item10_py_01 1 True 1 None
```

```
$ for f in /var/projects/_review_scratch/s15-cmp-r1-grok/*.outcome.json; do echo "$(basename $f) $(tr -d ' \n' < $f | grep -o -E '"exitCode":[0-9]+' | head -1) $(tr -d ' \n' < $f | grep -o -E '"limitsVerified":[a-z]+')"; done
c1c2.outcome.json "exitCode":0 "limitsVerified":true
form_ablate_02.outcome.json "exitCode":1 "limitsVerified":true
form_ablate_03.outcome.json "exitCode":0 "limitsVerified":true
form_ablate_04.outcome.json "exitCode":0 "limitsVerified":true
own_compare_05.outcome.json "exitCode":1 "limitsVerified":true
own_compare_06.outcome.json "exitCode":0 "limitsVerified":true
```

## M10. The legs' reports, scripts and stdout (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-cmp-r1-claude/REVIEW.md /var/projects/_review_scratch/s15-cmp-r1-claude/mine/make_copies.py /var/projects/_review_scratch/s15-cmp-r1-claude/mine/run_tests.py /var/projects/_review_scratch/s15-cmp-r1-claude/mine/ablation_effect.py /var/projects/_review_scratch/s15-cmp-r1-claude/mine/item10_py.py /var/projects/_review_scratch/s15-cmp-r1-claude/mine/item10_wl.wl /var/projects/_review_scratch/s15-cmp-r1-claude/mine/item10_join.py /var/projects/_review_scratch/s15-cmp-r1-claude/mine/item10_join.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-claude/mine/compare.py /var/projects/_review_scratch/s15-cmp-r1-claude/mine/compare_01.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-claude/runs/compare_B_base_vs_full07.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-claude/runs/whole_file_split.stdout.txt /var/projects/_review_scratch/s15-cmp-r1-grok/REVIEW.md /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2_probe.py /var/projects/_review_scratch/s15-cmp-r1-grok/c1c2.stdout /var/projects/_review_scratch/s15-cmp-r1-grok/own_compare.py /var/projects/_review_scratch/s15-cmp-r1-grok/own_compare.stdout /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate.py /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_02.stdout /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_rest.py /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_03.stdout /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_fix.py /var/projects/_review_scratch/s15-cmp-r1-grok/form_ablate_04.stdout | sed 's|/var/projects/_review_scratch/||'
28ac571ca0c7715591e0f3774c48da516333ec99b7370e263af853bfe584f79d  s15-cmp-r1-claude/REVIEW.md
e774c8bbe52ecf061eb4ce810b0482e45f8aefd35597ec6667026d51ff16fa6a  s15-cmp-r1-claude/mine/make_copies.py
d7a7fff6f4a071f99baa57acb3869cda9d2ae689b2aa58501a98580e30222d87  s15-cmp-r1-claude/mine/run_tests.py
7c478487d300eed74de305fb94c5c8d3e962a244186816daa64cc66e5001a721  s15-cmp-r1-claude/mine/ablation_effect.py
3e763d619d1b1dd6ac2164c7d4542648e991afa5442e11e1f5c623436540d2f3  s15-cmp-r1-claude/mine/item10_py.py
cc4d14247ae1c9eaa5405dd01024d824a8679f0bd4508118c28495d426e4645f  s15-cmp-r1-claude/mine/item10_wl.wl
573a56fd1832e38019e4870e24398f4eb0baf35e999e7606de5a74896357e1c8  s15-cmp-r1-claude/mine/item10_join.py
9dda49738a2c924183812f0deb7bba6b44b5df6d2fc72cd633ec9381bae6f933  s15-cmp-r1-claude/mine/item10_join.stdout.txt
239c0b079883f7ec0836345ff97891e064600cec8e140740baf6da6b5a17031f  s15-cmp-r1-claude/mine/compare.py
def4281a366caa80eec847da032d90fb8924932babe65d192b9788eadb7a33d2  s15-cmp-r1-claude/mine/compare_01.stdout.txt
000337563fe1fa8673a2daa0d339b955254a0f3d1f1eda84f315773d1bc308ac  s15-cmp-r1-claude/runs/compare_B_base_vs_full07.stdout.txt
d6b9b4fa5a7f08ac10a9f2e2e167151a27e2344609733827af3b9fa06c9a383c  s15-cmp-r1-claude/runs/whole_file_split.stdout.txt
7e9f1a9c0797fc96ada184e4d9fdbf465677b9673698695c4bdf08df1ff64001  s15-cmp-r1-grok/REVIEW.md
3581f78035358f3e0d533a676637869b6e4477cf673a73f55313beff558c4aa1  s15-cmp-r1-grok/c1c2_probe.py
d855c059859e3ed1257e4736f650249f03ba8ad66b5bc103541e1ec0c7cb53f8  s15-cmp-r1-grok/c1c2.stdout
2a0b18db0cc965e060e681970d434ce8ca02a47fcc4ddfbbbed7189e3570dedd  s15-cmp-r1-grok/own_compare.py
ec7627e0824344aa11c44b49dc0b6b16028584dfd279af3d8c188d8d67d13769  s15-cmp-r1-grok/own_compare.stdout
2bff99da838cb62855e9614acb545a6f55fd47b519a1ff794e3cc1c31c7af21e  s15-cmp-r1-grok/form_ablate.py
5afa1ed3932d498dac6dce37464ea27d68ef82546b4f9918386bf51975ec580b  s15-cmp-r1-grok/form_ablate_02.stdout
2f472e534541c35702cb36b96802f790fddcf9ed3f352a65441a510861e14a38  s15-cmp-r1-grok/form_ablate_rest.py
53667c658714f3f7f6781cfc7f76cae0990d3df38c79d015baa55981c977cbdc  s15-cmp-r1-grok/form_ablate_03.stdout
c1a706bbebf911339fa7fc3cd01b5a2a90cdfb41f0bec9038bc15286b3aa966c  s15-cmp-r1-grok/form_ablate_fix.py
6e97a7e3d8a53292af62ff66ad4f17d8e73ae8f505de2d94cffa9b63c080c9e1  s15-cmp-r1-grok/form_ablate_04.stdout
```
