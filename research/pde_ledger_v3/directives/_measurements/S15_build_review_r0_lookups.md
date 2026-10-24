# Measurements: S15 build, review round 0 lookups

Generator: `_scratch/s15/gen/s15_build_r0_lookups.sh` (GNU grep, sed, git, sha256sum, `LC_ALL=C`). Each block is a
command and its literal output: existence, verbatim retrieval by line, literal-match counts, hashes. The legs'
scripts were not re-run here; their filed stdout is retrieved.

## M0. Baseline, prompt, legs finished, tree unchanged

```
$ sha256sum -c _scratch/s15/build_review/s15_build_baseline_r0.sha256
research/pde_ledger_v3/scripts/S15_two_drain_region.py: OK
research/pde_ledger_v3/mathematica/S15_two_drain_region.wl: OK
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out: OK
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out: OK
research/pde_ledger_v3/directives/_measurements/S15_build_run_record.md: OK
```

```
$ sha256sum _scratch/s15/build_review/s15_build_review_prompt_r0.md
5ebb9c2bfd8364099f66a8a662ca62eb6bc35c18f609662a5c1cc62191da5984  _scratch/s15/build_review/s15_build_review_prompt_r0.md
```

```
$ cat /var/projects/_review_scratch/s15-build-r0-grok/DONE; ls -la /var/projects/_review_scratch/s15-build-r0-claude/REVIEW.md | cut -c1-80
exit 0
-rw-rw-r-- 1 trevnorris trevnorris 10034 Oct 24 17:05 /var/projects/_review_scra
```

```
$ git status --porcelain | /usr/bin/grep -v -F 'S15_build_review_r0_lookups.md' | diff - <(sed -n '1,6p' _scratch/s15/build_review/status_before_r0.txt) && echo 'status unchanged apart from this file'
2d1
< ?? research/pde_ledger_v3/directives/_measurements/S15_build_review_r0_disposition.md
```

```
$ sed -n '7p' _scratch/s15/build_review/status_before_r0.txt; git rev-parse HEAD
a9acfb1d8ff85a8b33893be0c2af4624efb9abe8
a9acfb1d8ff85a8b33893be0c2af4624efb9abe8
```

## M1. The verdicts

```
$ /usr/bin/grep -n -E '^\*\*Verdict' /var/projects/_review_scratch/s15-build-r0-claude/REVIEW.md /var/projects/_review_scratch/s15-build-r0-grok/REVIEW.md | cut -c1-200
/var/projects/_review_scratch/s15-build-r0-claude/REVIEW.md:3:**Verdict: NEEDS REVISION.** There is one finding. It changes what may be claimed from the transcripts, not any printed object. On a gener
/var/projects/_review_scratch/s15-build-r0-grok/REVIEW.md:3:**Verdict: CLEAR**
```

## M2. The final guarded runs that wrote the transcripts

```
$ for d in sympy_20261024_03 wolfram_20261024_02; do python3 -c "import json,sys;o=json.load(open('_scratch/s11c/s15/'+sys.argv[1]+'/outcome.json'));c=o.get('childOutcome',{});print(sys.argv[1],'exitCode',o.get('exitCode'),'limitsVerified',o.get('limitsVerified'),'child',c.get('exitCode'),'guardReason',c.get('guardReason'),'stderrBytes',o.get('stderrBytes'))" $d; done
sympy_20261024_03 exitCode 0 limitsVerified True child 0 guardReason None stderrBytes 0
wolfram_20261024_02 exitCode 0 limitsVerified True child 0 guardReason None stderrBytes 0
```

## M3. H1: the region balance in each engine, and the sentence each transcript prints

```
$ sed -n '284,299p;310p' research/pde_ledger_v3/scripts/S15_two_drain_region.py | cut -c1-260
        rate=S.diff(momentum,t)
        forcing=volume_integral(force[i])
        for key,obj in [('momentum',momentum),('momentum_rate',rate),('integrated_force',forcing)]:
            solo(tag+'/'+key,obj); dimension_audit(tag+'/'+key,obj)
        buckets=groups(raw)
        for support,part in buckets.items():
            group_momentum=volume_integral(part)
            solo(tag+'/momentum_labels_'+names[support],group_momentum)
            solo(tag+'/rate_labels_'+names[support],S.diff(group_momentum,t))
        pair(tag+'/momentum_group_partition_integrand',strip(raw),sum(buckets.values()))
        parts=[surface_parts[0][i],surface_parts[1][i],forcing]
        rhs=sum(obj for surface in surface_parts for obj in surface[i].values())+forcing
        emit(tag+'/balance_parts',parts)
        pair(tag+'/balance',rate,rhs)
        emit(tag+'/balance_equation',S.Eq(rate,rhs,evaluate=False))
        dimension_audit(tag+'/balance_rhs',rhs)
    emit('balance_construction','The rate differentiates the moving-domain momentum integral. Applying divergence and moving-domain transport to the supplied local balance produces the printed surface parts and undivided volume force. The balance residual is a
```

```
$ sed -n '160,176p' research/pde_ledger_v3/mathematica/S15_two_drain_region.wl | cut -c1-260
Do[
 tag="region/component_"<>ToString[i];
 mom=vint[untag[momentum[[i]]]]; lhs=rate[untag[momentum[[i]]]]; forcing=vint[force[[i]]];
 solo[tag<>"/momentum",mom];dimensionAudit[tag<>"/momentum",mom];
 solo[tag<>"/momentum_rate",lhs];dimensionAudit[tag<>"/momentum_rate",lhs];
 groups=partition[momentum[[i]]];
 KeyValueMap[Function[{key,value},solo[tag<>"/momentum_labels_"<>key,vint[value]];solo[tag<>"/rate_labels_"<>key,rate[value]]],groups];
 pair[tag<>"/momentum_group_partition_integrand",untag[momentum[[i]]],Total[Values[groups]]];
 solo[tag<>"/integrated_force",forcing];dimensionAudit[tag<>"/integrated_force",forcing];
 rhs=Total[Flatten[Table[Values[surfaceResults[[j,i]]],{j,2}]]]+forcing;
 emit[tag<>"/balance_parts",{surfaceResults[[1,i]],surfaceResults[[2,i]],forcing}];
 pair[tag<>"/balance",lhs,rhs];emit[tag<>"/balance_equation",Equal[lhs,rhs]];
 dimensionAudit[tag<>"/balance_rhs",rhs];
 emit[tag<>"/balance_dimension_operand_a",dim[lhs]];emit[tag<>"/balance_dimension_operand_b",dim[rhs]];
 emit[tag<>"/balance_dimension_residuals",Flatten[Table[da-db,{da,dim[lhs]},{db,dim[rhs]}],1]],
 {i,3}];
emit["balance_construction", "Leibniz differentiation of the moving spherical-chart momentum integral supplies the rate. Integrating the supplied local balance and applying divergence and moving-domain transport gives the separately printed surface parts plus 
```

```
$ /usr/bin/grep -n -F 'balance_construction' research/pde_ledger_v3/scripts/out/S15_two_drain_region.out research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out | cut -c1-700
research/pde_ledger_v3/scripts/out/S15_two_drain_region.out:2241:balance_construction: The rate differentiates the moving-domain momentum integral. Applying divergence and moving-domain transport to the supplied local balance produces the printed surface parts and undivided volume force. The balance residual is an equation operand, not imposed on arbitrary independent live fields.
research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out:1365:balance_construction: "Leibniz differentiation of the moving spherical-chart momentum integral supplies the rate. Integrating the supplied local balance and applying divergence and moving-domain transport gives the separately printed surface parts plus undivided force. The residual is printed as an equation operand, not imposed on arbitrary independent live fields."
```

```
$ /usr/bin/grep -c -E '^region/component_[123]/balance/residual' research/pde_ledger_v3/scripts/out/S15_two_drain_region.out; /usr/bin/grep -c -E '^region/component_[123]/balance/residual' research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out
3
3
```

```
$ /usr/bin/grep -n -F 'From the supplied balance' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-200; /usr/bin/grep -n -F 'never states' research/pde_ledger_v3/directives/S15_decision_list.md | cut -c1-200
115:2. **The region's momentum balance.** From the supplied balance: the rate of the region's momentum, the integral of
126:A name may name an object, never its value or a comparison outcome. The script prints objects and never states
```

```
$ /usr/bin/grep -n -E 'BALANCE_x' /var/projects/_review_scratch/s15-build-r0-claude/check/eval_engine_sympy.stdout /var/projects/_review_scratch/s15-build-r0-claude/check/compare_wolfram_to_derive.stdout | sed 's|/var/projects/_review_scratch/||' | cut -c1-260
s15-build-r0-claude/check/eval_engine_sympy.stdout:310:engine_py	BALANCE_x1	rate=1.56131353576873	surface_parts=-2.41655584222939	int_f=3.97786937799814	rate_minus_parts_minus_f=-1.288e-14
s15-build-r0-claude/check/eval_engine_sympy.stdout:312:engine_py	BALANCE_x2	rate=0.190201430853231	surface_parts=-4.53076142222279	int_f=4.72096285307602	rate_minus_parts_minus_f=0.000e+00
s15-build-r0-claude/check/eval_engine_sympy.stdout:314:engine_py	BALANCE_x3	rate=0.0362683734792547	surface_parts=2.105392153478	int_f=-2.06912377999876	rate_minus_parts_minus_f=1.377e-14
s15-build-r0-claude/check/compare_wolfram_to_derive.stdout:310:engine_wl	BALANCE_x1	rate=1.5613135357687278	surface_parts=-2.416555842229388	int_f=3.9778693779981253	rate_minus_parts_minus_f=-9.769962616701378e-15
s15-build-r0-claude/check/compare_wolfram_to_derive.stdout:312:engine_wl	BALANCE_x2	rate=0.19020143085322927	surface_parts=-4.530761422222762	int_f=4.720962853075996	rate_minus_parts_minus_f=-4.440892098500626e-15
s15-build-r0-claude/check/compare_wolfram_to_derive.stdout:314:engine_wl	BALANCE_x3	rate=0.03626837347925446	surface_parts=2.1053921534779825	int_f=-2.0691237799987467	rate_minus_parts_minus_f=1.865174681370263e-14
```

```
$ cat /var/projects/_review_scratch/s15-build-r0-claude/ablate/W7*/summary.txt | /usr/bin/grep -E 'BALANCE' | cut -c1-260
  BALANCE_x1	rate=1.5613135357687278	surface_parts=4.136236261748875	int_f=3.9778693779981253	rate_minus_parts_minus_f=-6.552792103978273
  BALANCE_x2	rate=0.19020143085322927	surface_parts=5.963638363137003	int_f=4.720962853075996	rate_minus_parts_minus_f=-10.494399785359768
  BALANCE_x3	rate=0.03626837347925446	surface_parts=-3.1155307831691323	int_f=-2.0691237799987467	rate_minus_parts_minus_f=5.220922936647133
```

```
$ for a in W1_inner_normal_away W2_drop_X2_motion_in_vS W7_corrupt_balance_stress_sign; do cmp -s <(/usr/bin/grep -F 'balance_construction' /var/projects/_review_scratch/s15-build-r0-claude/ablate/$a/out/S15_two_drain_region.out) <(/usr/bin/grep -F 'balance_construction' research/pde_ledger_v3/mathematica/out/S15_two_drain_region.out) && echo "$a (wl): balance_construction line identical to the working-tree transcript" || echo "$a (wl): differs"; done; for a in P1_inner_normal_away P2_drop_X2_motion_in_vS; do cmp -s <(/usr/bin/grep -F 'balance_construction' /var/projects/_review_scratch/s15-build-r0-claude/ablate/$a/out/S15_two_drain_region.out) <(/usr/bin/grep -F 'balance_construction' research/pde_ledger_v3/scripts/out/S15_two_drain_region.out) && echo "$a (py): balance_construction line identical to the working-tree transcript" || echo "$a (py): differs"; done
W1_inner_normal_away (wl): balance_construction line identical to the working-tree transcript
W2_drop_X2_motion_in_vS (wl): balance_construction line identical to the working-tree transcript
W7_corrupt_balance_stress_sign (wl): balance_construction line identical to the working-tree transcript
P1_inner_normal_away (py): balance_construction line identical to the working-tree transcript
P2_drop_X2_motion_in_vS (py): balance_construction line identical to the working-tree transcript
```

```
$ /usr/bin/grep -n -F 'equation operand' /var/projects/_review_scratch/s15-build-r0-grok/REVIEW.md | cut -c1-260
45:**Item 2.** The rate of `∫ G^i`, both surfaces’ parts, and `∫ f^i` are printed apart (SymPy 284–297; Wolfram 87–90 and 162–171). `f` is the supplied vector, not decomposed. The region residual is the unevaluated difference of those objects. It i
```

## M4. Note: the free c_comp symbol in each engine

```
$ /usr/bin/grep -n -E 'c_comp|ccomp' research/pde_ledger_v3/scripts/S15_two_drain_region.py research/pde_ledger_v3/mathematica/S15_two_drain_region.wl | cut -c1-200
research/pde_ledger_v3/scripts/S15_two_drain_region.py:76:    solo('c_comp_not_present_in_supplied_equations', S.Function('c_comp')(*args))
research/pde_ledger_v3/scripts/S15_two_drain_region.py:170:    register(S.Function('c_comp'),sub(L,T),(L,L,L,T))
research/pde_ledger_v3/mathematica/S15_two_drain_region.wl:45:emit["LIVE_PARAMETER_NOT_IN_SUPPLIED_BALANCE", ccomp[x,y,z,t]];
research/pde_ledger_v3/mathematica/S15_two_drain_region.wl:105:register[ccomp,lenDim-timeDim,Join[ConstantArray[lenDim,3],{timeDim}]];
```

```
$ sed -n '367,370p' research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md
c_comp(ρ_br)² ≡ dp_br(ρ_br)/dρ_br .
```

`p_br` is a general function of `ρ_br` only; `c_comp` remains live wherever `ρ_br` varies.
```

## M5. Reproduction and controls the legs report (their stdout)

```
$ /usr/bin/grep -n -E '85d99aed|8774462d' /var/projects/_review_scratch/s15-build-r0-claude/REVIEW.md /var/projects/_review_scratch/s15-build-r0-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||' | cut -c1-220
s15-build-r0-claude/REVIEW.md:56:- **Item 8.** The run record names the command, guard outcome, transcript and sha256 for each engine, and these match the guard directories. Guarded reruns of unmodified copies reproduce 
s15-build-r0-grok/REVIEW.md:57:**Item 8.** The run record’s commands, exit codes, wall times, and sha256 values match `outcome.json`, `child-outcome.json`, and `invocation.json` in `sympy_20261024_03` and `wolfram_2026
```

```
$ tail -1 /var/projects/_review_scratch/s15-build-r0-claude/check/eval_engine_sympy.stdout; tail -1 /var/projects/_review_scratch/s15-build-r0-claude/check/compare_wolfram_to_derive.stdout
engine_py	N_OBJECTS_EVALUATED	309
engine_wl	N_OBJECTS_EVALUATED	309
```

## M6. The legs' scripts and stdout files exist (hashes)

```
$ sha256sum /var/projects/_review_scratch/s15-build-r0-claude/derive/derive_s15_v2.py /var/projects/_review_scratch/s15-build-r0-claude/derive/derive_s15_v2.stdout /var/projects/_review_scratch/s15-build-r0-claude/check/eval_engine_sympy.py /var/projects/_review_scratch/s15-build-r0-claude/check/eval_engine_sympy.stdout /var/projects/_review_scratch/s15-build-r0-claude/check/compare_wolfram_to_derive.stdout /var/projects/_review_scratch/s15-build-r0-grok/derive/derive_s15.py /var/projects/_review_scratch/s15-build-r0-grok/derive/derive_s15_stdout.txt /var/projects/_review_scratch/s15-build-r0-grok/evidence/py_normal_diff.txt /var/projects/_review_scratch/s15-build-r0-grok/evidence/wl_normal_diff.txt /var/projects/_review_scratch/s15-build-r0-grok/evidence/py_balance_diff_short.txt | sed 's|/var/projects/_review_scratch/||'
c2810f124e8a3f526fc3de617cabc4aea4df21c8ef6fdb9f44c8f2af3b004e02  s15-build-r0-claude/derive/derive_s15_v2.py
311c6c3274958482273d6a809c04184ae5945a828b207e1c66907fe764c82d19  s15-build-r0-claude/derive/derive_s15_v2.stdout
49e24587a413ffe51a3c393ca59f1cf221aacd21797038892b09f1bf085eb804  s15-build-r0-claude/check/eval_engine_sympy.py
3ce758d223685c476ec5f87831f4f3810c6a23bd319c16a0cc255ee306dfdc8f  s15-build-r0-claude/check/eval_engine_sympy.stdout
4044cdce51df6c7253d74e973fd7d5dc1551f2797bb8345d675471ec9f9eca91  s15-build-r0-claude/check/compare_wolfram_to_derive.stdout
07eaf2050fb2ad5b0977cef6845c08ac39f74a1bb01422fd8ce1046e5c7349ba  s15-build-r0-grok/derive/derive_s15.py
cf257f753645b8247752b3ca8ce01b72e68c50214f55ed8d25a042c53c8757e2  s15-build-r0-grok/derive/derive_s15_stdout.txt
534ee63fd91831ddcacb62827a2e39524e6ca3065f6d4bd2cabfaca5db09eff6  s15-build-r0-grok/evidence/py_normal_diff.txt
61ec6581c4035f2776857f671b2131f4fceb43bc6c1e37d62189474c2144589f  s15-build-r0-grok/evidence/wl_normal_diff.txt
71f9dc837ae9da64c86bcfca74863cefce0f64dbe78643904375533ca05a39fd  s15-build-r0-grok/evidence/py_balance_diff_short.txt
```
