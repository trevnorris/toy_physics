# Measurements: S15 reading build, review round 0 lookups

Generator: `_scratch/s15/gen/s15_reading_r0_lookups.sh` (GNU grep, sed, cut, cat, find, git, systemctl, sha256sum,
`LC_ALL=C`). Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts and
hashes. The legs' scripts were not re-run here; their filed stdout is retrieved. `C/` is `/var/projects/_review_scratch/s15-reading-r0-claude/`; `G/` is
`/var/projects/_review_scratch/s15-reading-r0-grok/`. Long lines are cut at a fixed column.

## M0. Baseline, prompt, legs finished, reports fixed, tree unchanged, no new files in the repository root

```
$ sha256sum -c _scratch/s15/reading_review/s15_reading_build_baseline_r0.sha256
research/pde_ledger_v3/scripts/S15_separation_departure.py: OK
research/pde_ledger_v3/mathematica/S15_separation_departure.wl: OK
research/pde_ledger_v3/scripts/out/S15_separation_departure.out: OK
research/pde_ledger_v3/mathematica/out/S15_separation_departure.out: OK
research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md: OK
```

```
$ sha256sum _scratch/s15/reading_review/s15_reading_build_review_prompt_r0.md
e82ea291cb5236bf6bbcec2b59d506efe891a42a6d21608918e4132ee7f6e0d9  _scratch/s15/reading_review/s15_reading_build_review_prompt_r0.md
```

```
$ cat /var/projects/_review_scratch/s15-reading-r0-grok/DONE; ls -la /var/projects/_review_scratch/s15-reading-r0-claude/REVIEW.md /var/projects/_review_scratch/s15-reading-r0-grok/REVIEW.md | cut -c24-
exit 0
 trevnorris 10859 Oct 28 18:30 /var/projects/_review_scratch/s15-reading-r0-claude/REVIEW.md
 trevnorris 11264 Oct 28 15:21 /var/projects/_review_scratch/s15-reading-r0-grok/REVIEW.md
```

```
$ (cd /var/projects/_review_scratch && sha256sum -c /var/projects/toy_physics/_scratch/s15/reading_review/r0_reports.sha256)
s15-reading-r0-claude/REVIEW.md: OK
s15-reading-r0-grok/REVIEW.md: OK
```

```
$ diff <(git status --porcelain | grep -v -F -e 'S15_reading_build_review_r0_lookups.md' -e 'S15_reading_build_review_r0_disposition.md'; git rev-parse HEAD) _scratch/s15/reading_review/status_before_r0.txt && echo 'status unchanged apart from the r0 disposition and lookups files'
status unchanged apart from the r0 disposition and lookups files
```

```
$ find . -maxdepth 1 -newer _scratch/s15/reading_review/s15_reading_build_review_prompt_r0.md ! -name .git | wc -l
0
```

```
$ systemctl --user is-active s11c-guard-44fc98082f0b.service; grep -E '"(unit|cwd)"|d01_general' _scratch/s15/reading-review-r0-guard/claude_d01/invocation.json; ls _scratch/s15/reading-review-r0-guard/claude_d01/
inactive
    "d01_general_profiles.py"
  "cwd": "/var/projects/_review_scratch/s15-reading-r0-claude",
  "unit": "s11c-guard-44fc98082f0b",
duration-validation.json
effective-limits.json
invocation.json
limit-validation.json
resource-samples.jsonl
stderr
stdout
```

## M1. The verdicts

```
$ grep -n -F '**Verdict' /var/projects/_review_scratch/s15-reading-r0-claude/REVIEW.md /var/projects/_review_scratch/s15-reading-r0-grok/REVIEW.md | sed 's|/var/projects/_review_scratch/||'
s15-reading-r0-claude/REVIEW.md:3:**Verdict: NEEDS REVISION**
s15-reading-r0-grok/REVIEW.md:3:**Verdict: NEEDS REVISION**
```

## M1a. The contract: Part B items 3 and 4 and the exactness paragraph (the reading list at 7f52b4d2)

```
$ git diff --quiet 7f52b4d2 -- research/pde_ledger_v3/directives/S15_reading_decision_list.md && echo 'list unchanged since 7f52b4d2'; sed -n '32,33p;75,82p;86,88p' research/pde_ledger_v3/directives/S15_reading_decision_list.md
list unchanged since 7f52b4d2
- Each engine constructs them again from the supplied equations. The `.wl` imports nothing; it may reuse the text of
  `mathematica/S15_two_drain_region.wl`. Neither engine reads the other engine's code or output.
3. `E_Q^geo` on the profile family, with its dependence on `d` explicit and exact. Where the integrand holds no general
   field besides the family's amplitudes, the region's quantities and the velocities, the integrals are evaluated in
   closed form. Otherwise the object is printed with the family substituted and left unevaluated.
4. For each `Q` whose item-3 object is in closed form:
   - the conditions under which that object vanishes for every `D` on the domain, exact, with every special exponent
     value treated. A condition is on the exponents and, where they enter, on the amplitudes, the radii and their
     rates, and the velocities. It is printed as equations and inequalities in the named quantities;
   - `Q` on the family under each condition.
Every object is exact on S15's domain under the conditions it carries: S15's pressure conditions, and for each
departure object a condition under which it exists, as above. No limit or expansion is taken in place of an exact
object. Where an engine does not obtain an object exactly, it prints a label saying so, beside what it did print.
```

## M2. F1, the engines: the item-3 route and its label

```
$ grep -n -F -e 'if not angular_polynomial' -e 'return S.Integral(f,*expr.limits)' research/pde_ledger_v3/scripts/S15_separation_departure.py
433:            if not angular_polynomial:
434:                return S.Integral(f,*expr.limits)
```

```
$ sed -n '484,487p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
        else:
            result=closed_route(family_geo)
            status='Closed form' if not result.has(S.Integral) else 'NOT OBTAINED IN CLOSED FORM: angular integrand is outside this engine\'s exact polynomial integration route.'
            if result.has(S.Integral): gaps.append((name,'item3',status))
```

```
$ sed -n '266,269p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200
 poly=PolynomialQ[angularInput,angularVars] && FreeQ[angularInput,theta|phi];
 answer=If[poly,
  Fold[Integrate[#1,Sequence@@{#2},Assumptions->(assumptions && Element[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}]],Reals]),GenerateConditions->True]&,integrand,SortBy[rangeList,If[First[#]===s,1,0]&]
  Inactive[Integrate][integrand,ranges]];
```

```
$ sed -n '307,312p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200
 familyResult=If[hasGeneral[familyGeo],familyGeo,angularRoute[familyGeo]];
 solo[readingTag<>"/item3/object",familyResult];
 item3Status=Which[hasGeneral[familyGeo],"Exact family integral with general fields retained, as permitted by item 3.",
  !FreeQ[familyResult,Inactive[Integrate]|Integrate],"NOT OBTAINED IN CLOSED FORM: non-polynomial angular factors remain outside this engine's exact integration route.",
  !FreeQ[familyResult,ConditionalExpression],"NOT OBTAINED ON ALL EXPONENT BRANCHES: Integrate returned a ConditionalExpression; the complement is unresolved.",
  True,"Closed form"];
```

## M3. F1, the transcripts: item-3 statuses, the gap index, and the d-dependent item-4 branch

```
$ grep -c -E '/item3/status: "?NOT OBTAINED IN CLOSED FORM' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -c -E '/item3/status: "?Closed form' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -c -E '/item3/status: "?Exact family' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -c -E '/item3/status: "?NOT OBTAINED ON ALL' research/pde_ledger_v3/scripts/out/S15_separation_departure.out
51
93
192
0
```

```
$ grep -c '^EXACT_OBJECT_GAP' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -c -F 'item4/equation_before_quantifier_elimination' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep -c -F '/item4/condition: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out
51
0
93
```

```
$ grep -c -E '/item3/status: "?NOT OBTAINED IN CLOSED FORM' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep -c -E '/item3/status: "?Closed form' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep -c -E '/item3/status: "?Exact family' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep -c -E '/item3/status: "?NOT OBTAINED ON ALL' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
48
96
192
0
```

```
$ grep -c '^EXACT_OBJECT_GAP' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep -c -F 'item4/equation_before_quantifier_elimination' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep -c -F '/item4/condition: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
48
0
96
```

```
$ grep -n -F 'reading/outer/transport_radius/component_1/labels_drain1/item3/' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-150
13108:reading/outer/transport_radius/component_1/labels_drain1/item3/family_integral: 2*Integral((-(N_1*d + X_2_1(t))*((N_1*d + X_2_1(t))**2 - 2*(N_1*
13109:reading/outer/transport_radius/component_1/labels_drain1/item3/family_integral/comparison: There is nothing to compare it with in this engine.
13110:reading/outer/transport_radius/component_1/labels_drain1/item3/existence_condition: SUFFICIENT existence condition (not necessary or complete): 
13111:reading/outer/transport_radius/component_1/labels_drain1/item3/object: 2*Integral((-(N_1*d + X_2_1(t))*((N_1*d + X_2_1(t))**2 - 2*(N_1*d + X_2_1
13112:reading/outer/transport_radius/component_1/labels_drain1/item3/object/comparison: There is nothing to compare it with in this engine.
13113:reading/outer/transport_radius/component_1/labels_drain1/item3/status: NOT OBTAINED IN CLOSED FORM: angular integrand is outside this engine's e
```

```
$ grep -n -F 'reading/outer/transport_radius/component_1/labels_{drain1}/item3/' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-150
12698:reading/outer/transport_radius/component_1/labels_{drain1}/item3/family_integral: 2*Inactive[Integrate][rout[t]^2*Sin[theta]*(-(Cos[theta]^2*((-
12699:reading/outer/transport_radius/component_1/labels_{drain1}/item3/family_integral/comparison: "There is nothing to compare it with in this engine
12700:reading/outer/transport_radius/component_1/labels_{drain1}/item3/existence_condition: "SUFFICIENT condition under which the departure exists, no
12701:reading/outer/transport_radius/component_1/labels_{drain1}/item3/object: 2*Inactive[Integrate][rout[t]^2*Sin[theta]*(-(Cos[theta]^2*((-(dist*dir
12702:reading/outer/transport_radius/component_1/labels_{drain1}/item3/object/comparison: "There is nothing to compare it with in this engine."
12703:reading/outer/transport_radius/component_1/labels_{drain1}/item3/status: "NOT OBTAINED IN CLOSED FORM: non-polynomial angular factors remain out
```

```
$ grep -n -F 'reading/region/component_1/rate_labels_drain2/item3/' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-150
14426:reading/region/component_1/rate_labels_drain2/item3/family_integral: 2*Integral(((R(t)*sin(theta)*cos(phi) + X_2_1(t))*((R(t)*cos(theta) + X_2_3
14427:reading/region/component_1/rate_labels_drain2/item3/family_integral/comparison: There is nothing to compare it with in this engine.
14428:reading/region/component_1/rate_labels_drain2/item3/existence_condition: SUFFICIENT existence condition (not necessary or complete): on an open 
14429:reading/region/component_1/rate_labels_drain2/item3/object: 2*Integral(((R(t)*sin(theta)*cos(phi) + X_2_1(t))*(R(t)**2)**(a_2/2)*(R(t)**2)**(b_2
14430:reading/region/component_1/rate_labels_drain2/item3/object/comparison: There is nothing to compare it with in this engine.
14431:reading/region/component_1/rate_labels_drain2/item3/status: NOT OBTAINED IN CLOSED FORM: angular integrand is outside this engine's exact polyno
```

```
$ grep -n -F 'reading/region/component_1/rate_labels_{drain2}/item3/' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-150
14028:reading/region/component_1/rate_labels_{drain2}/item3/family_integral: Inactive[Integrate][0, {s, rin[t], rout[t]}, {theta, 0, Pi}, {phi, 0, 2*P
14029:reading/region/component_1/rate_labels_{drain2}/item3/family_integral/comparison: "There is nothing to compare it with in this engine."
14030:reading/region/component_1/rate_labels_{drain2}/item3/existence_condition: "SUFFICIENT condition under which the departure exists, not a necessa
14031:reading/region/component_1/rate_labels_{drain2}/item3/object: 0
14032:reading/region/component_1/rate_labels_{drain2}/item3/object/comparison: "There is nothing to compare it with in this engine."
14033:reading/region/component_1/rate_labels_{drain2}/item3/status: "Closed form"
```

```
$ sed -n '1,3p;14,16p;22,24p' /var/projects/_review_scratch/s15-reading-r0-claude/lookups.stdout | cut -c1-160
========== py
L1 item-1 object tags 309  in accepted 309  byte-identical 309
L2 item3 status tally: {'Exact family object retained unevaluated: integrand contains a general field, as permitted by item 3.': 192, "NOT OBTAINED IN CLOSED FO
     3 reading/outer/transport_centre/component_i/labels_drain1+drain2
     3 reading/outer/transport_radius/component_i/labels_drain1
     3 reading/outer/transport_radius/component_i/labels_drain1+drain2
L4 closed-form objects 93  nonzero 6 ['reading/inner/transport_centre/component_1/labels_drain2', 'reading/inner/transport_centre/component_2/labels_drain2', 'r
L5 item4/equation_before_quantifier_elimination count 0
L6 item4 conditions 93  nontrivial 6 ['reading/inner/transport_centre/component_1/labels_drain2/item4/condition', 'reading/inner/transport_centre/component_2/la
```

## M4. F1, the legs: closed forms reached, checked against a second route

Residual lines and existence counts only: the legs' closed forms and conditions stay in their scratch files and are
not retrieved here.

```
$ grep -n -E '^(====|alpha=.*residual)' /var/projects/_review_scratch/s15-reading-r0-claude/d02_family.stdout | sed -n '1,8p' | cut -c1-140
1:==================== (a) carried group {drain1}, outer surface ====================
4:alpha=7/10: closed form e_y*q = -13.502477745768   quadrature = -13.502477745768   residual = 0
5:alpha=-1: closed form e_y*q = -1.3425266303147   quadrature = -1.3425266303147   residual = -2.22e-16
6:alpha=-2: closed form e_y*q = 0   quadrature = 3.6082248300318e-16   residual = -3.61e-16
7:alpha=-4: closed form e_y*q = 0.47602564013744   quadrature = 0.47602564013744   residual = -5e-16
8:alpha=0: closed form e_y*q = -5.6658951817358   quadrature = -5.6658951817358   residual = 0
9:alpha=-33/10: closed form e_y*q = 0.43281076733122   quadrature = 0.43281076733122   residual = -3.89e-16
10:alpha=1: closed form e_y*q = -19.340438636314   quadrature = -19.340438636314   residual = 0
```

```
$ grep -n -E '^(====|gamma=.*residual)' /var/projects/_review_scratch/s15-reading-r0-claude/d02_family.stdout | sed -n '2,11p' | cut -c1-140
40:==================== (b) region momentum group {drain1} ====================
43:gamma=-13/10: closed form e_z*phi = -1.15222140499  quadrature = -1.15222140499  residual = -5.15e-13
44:gamma=-2: closed form e_z*phi = -0.690044813038  quadrature = -0.690044813038  residual = -5.26e-13
45:gamma=-1: closed form e_z*phi = -1.45314880502  quadrature = -1.45314880502  residual = -4.81e-13
46:gamma=-3: closed form e_z*phi = -0.357084198626  quadrature = -0.357084198625  residual = -4.82e-13
47:gamma=-5: closed form e_z*phi = -0.124125890885  quadrature = -0.124125890885  residual = -3.51e-13
48:gamma=0: closed form e_z*phi = -3.30725959089  quadrature = -3.30725959089  residual = 3.99e-13
49:gamma=1: closed form e_z*phi = -8.04122055613  quadrature = -8.04122055613  residual = 2.37e-12
63:==================== (c) rate of region momentum, group {drain1} ====================
64:gamma=-13/10: time finite difference of P1_x = -1.32820205; assembled rate_x = -1.32820205; residual = -4.96e-11
```

```
$ grep -n -E '^(==========|  component)' /var/projects/_review_scratch/s15-reading-r0-claude/d03_more_closed_forms.stdout | cut -c1-140
1:========== (e) transport_centre {drain1} outer
6:  component 0: closed form -0.3787371481247  quadrature -0.3787371481247  residual -1.11e-16
7:  component 1: closed form -0.2419619894444  quadrature -0.2419619894444  residual 5.55e-17
8:  component 2: closed form 0.6452319718518  quadrature 0.6452319718518  residual -2.22e-16
12:========== (f) carried {drain1,drain2} outer
17:  component 0: closed form 19.7099477232  quadrature 19.7099477232  residual -1.07e-14
18:  component 1: closed form 0  quadrature -6.699386311061e-16  residual 6.7e-16
19:  component 2: closed form 0  quadrature -8.358010617365e-16  residual 8.36e-16
```

```
$ grep -E '^(Q_TR_OUTER_GROUP1_FORMULA|Q_TR_OUTER_GROUP1_PIECEWISE|EGEO_TR_GROUP1_GENERIC|P_REGION_GROUP1_GENERIC|Q_REGION_MOMENTUM_GROUP2):' /var/projects/_review_scratch/s15-reading-r0-grok/derive_representative.stdout | cut -d: -f1; grep -E '^(Q_TR_OUTER_GROUP1_FORMULA|Q_TR_OUTER_GROUP1_PIECEWISE|EGEO_TR_GROUP1_GENERIC|P_REGION_GROUP1_GENERIC):' /var/projects/_review_scratch/s15-reading-r0-grok/derive_representative.stdout | grep -c -F 'Integral'; grep -o -F 'Eq(p, -' /var/projects/_review_scratch/s15-reading-r0-grok/derive_representative.stdout | wc -l
Q_TR_OUTER_GROUP1_FORMULA
Q_TR_OUTER_GROUP1_PIECEWISE
EGEO_TR_GROUP1_GENERIC
Q_REGION_MOMENTUM_GROUP2
P_REGION_GROUP1_GENERIC
0
2
```

```
$ grep -E '^(Q_P0_RESIDUAL_ON_DOMAIN|SHELL_GAMMA2_RESIDUAL|P_GAMMA2_DIV_MINUS_DIRECT|DONE)' /var/projects/_review_scratch/s15-reading-r0-grok/verify_domain_identities.stdout
Q_P0_RESIDUAL_ON_DOMAIN: 0
SHELL_GAMMA2_RESIDUAL: -3.308722450212110699485635e-24
P_GAMMA2_DIV_MINUS_DIRECT: 0
DONE
```

```
$ cut -c1-140 /var/projects/_review_scratch/s15-reading-r0-claude/e01_check_py_engine.stdout
found tags: ['F1', 'F2', 'F3', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6']
point: d=2.26715680975 R=1.3 r_in=0.4 a1=-0.45 b1=-0.85 (gamma=a1+b1=-1.3, alpha=b1+2a1=-1.75)
[general carried comp 1] reference Q = 0.486935839714
   E_geo : engine -0.956359950434   reference -0.95635995101   residual 5.77e-10
   E_resp: engine -1.63823614775   reference -1.63823614507   residual -2.67e-09
[general carried comp 3] reference Q = 0.365201879785
   E_geo : engine -0.717269962825   reference -0.717269963258   residual 4.33e-10
   E_resp: engine -1.22867711081   reference -1.2286771088   residual -2.01e-09
[general transport comp 1] reference Q = -1.32145704278
   E_geo : engine -2.55541926292   reference -2.55541926293   residual 6.07e-12
   E_resp: engine -0.436150486471   reference -0.436150486648   residual 1.76e-10
[family carried drain1 comp 3] E_geo: engine -0.12258291583  reference -0.122582915536  residual -2.94e-10
[family momentum drain1 comp 3] E_geo: engine -1.21142067436  reference -1.21142067429  residual -6.84e-11
[family rate drain1 comp 1] E_geo: engine -0.332933707157  reference -0.332933707133  residual -2.4e-11
```

```
$ grep -v -E '^ *$|^ +-?[0-9]+ *$|NIntegrate|General::stop' /var/projects/_review_scratch/s15-reading-r0-claude/e02_check_wl_engine.stdout; grep -c -F 'NIntegrate::eincr' /var/projects/_review_scratch/s15-reading-r0-claude/e02_check_wl_engine.stdout
found: 9
reading/outer/carried/component_1/labels_{drain1}/E_geo -> -0.9563599504334704
reading/outer/carried/component_1/labels_{drain1}/E_resp -> -1.6382361477471432
reading/outer/carried/component_3/labels_{drain1}/E_geo -> -0.7172699628244179
reading/outer/carried/component_3/labels_{drain1}/E_resp -> -1.2286771107584031
reading/outer/sum_transport/component_1/E_geo -> -2.5554192629234
reading/outer/sum_transport/component_1/E_resp -> -0.4361504864711962
reading/outer/carried/component_3/labels_{drain1}/item3/object -> -0.12258291582950948
reading/region/component_3/momentum_labels_{drain1}/item3/object -> -1.2114203142791047
reading/region/component_1/rate_labels_{drain1}/item3/object -> -0.332933421835256
4
```

```
$ wc -c < /var/projects/_review_scratch/s15-reading-r0-claude/d02_family.stderr; wc -c < /var/projects/_review_scratch/s15-reading-r0-claude/d03_more_closed_forms.stderr; wc -c < /var/projects/_review_scratch/s15-reading-r0-claude/e01_check_py_engine.stderr
0
0
0
```

## M5. F2, the item-4 branch: engines, transcripts, run record, legs

```
$ sed -n '495,512p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
        elif result.has(d,*direction_D):
            emit(tag+'/item4/status','NOT OBTAINED EXACTLY: elimination of the universal separation and direction quantifiers is not supplied by this integration route.')
            solo(tag+'/item4/equation_before_quantifier_elimination',S.Eq(result,0,evaluate=False))
            solo(tag+'/item4/Q_on_family',family_q)
            gaps.append((name,'item4','Universal separation/direction elimination unresolved.'))
        else:
            # No D occurs, so the exact equation is already a condition in the
            # named parameters. Piecewise integration branches retain every
            # special exponent value; no solve/division discards a branch.
            condition=S.Eq(result,0,evaluate=False)
            solo(tag+'/item4/condition',condition)
            emit(tag+'/item4/condition_domain',(rout>rin,rin>0,*[S.Contains(a,S.S.Reals) for a in (*exponents_u,*exponents_rho)]))
            folded=S.piecewise_fold(result)
            solo(tag+'/item4/exponent_branches',folded)
            qclosed=closed_route(family_q) if not has_general_field(family_q) else family_q
            solo(tag+'/item4/Q_under_condition',qclosed)
            emit(tag+'/item4/Q_condition',condition)
            emit(tag+'/item4/status','Exact implicit equation in the named parameters; its Piecewise branches include exceptional exponent values. Q is printed restricted to that equation.')
```

```
$ sed -n '318,330p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-200
  If[!FreeQ[familyResult,dist|directionSeparation],
   solo[readingTag<>"/item4/equation_before_quantifier_elimination",familyResult==0];
   solo[readingTag<>"/item4/Q_on_family",familyQ];
   emit[readingTag<>"/item4/status","NOT OBTAINED EXACTLY: universal separation and direction elimination is unresolved."];
   AppendTo[exactGaps,{objectName,"item4","Universal separation/direction elimination unresolved."}],
   parameterCondition=familyResult==0;
   solo[readingTag<>"/item4/condition",parameterCondition];
   emit[readingTag<>"/item4/condition_domain",{rout[t]>rin[t]>0,Element[Join[Table[aa[k],{k,2}],Table[bb[k],{k,2}]],Reals]}];
   solo[readingTag<>"/item4/exponent_branches",PiecewiseExpand[familyResult]];
   conditionedQ=If[hasGeneral[familyQ],familyQ,angularRoute[familyQ]];
   solo[readingTag<>"/item4/Q_under_condition",conditionedQ];
   emit[readingTag<>"/item4/Q_condition",parameterCondition];
   emit[readingTag<>"/item4/status","Exact implicit equation in named parameters, retaining every Piecewise exponent branch; Q is restricted to the printed equation."]
```

```
$ grep -F 'reading/outer/transport_centre/component_1/labels_drain2/item' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | grep -E '/item(3/object|4/(condition|exponent_branches|Q_under_condition|status)): ' | cut -c1-230
reading/outer/transport_centre/component_1/labels_drain2/item3/object: 8*pi*(R(t)**2)**(a_2/2)*(R(t)**2)**(b_2/2)*A_2(t)*B_2(t)*R(t)**3*Derivative(X_2_1(t), t)/(3*sqrt(R(t)**2))
reading/outer/transport_centre/component_1/labels_drain2/item4/condition: Eq(8*pi*(R(t)**2)**(a_2/2)*(R(t)**2)**(b_2/2)*A_2(t)*B_2(t)*R(t)**3*Derivative(X_2_1(t), t)/(3*sqrt(R(t)**2)), 0)
reading/outer/transport_centre/component_1/labels_drain2/item4/exponent_branches: 8*pi*(R(t)**2)**(a_2/2)*(R(t)**2)**(b_2/2)*A_2(t)*B_2(t)*R(t)**3*Derivative(X_2_1(t), t)/(3*sqrt(R(t)**2))
reading/outer/transport_centre/component_1/labels_drain2/item4/Q_under_condition: 4*pi*(R(t)**2)**(a_2/2)*(R(t)**2)**(b_2/2)*A_2(t)*B_2(t)*R(t)**3*Derivative(X_2_1(t), t)/(3*sqrt(R(t)**2))
reading/outer/transport_centre/component_1/labels_drain2/item4/status: Exact implicit equation in the named parameters; its Piecewise branches include exceptional exponent values. Q is printed restricted to that equation.
```

```
$ grep -F 'reading/outer/transport_centre/component_1/labels_{drain2}/item' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -E '/item(3/object|4/(condition|exponent_branches|Q_under_condition|status)): ' | cut -c1-230
reading/outer/transport_centre/component_1/labels_{drain2}/item3/object: (8*Pi*rout[t]^(2 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3
reading/outer/transport_centre/component_1/labels_{drain2}/item4/condition: (8*Pi*rout[t]^(2 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3 == 0
reading/outer/transport_centre/component_1/labels_{drain2}/item4/exponent_branches: (8*Pi*rout[t]^(2 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3
reading/outer/transport_centre/component_1/labels_{drain2}/item4/Q_under_condition: (4*Pi*rout[t]^(2 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3
reading/outer/transport_centre/component_1/labels_{drain2}/item4/status: "Exact implicit equation in named parameters, retaining every Piecewise exponent branch; Q is restricted to the printed equation."
```

```
$ grep -c -F 'Piecewise' <(grep -F '/item4/exponent_branches: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out); grep -c -F 'Piecewise' <(grep -F '/item4/exponent_branches: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out)
0
0
```

```
$ grep -n -F 'Q restricted to the condition' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md | cut -c1-80; grep -o -F 'alongside Q restricted to the condition.' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md
47:General-field integrals retained under item 3’s explicit exception are not 
alongside Q restricted to the condition.
```

```
$ sed -n '105,106p' /var/projects/_review_scratch/s15-reading-r0-claude/d02_family.stdout | cut -c1-200
==================== (d) transport centre group {drain2}, outer ====================
Qt2 = Matrix([[4*pi*A2*B2*R**(a2 + b2 + 2)*Xd2x/3, 4*pi*A2*B2*R**(a2 + b2 + 2)*Xd2y/3, 4*pi*A2*B2*R**(a2 + b2 + 2)*Xd2z/3]]) ;  Egeo = 2*Qt2 (no X_1 dependence)
```

```
$ sed -n '22p;31p' /var/projects/_review_scratch/s15-reading-r0-grok/derive_representative.stdout | cut -c1-200
Q_TC_OUTER_GROUP2_PAR: 4*pi*A_2*B_2*R**(a_2 + b_2 + 2)*v2_par/3
EGEO_TRANSPORT_GROUP2_PAR: 8*pi*A_2*B_2*v2_par*(R**(a_2 + b_2 + 2) - r_in**(a_2 + b_2 + 2))/3
```

## M6. Reruns (examined, not filed)

```
$ cut -c1-150 /var/projects/_review_scratch/s15-reading-r0-claude/rerun_sha.stdout
651c200ffb36196fb3241e56dc9e9f85432d6f224e06c7664699577fbfcc56b7  runs/rerun/mathematica/out/S15_separation_departure.out
3220384e4b113fc029846aadd316e2fd3ee1d175d9a32a3d1c3f3eece81c5315  runs/rerun/scripts/out/S15_separation_departure.out
651c200ffb36196fb3241e56dc9e9f85432d6f224e06c7664699577fbfcc56b7  /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S15_separation_depa
3a0375e8ec2c8bbd4811cf32a26ecae5a090acb9796384a1edae36ba8e5d836b  /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S15_separation_departur
72d694f5cf9de6480127e6bd50a87c95d086052df7c8e3cdfb919b563caa23b7  /var/projects/toy_physics/research/pde_ledger_v3/scripts/S15_separation_departure.py
5c094e8e8d2e459c3cf5c1de996d46a75ce5c6dcce3d590d3fde720205e9a980  /var/projects/toy_physics/research/pde_ledger_v3/mathematica/S15_separation_departur
```

```
$ sed -n '1p' /var/projects/_review_scratch/s15-reading-r0-claude/e03_rerun_equivalence.stdout; grep -c -F 'residual 0' /var/projects/_review_scratch/s15-reading-r0-claude/e03_rerun_equivalence.stdout; grep -c -E 'residuals \[[^]]*[1-9]' /var/projects/_review_scratch/s15-reading-r0-claude/e03_rerun_equivalence.stdout
changed tags: 17
17
0
```

```
$ cat /var/projects/_review_scratch/s15-reading-r0-grok/rerun_residual.stdout /var/projects/_review_scratch/s15-reading-r0-grok/rerun_pressure_residual.stdout
OUTER_C1_ITEM3_EXPANDED 0
OUTER_C1_ITEM3_SIMPLIFIED 0
OUTER_C1_ITEM3_ON_DOMAIN 0
OUTER_C1_Q_EXPANDED 0
OUTER_C1_Q_SIMPLIFIED 0
OUTER_C1_Q_ON_DOMAIN 0
OUTER_C2_ITEM3_EXPANDED 0
OUTER_C2_ITEM3_SIMPLIFIED 0
OUTER_C2_ITEM3_ON_DOMAIN 0
INNER_C1_ITEM3_EXPANDED 0
INNER_C1_ITEM3_SIMPLIFIED 0
INNER_C1_ITEM3_ON_DOMAIN 0
DONE
LINE_3433_RESIDUAL 0
LINE_3465_RESIDUAL 0
LINE_9361_RESIDUAL 0
DONE
```

```
$ grep -n -F -e 'SHA-256' -e 'reproduc' research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md | cut -c1-160
15:- [Script](../../scripts/S15_separation_departure.py); SHA-256 `72d694f5cf9de6480127e6bd50a87c95d086052df7c8e3cdfb919b563caa23b7`.
16:- [Transcript](../../scripts/out/S15_separation_departure.out); SHA-256 `3a0375e8ec2c8bbd4811cf32a26ecae5a090acb9796384a1edae36ba8e5d836b`.
28:- [Script](../../mathematica/S15_separation_departure.wl); SHA-256 `5c094e8e8d2e459c3cf5c1de996d46a75ce5c6dcce3d590d3fde720205e9a980`.
29:- [Transcript](../../mathematica/out/S15_separation_departure.out); SHA-256 `651c200ffb36196fb3241e56dc9e9f85432d6f224e06c7664699577fbfcc56b7`.
```

## M7. Ablations (diff heads)

```
$ grep -n -E '^(== |[<>] )' /var/projects/_review_scratch/s15-reading-r0-claude/runs/ablation_diffs.txt | grep -v -F 'objects.items()' | cut -c1-170
1:== abl_py_corrupt_G
3:<     G = density*velocity
5:>     G = density**2*velocity
10:== abl_py_freeze_u
12:<     changing = [*u_heads,*rho_heads,ds[-1].func,*[q.func for q in vs[-1]],
14:>     changing = [*rho_heads,ds[-1].func,*[q.func for q in vs[-1]],
19:== abl_py_geo_holds_change
25:<         geo_derivative=variation(q,False)
27:>         geo_derivative=variation(q,True)
28:== abl_wl_corrupt_family
30:< profileDefinitions=Association[Join[Table[uu[k]->(AA[k][profileTime] profileRadius^aa[k]),{k,2}],Table[rr[k]->(BB[k][profileTime] profileRadius^bb[k]),{k,2}]]];
32:> profileDefinitions=Association[Join[Table[uu[k]->(AA[k][profileTime] profileRadius^aa[k] (1+profileRadius)),{k,2}],Table[rr[k]->(BB[k][profileTime] profileRadius^bb[
33:== abl_wl_freeze_X1dot
35:<   Cases[e,f:Derivative[ns__][h_][args___]/;MemberQ[changingHeads,h]:>f,{0,Infinity}],
36:<   Select[centreVelocityAtoms,!FreeQ[e,#]&]]];
38:>   Cases[e,f:Derivative[ns__][h_][args___]/;MemberQ[changingHeads,h]:>f,{0,Infinity}]]];
39:== abl_wl_move_X2
41:< positionForward=Thread[centres[[1]]->positionSymbols];
42:< positionBackward=Thread[positionSymbols->centres[[1]]];
44:> positionForward=Thread[centres[[2]]->positionSymbols];
45:> positionBackward=Thread[positionSymbols->centres[[2]]];
```

```
$ for f in /var/projects/_review_scratch/s15-reading-r0-claude/cmp_abl_*.stdout; do echo "$(basename $f): $(head -1 $f | cut -c1-140)"; done
cmp_abl_py_corrupt_G.stdout: tags: base 178  new 178
cmp_abl_py_corrupt_G_byobject.stdout: reading/outer/carried/component_1/labels_drain1/: identical tags: 43  changed tags: 0
cmp_abl_py_corrupt_G_item1.stdout: tags: base 1194  new 1194
cmp_abl_py_freeze_u.stdout: tags: base 185  new 185
cmp_abl_py_geo_holds_change.stdout: tags: base 178  new 172
cmp_abl_wl_corrupt_family.stdout: tags: base 15696  new 15696
cmp_abl_wl_freeze_X1dot.stdout: tags: base 15696  new 15696
cmp_abl_wl_move_X2.stdout: tags: base 15696  new 15696
```

```
$ cut -c1-170 /var/projects/_review_scratch/s15-reading-r0-grok/ablation_sy_form_diff.txt /var/projects/_review_scratch/s15-reading-r0-grok/ablation_wl_form_diff.txt /var/projects/_review_scratch/s15-reading-r0-grok/ablation_wl_family_diff.txt
TAG reading/outer/transport_centre/component_1/labels_drain2/item3/object:
ORIG L12063: 8*pi*(R(t)**2)**(a_2/2)*(R(t)**2)**(b_2/2)*A_2(t)*B_2(t)*R(t)**3*Derivative(X_2_1(t), t)/(3*sqrt(R(t)**2))
FORM L11991: 2*Integral(((R(t)*sin(theta)*cos(phi) + X_2_1(t))*((R(t)*cos(theta) + X_2_3(t))**2 - 2*(R(t)*cos(theta) + X_2_3(t))*X_2_3(t) + (R(t)*sin(phi)*sin(theta) + X_
TAG reading/outer/transport_centre/component_1/labels_drain2/item3/status:
ORIG L12065: Closed form
FORM L11993: Exact family object retained unevaluated: integrand contains a general field, as permitted by item 3.
TAG reading/outer/transport_centre/component_1/labels_drain2/E_resp:
ORIG L12033 LEN 12365 change_count 12 window **2 + X_2_2(t)**2 + X_2_3(t)**2), t)*u_2_change(sqrt((R(t)*cos(theta) + X_2_3(t))**2 - 2*(R(t)*cos(theta) + X_2_3(t))*X_
FORM L11961 LEN 1 change_count 0 window NO_CHANGE
TAG reading/outer/transport_centre/component_1/labels_{drain2}/item3/object
ORIG L11607:  (8*Pi*rout[t]^(2 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3
FORM L11535:  Inactive[Integrate][0, {theta, 0, Pi}, {phi, 0, 2*Pi}] + 2*Inactive[Integrate][rout[t]^2*Sin[theta]*(-(Cos[phi]*Sin[theta]*(Cos[theta]^2*rout[t]^2 + Cos[phi
TAG reading/outer/transport_centre/component_1/labels_{drain2}/item3/status
ORIG L11609:  "Closed form"
FORM L11537:  "Exact family integral with general fields retained, as permitted by item 3."
FORM_HAS_change False
FORM_HAS_centreVelocityChange False

FORM_OBJECT_CONTAINS change[rr[2]] (2 occurrences of change[). Window:
Inactive[Integrate][rout[t]^2*Sin[theta]*change[rr[2]][Sqrt[...], t]*...
ORIG L11607: (8*Pi*rout[t]^(2 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3
FAMILY L11607: (8*Pi*rout[t]^(3 + aa[2] + bb[2])*AA[2][t]*BB[2][t]*Derivative[1][xx[2, 1]][t])/3
```

## M8. The legs' reports, scripts and stdout (hashes)

```
$ (cd /var/projects/_review_scratch && sha256sum s15-reading-r0-claude/{REVIEW.md,d02_family.py,d02_family.stdout,d03_more_closed_forms.py,d03_more_closed_forms.stdout,e01_check_py_engine.py,e01_check_py_engine.stdout,e02_check_wl_engine.wl,e02_check_wl_engine.stdout,e03_rerun_equivalence.py,e03_rerun_equivalence.stdout,lookups.py,lookups.stdout,rerun_sha.stdout,runs/ablation_diffs.txt} s15-reading-r0-grok/{REVIEW.md,derive_representative.py,derive_representative.stdout,verify_domain_identities.py,verify_domain_identities.stdout,rerun_residual.py,rerun_residual.stdout,ablation_sy_form_diff.txt,ablation_wl_form_diff.txt,ablation_wl_family_diff.txt})
1edf6040123b34ae56a9e53f7d3308c652d2cd16b34fa9336fe83fb55550c1a7  s15-reading-r0-claude/REVIEW.md
933850c8936adfe750462d2910cdec760c87016cbeaf058942f18e41e8f7549f  s15-reading-r0-claude/d02_family.py
5c84293b91e95e3cce7e9b16ce5dcc3fa4a527a48dc6eeb37305a9ca772d5144  s15-reading-r0-claude/d02_family.stdout
e8d8d5a74cdc516b60a1783265b0fb19a8b3a08e08e6f01b0b4222b8b8c56388  s15-reading-r0-claude/d03_more_closed_forms.py
30dbaf33bbd910419df7d67120bb132b4bf7bf1e665134ad5c6f34006366bfc0  s15-reading-r0-claude/d03_more_closed_forms.stdout
d88932729d0843fb68b8a0f5a9ee5415a6125d32233c26e05d836d71cdb231b7  s15-reading-r0-claude/e01_check_py_engine.py
b77637b31e9373e71b7ef7f0b041b98030fc10701aac02b7e997efbf549ee833  s15-reading-r0-claude/e01_check_py_engine.stdout
73ad2f542fa93be03a3f2abc6638d52d4c3fb55803348abcf7172fad024b6fae  s15-reading-r0-claude/e02_check_wl_engine.wl
ead198d1f9de8651a329f744204dd1ab84445c0f8329e6b95f1319f4c798535d  s15-reading-r0-claude/e02_check_wl_engine.stdout
77614d5b68119b25d95a8d4b59a4830cb634af684add93a028de2fdfdc2a76ea  s15-reading-r0-claude/e03_rerun_equivalence.py
4d404851aa671322cf041149718b307d410fd464a144a0485f58d3ea0bdd88ff  s15-reading-r0-claude/e03_rerun_equivalence.stdout
fff691366cd39a6c05842b3d9343d8947e91b26ada3719f2e5cca268f6929920  s15-reading-r0-claude/lookups.py
58b6e5f71e3c301eaab8fbaea6055152846547d1bb45a30f612b89652b935145  s15-reading-r0-claude/lookups.stdout
0f0e9ccceaf73821e819b5736026c382681b33579883bb0b0a6f65bb629536e7  s15-reading-r0-claude/rerun_sha.stdout
467ec1e4741087223226f5a4a4d3023745db1df07a3cb46fc4bfb80fd83c7241  s15-reading-r0-claude/runs/ablation_diffs.txt
562ba695dd724a6d8e92061aeda24e9ba2cd89723a1992642e85a7e1c7a4560a  s15-reading-r0-grok/REVIEW.md
af6621d1ad9a0da96af3360cf509e0548e828d3e2623d5755d35f690e8f59d13  s15-reading-r0-grok/derive_representative.py
4bf2bc9936f5f55067a8af0ccc172ea0a8590dd540e531b06668d1b913665818  s15-reading-r0-grok/derive_representative.stdout
25e9c7c703b5f33df39e65a1b0f4c5403a4d265b6811d0cd55808d7b29a3e8d4  s15-reading-r0-grok/verify_domain_identities.py
ce386bea108c59c1202ae6d166722a2d4486ece27b7adac86ef6445994af2fd2  s15-reading-r0-grok/verify_domain_identities.stdout
a1917138afd0a1d601fb620632f552c0974aa23264b0931119af0a0c3b7b1947  s15-reading-r0-grok/rerun_residual.py
b70be3b279283f546aaf97ca50aa5fb949b6ce4668505aab0165c53798d422ef  s15-reading-r0-grok/rerun_residual.stdout
dd345f117c2f432c2276387ae29d8a47d3f3887f07dddb0f40f9d3da066998e3  s15-reading-r0-grok/ablation_sy_form_diff.txt
8da899d37a06cdfed27b8483255670c2988e13beace4b8fecc73ba88ee06ca0a  s15-reading-r0-grok/ablation_wl_form_diff.txt
c8578c2778b0aff996e2429fbfb75cd95d5c0e1c703d549136643bf4d43ffc15  s15-reading-r0-grok/ablation_wl_family_diff.txt
```
