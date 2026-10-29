# Measurements: S15 reading build review round 1, lookups

Generator: `_scratch/s15/gen/s15_reading_r1_lookups.sh` (GNU grep, sed, awk field selection, git, sha256sum,
`LC_ALL=C`). Each block is a command and its literal output: existence, verbatim retrieval, literal-match counts,
hashes. Computed results are quoted from the legs' saved stdout; their scripts were not re-run. The legs' closed
forms and conditions are not quoted here.

## M0. The round

```
$ sha256sum _scratch/s15/reading_review/s15_reading_build_review_prompt_r1.md
dc3e4cdee44428af6eec345854989e660fd345ae4f38998a149e5547b244644d  _scratch/s15/reading_review/s15_reading_build_review_prompt_r1.md
```

```
$ sha256sum -c _scratch/s15/reading_review/s15_reading_build_baseline_r1.sha256
research/pde_ledger_v3/scripts/S15_separation_departure.py: OK
research/pde_ledger_v3/mathematica/S15_separation_departure.wl: OK
research/pde_ledger_v3/scripts/out/S15_separation_departure.out: OK
research/pde_ledger_v3/mathematica/out/S15_separation_departure.out: OK
research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md: OK
```

```
$ cat _scratch/s15/reading_review/r1_reports.sha256; sha256sum -c _scratch/s15/reading_review/r1_reports.sha256
ac3c04ead2f8d5e206550f588629669b2e3e777cce73a148eef0b8fd44ddef96  /var/projects/_review_scratch/s15-reading-r1-grok/REVIEW.md
fb9a66c0e5c101803d96b4c2d11d2ca2c9e2294149823d8bcd097543d3dbeef4  /var/projects/_review_scratch/s15-reading-r1-claude/REVIEW.md
/var/projects/_review_scratch/s15-reading-r1-grok/REVIEW.md: OK
/var/projects/_review_scratch/s15-reading-r1-claude/REVIEW.md: OK
```

```
$ cat /var/projects/_review_scratch/s15-reading-r1-grok/DONE; wc -w /var/projects/_review_scratch/s15-reading-r1-claude/REVIEW.md /var/projects/_review_scratch/s15-reading-r1-grok/REVIEW.md
exit 0
 1753 /var/projects/_review_scratch/s15-reading-r1-claude/REVIEW.md
 1290 /var/projects/_review_scratch/s15-reading-r1-grok/REVIEW.md
 3043 total
```

```
$ git rev-parse HEAD; git status --short
5f0fc03a66226d9e95c18e307287cbc47c87bfbc
 M .claude/skills/prep-for-compact/SKILL.md
 M research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md
 M research/pde_ledger_v3/mathematica/S15_separation_departure.wl
 T research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
 M research/pde_ledger_v3/scripts/S15_separation_departure.py
 T research/pde_ledger_v3/scripts/out/S15_separation_departure.out
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_repair1_brief_lookups.md
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r1_disposition.md
?? research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r1_lookups.md
```

```
$ diff <(git status --short | grep -v -F S15_reading_build_review_r1_; git rev-parse HEAD) _scratch/s15/reading_review/status_before_r1.txt && echo 'status unchanged apart from the r1 disposition files'
status unchanged apart from the r1 disposition files
```

```
$ find . -path ./.git -prune -o -path ./_scratch -prune -o -newer _scratch/s15/reading_review/s15_reading_build_review_prompt_r1.md -type f -print | grep -v -F S15_reading_build_review_r1_
```

```
$ for j in claude_d01_a claude_e02 grok-d01; do echo "$j: $(tr -d ' \n' < _scratch/s11c/s15-reading-review-r1-guard/$j/outcome.json | cut -c1-150)"; done
claude_d01_a: {"exitCode":1,"wallSeconds":18044.661297265906,"unit":"s11c-guard-a954bd5b63c5","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
claude_e02: {"exitCode":1,"wallSeconds":13401.464203141979,"unit":"s11c-guard-faab7c080921","stderrBytes":0,"limitsVerified":true,"childOutcome":{"exitCode":1,"gu
grok-d01: {"exitCode":1,"wallSeconds":1336.5090611090418,"unit":"s11c-guard-8c28c192315a","stderrBytes":0,"limitsVerified":true,"childOutcome":null}
```

```
$ cat ~/.local/state/s11c-guard/resource-pool/reservations.json
{}
```

## M1. The contract and the round 0 obligations

```
$ git show 7f52b4d2:research/pde_ledger_v3/directives/S15_reading_decision_list.md | sed -n '75,82p;86,88p'
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

```
$ sed -n '82,90p;113,122p' research/pde_ledger_v3/directives/_measurements/S15_reading_build_review_r0_disposition.md
- **What must be true.**
  - Part B item 3 holds for every object in the gap index: each is printed in closed form, exact on the domain, with
    its dependence on `d` explicit and every special exponent value treated.
  - An object is labelled as not obtained only where that engine attempted it exactly and the attempt did not return
    it. The label says what was attempted.
  - Each engine reaches its closed forms by its own construction (Part B, "Supplied"; M1a).
  - No item-1, item-2 or item-5 object changes.

**V2. Item 4's condition is the closed object set equal to zero, and "Q under the condition" is the family's `Q`.**
- **What must be true.**
  - Part B item 4 holds for every closed form that item 3 prints, including those that depend on the separation. The
    conditions are reached by computation from the closed form and printed as equations and inequalities in the
    named quantities, and `Q` on the family is printed under each.
  - Where an engine does not obtain the conditions exactly after attempting them, it prints a label saying so beside
    what it did print (M1a).
  - No status line or run-record sentence states a property of a printed object that no computation in that engine
    produced.

## Examined and not filed
```

## M2. The engines at the r1 baseline

```
$ sed -n '771,788p;820p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-200
        attempted=[]
        for powers,coefficient in reduced_coefficients.items():
            emit(tag+'/scalar_identity_operand',(powers,coefficient))
            print(f'WORKING_ON: {current_object} / kernel abbreviation {powers}',flush=True)
            emit('WORKING_ON',current_object+' / kernel abbreviation '+str(powers))
            solver_operand,definitions=abbreviate_identity(coefficient,slots) if coefficient.has(d) else (coefficient,{})
            emit(tag+'/kernel_definitions',definitions)
            emit(tag+'/scalar_solver_operand',(powers,solver_operand))
            pair(tag+'/kernel_definition_reconstruction',coefficient,solver_operand.xreplace(definitions))
            emit(tag+'/kernel_definition_scope','Exact abbreviations with every free argument retained. Coefficient matching does not establish analytic independence of these kernels or classify their
            try:
                result_for_coefficient=cas_call('solve_undetermined_coeffs '+str(powers),S.solve_undetermined_coeffs,solver_operand,slots,d,dict=True)
            except (NotImplementedError,S.PolynomialError) as exc:
                result_for_coefficient=type(exc).__name__+': '+str(exc)
            attempted.append((powers,result_for_coefficient))
            emit(tag+'/scalar_identity_solver_return',(powers,result_for_coefficient))
        emit(tag+'/identity_solver_return',attempted)
        complete=(not reduced.has(d,*direction_D,S.Piecewise,S.appellf1,S.polylog))
        else: solutions=[]
```

```
$ grep -n -E 'result_for_coefficient|attempted' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-120
771:        attempted=[]
782:                result_for_coefficient=cas_call('solve_undetermined_coeffs '+str(powers),S.solve_undetermined_coeffs
784:                result_for_coefficient=type(exc).__name__+': '+str(exc)
785:            attempted.append((powers,result_for_coefficient))
786:            emit(tag+'/scalar_identity_solver_return',(powers,result_for_coefficient))
787:        emit(tag+'/identity_solver_return',attempted)
840:            reason='NOT OBTAINED EXACTLY: Unit-sphere monomial reduction, exact kernel abbreviation retaining all ar
892:            status='NOT OBTAINED IN CLOSED FORM: exact azimuth, squared-distance and Appell-kernel integration was a
```

```
$ sed -n '492p;497,512p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-220
 query=Implies[dist>rout[t]/.rules,directionAnswer];
     domain as an implication. Its unsupported return is kept, not promoted
     to conditions. No numerical separation or exponent is substituted. *)
  Block[{$MessageList={}},identityReturn=Quiet[trackedSolveAlways[query,dist]];messages=$MessageList];
  solo[prefix<>"/identity_solver_return",identityReturn];emit[prefix<>"/identity_solver_messages",messages];
  If[MatchQ[identityReturn,{{(_Rule|_RuleDelayed)...}...}] && messages==={},
   identityConditions=Or@@Map[And@@(Equal@@#& /@ (List@@@#))&,identityReturn];
   solo[prefix<>"/identity_candidate_conditions",identityConditions];
   (* SolveAlways supplies generic solutions; a returned list is not a proof
      of completeness on all exceptional parameter strata. Keep that gap. *)
   answer=identityReturn,
   answer=identityReturn]
 ];
 solo[prefix<>"/solver_return",answer];emit[prefix<>"/solver_messages",messages];
 If[!FreeQ[answer,Reduce|SolveAlways|ForAll|Exists] || !FreeQ[answer,dist] || ListQ[answer],
  emit[prefix<>"/status","NOT OBTAINED EXACTLY: indicator-polynomial direction elimination was computed; SolveAlways on the separation-domain implication, or Reduce for a separation-independent condition, did not return 
  AppendTo[exactGaps,{name,"item4","Indicator-polynomial direction elimination; SolveAlways on the separation-domain implication (Reduce for separation-independent conditions) returned unresolved."}],
```

## M3. Item 3 and item 4 in the transcripts

```
$ grep -E '^reading/.*/item3/status: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | sed 's/^[^:]*: //' | cut -c1-70 | sort | uniq -c
    144 Closed form from exact azimuth moments, squared-distance integration a
    192 Exact family object retained unevaluated: integrand contains a general
```

```
$ grep -E '^reading/.*/item3/status: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | sed 's/^[^:]*: //' | cut -c1-70 | sort | uniq -c
    144 "Closed form from generating-polynomial angular moments, distance subs
    192 "Exact family integral with general fields retained, as permitted by i
```

```
$ grep -E '^reading/.*/item4/status: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | sed 's/^[^:]*: //' | cut -c1-70 | sort | uniq -c
     96 Computed factor zero-set conditions and Q after their substitutions ar
     48 NOT OBTAINED EXACTLY: Unit-sphere monomial reduction, exact kernel abb
    192 Not applicable: the printed item-3 form retains a general field or is 
```

```
$ grep -E '^reading/.*/item4/status: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | sed 's/^[^:]*: //' | cut -c1-70 | sort | uniq -c
     48 "NOT OBTAINED EXACTLY: indicator-polynomial direction elimination was 
    192 "Not applicable to the general-field family integral."
     96 "Reduce conditions and Q simplified under each computed condition are 
```

```
$ grep -n -m1 -F 'reading/outer/carried/component_1/labels_drain1/item4/status' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | cut -c1-420
12967:reading/outer/carried/component_1/labels_drain1/item4/status: NOT OBTAINED EXACTLY: Unit-sphere monomial reduction, exact kernel abbreviation retaining all arguments, and solve_undetermined_coeffs on each scalar coefficient were attempted. Kernel identities and exceptional exponent strata remain unclassified; complete conditions and all corresponding restricted Q objects were not obtained.
```

```
$ grep -n -m1 -F 'reading/outer/carried/component_1/labels_{drain1}/item4/status' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | cut -c1-420
11548:reading/outer/carried/component_1/labels_{drain1}/item4/status: "NOT OBTAINED EXACTLY: indicator-polynomial direction elimination was computed; SolveAlways on the separation-domain implication, or Reduce for a separation-independent condition, did not return a complete exact parameter classification; SolveAlways lists, if present, are only candidate conditions without a completeness proof."
```

```
$ echo "py scalar_identity_solver_return lines: $(grep -c '/item4/scalar_identity_solver_return: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out)"; echo "  of which end ', None)': $(grep '/item4/scalar_identity_solver_return: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | grep -c ', None)$')"; echo "  of which in a labels_drain1 group: $(grep '/item4/scalar_identity_solver_return: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | grep -c 'labels_drain1')"; echo "  in a labels_drain1 group and ending ', None)': $(grep '/item4/scalar_identity_solver_return: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | grep 'labels_drain1' | grep -c ', None)$')"
py scalar_identity_solver_return lines: 258
  of which end ', None)': 162
  of which in a labels_drain1 group: 216
  in a labels_drain1 group and ending ', None)': 162
```

```
$ echo "py condition_union lines: $(grep -c '/item4/condition_union: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out); in a labels_drain1 group: $(grep '/item4/condition_union: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | grep -c 'labels_drain1')"
py condition_union lines: 96; in a labels_drain1 group: 54
```

```
$ echo "wl identity_solver_return lines: $(grep -c '/item4/identity_solver_return: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out); containing 'SolveAlways[Implies[': $(grep '/item4/identity_solver_return: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -c -F 'SolveAlways[Implies[')"; echo "wl identity_candidate_conditions lines: $(grep -c '/item4/identity_candidate_conditions' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out)"; echo "wl /item4/condition lines in a {drain1} group: $(grep '/item4/condition: ' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -c -F 'drain1')"
wl identity_solver_return lines: 48; containing 'SolveAlways[Implies[': 48
wl identity_candidate_conditions lines: 0
wl /item4/condition lines in a {drain1} group: 54
```

```
$ grep -c 'EXACT_OBJECT_GAP' research/pde_ledger_v3/scripts/out/S15_separation_departure.out; grep 'EXACT_OBJECT_GAP' research/pde_ledger_v3/scripts/out/S15_separation_departure.out | grep -o -E "'item[34]'|\"item[34]\"" | sort | uniq -c
48
     48 'item4'
```

```
$ grep -c 'EXACT_OBJECT_GAP' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out; grep 'EXACT_OBJECT_GAP' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out | grep -o -E '"item[34]"' | sort | uniq -c
48
     48 "item4"
```

## M4. F1: the legs' evidence (probes in full; derivations as residuals and existence counts only)

```
$ cat /var/projects/_review_scratch/s15-reading-r1-claude/f01_solver_probe.stdout | cut -c1-220
A*(p+1)*d**p  [d**p abbreviated]: operand A*(p + 1)*condition_kernel_0(d, p); solve_undetermined_coeffs return: None
A*(p+1)*d**p  [d**p abbreviated]: engine gate complete=(not operand.has(d,...)) -> False
A*(p+1)*d**p  [not abbreviated]: operand A*d**p*(p + 1); solve_undetermined_coeffs return: None
A*(p+1)*d**p  [not abbreviated]: engine gate complete=(not operand.has(d,...)) -> False
A*(p+1)/d**2 + A*(p+1)*s/d**3  [rational in d]: operand A*(p + 1)/d**2 + A*s*(p + 1)/d**3; solve_undetermined_coeffs return: [{A: 0}, {p: -1}]
A*(p+1)/d**2 + A*(p+1)*s/d**3  [rational in d]: engine gate complete=(not operand.has(d,...)) -> False
A*(p+1)*K(d,p) + A*(p+1)*s/d**3  [mixed]: operand A*(p + 1)*condition_kernel_0(d, p) + A*s*(p + 1)/d**3; solve_undetermined_coeffs return: None
A*(p+1)*K(d,p) + A*(p+1)*s/d**3  [mixed]: engine gate complete=(not operand.has(d,...)) -> False
```

```
$ cat /var/projects/_review_scratch/s15-reading-r1-claude/f02_solvealways_probe.stdout | cut -c1-260
Implies, polynomial in dist: HoldForm[SolveAlways[Implies[dist > r, a*((p + 1)/dist^2) == 0], dist]] -> SolveAlways[Implies[dist > r, (a*(1 + p))/dist^2 == 0], dist]
equation, polynomial in dist: HoldForm[SolveAlways[a*((p + 1)/dist^2) == 0, dist]] -> {{a -> 0}, {p -> -1}}
Implies, Boole stratum: HoldForm[SolveAlways[Implies[dist > r, a*(Boole[p == -2]/dist^2) + a*((p + 1)/dist^3) == 0], dist]] -> SolveAlways[Implies[dist > r, (a*(1 + p))/dist^3 + (a*Boole[p == -2])/dist^2 == 0], dist]
Implies, power of dist: HoldForm[SolveAlways[Implies[dist > r, a*(p + 1)*dist^p == 0], dist]] -> SolveAlways[Implies[dist > r, a*dist^p*(1 + p) == 0], dist]
```

```
$ ls -la /var/projects/_review_scratch/s15-reading-r1-claude/d02_derive.py /var/projects/_review_scratch/s15-reading-r1-claude/d02_derive.stdout /var/projects/_review_scratch/s15-reading-r1-claude/d03_quadrature.py /var/projects/_review_scratch/s15-reading-r1-claude/d03_quadrature.stdout | awk '{print $5, $9}'; echo "EXPONENT_SET lines in d02: $(grep -c 'EXPONENT_SET' /var/projects/_review_scratch/s15-reading-r1-claude/d02_derive.stdout)"; grep -o 'max|res|: [^ ]*' /var/projects/_review_scratch/s15-reading-r1-claude/d03_quadrature.stdout | sort | uniq -c
10674 /var/projects/_review_scratch/s15-reading-r1-claude/d02_derive.py
21899 /var/projects/_review_scratch/s15-reading-r1-claude/d02_derive.stdout
13268 /var/projects/_review_scratch/s15-reading-r1-claude/d03_quadrature.py
15633 /var/projects/_review_scratch/s15-reading-r1-claude/d03_quadrature.stdout
EXPONENT_SET lines in d02: 6
      1 max|res|: 1.21e-13
      1 max|res|: 1.31e-15
      1 max|res|: 1.33e-14
      1 max|res|: 1.35e-15
      1 max|res|: 1.37e-14
      1 max|res|: 1.44e-15
      1 max|res|: 1.67e-15
      1 max|res|: 1.71e-16
      1 max|res|: 1.73e-18
      1 max|res|: 1.78e-14
      2 max|res|: 1.78e-15
      1 max|res|: 1.80e-16
      1 max|res|: 1.87e-16
      1 max|res|: 1.88e-15
      2 max|res|: 1.95e-14
      1 max|res|: 2.08e-17
      1 max|res|: 2.13e-14
      1 max|res|: 2.16e-15
      1 max|res|: 2.22e-15
      1 max|res|: 2.22e-16
      1 max|res|: 2.29e-16
      1 max|res|: 2.39e-17
      1 max|res|: 2.57e-18
      1 max|res|: 2.66e-15
      1 max|res|: 2.78e-15
      1 max|res|: 2.86e-17
      1 max|res|: 3.12e-17
      1 max|res|: 3.20e-14
      1 max|res|: 3.22e-15
      1 max|res|: 3.30e-17
      1 max|res|: 3.33e-16
      1 max|res|: 3.47e-18
      1 max|res|: 3.55e-14
      1 max|res|: 3.55e-15
      1 max|res|: 3.66e-15
      2 max|res|: 3.77e-15
      1 max|res|: 3.89e-15
      1 max|res|: 3.89e-16
      1 max|res|: 3.90e-18
      1 max|res|: 4.16e-16
      1 max|res|: 4.44e-15
      1 max|res|: 4.50e-17
      1 max|res|: 4.70e-16
      1 max|res|: 4.72e-16
      1 max|res|: 4.97e-14
      2 max|res|: 5.00e-16
      1 max|res|: 5.01e-16
      1 max|res|: 5.33e-15
      1 max|res|: 6.11e-16
      2 max|res|: 7.11e-15
      1 max|res|: 7.33e-15
      1 max|res|: 7.59e-16
      1 max|res|: 7.99e-15
      1 max|res|: 8.33e-16
      1 max|res|: 8.53e-14
      1 max|res|: 8.67e-16
      1 max|res|: 8.88e-15
      2 max|res|: 8.88e-16
```

```
$ ls -la /var/projects/_review_scratch/s15-reading-r1-grok/d02_conditions.py /var/projects/_review_scratch/s15-reading-r1-grok/d03_vanishing.py /var/projects/_review_scratch/s15-reading-r1-grok/d06_numeric.py /var/projects/_review_scratch/s15-reading-r1-grok/d07_specials.py /var/projects/_review_scratch/s15-reading-r1-grok/d08_wolfram_numeric.wl /var/projects/_review_scratch/s15-reading-r1-grok/d08b_wolfram_specials.wl | awk '{print $5, $9}'; grep -h -o 'residual: .*' /var/projects/_review_scratch/s15-reading-r1-grok/d03_vanishing.stdout /var/projects/_review_scratch/s15-reading-r1-grok/d06_numeric.stdout /var/projects/_review_scratch/s15-reading-r1-grok/d07_specials.stdout /var/projects/_review_scratch/s15-reading-r1-grok/d08_wolfram_numeric.stdout /var/projects/_review_scratch/s15-reading-r1-grok/d08b_wolfram_specials.stdout | sort | uniq -c
13624 /var/projects/_review_scratch/s15-reading-r1-grok/d02_conditions.py
4618 /var/projects/_review_scratch/s15-reading-r1-grok/d03_vanishing.py
10017 /var/projects/_review_scratch/s15-reading-r1-grok/d06_numeric.py
8345 /var/projects/_review_scratch/s15-reading-r1-grok/d07_specials.py
2040 /var/projects/_review_scratch/s15-reading-r1-grok/d08_wolfram_numeric.wl
1712 /var/projects/_review_scratch/s15-reading-r1-grok/d08b_wolfram_specials.wl
     24 residual: 0
      1 residual: 0``65.15703876762568
      1 residual: 0``65.94954503691338
      1 residual: 0``67.36553031334692
      1 residual: 0``68.77751658564723
      2 residual: 1.31208517726e-141
      1 residual: 2.14972035442e-137
```

## M5. F2: the Wolfram region rate {drain1}

```
$ grep -n -o '^reading/region/component_[123]/rate_labels_{drain1}/item3/[A-Za-z_]*' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
21454:reading/region/component_1/rate_labels_{drain1}/item3/family_integral
21455:reading/region/component_1/rate_labels_{drain1}/item3/family_integral
21456:reading/region/component_1/rate_labels_{drain1}/item3/existence_condition
21458:reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand
21459:reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand
21461:reading/region/component_1/rate_labels_{drain1}/item3/object
21462:reading/region/component_1/rate_labels_{drain1}/item3/object
21463:reading/region/component_1/rate_labels_{drain1}/item3/status
22849:reading/region/component_2/rate_labels_{drain1}/item3/family_integral
22850:reading/region/component_2/rate_labels_{drain1}/item3/family_integral
22851:reading/region/component_2/rate_labels_{drain1}/item3/existence_condition
22853:reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand
22854:reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand
22856:reading/region/component_2/rate_labels_{drain1}/item3/object
22857:reading/region/component_2/rate_labels_{drain1}/item3/object
22858:reading/region/component_2/rate_labels_{drain1}/item3/status
24219:reading/region/component_3/rate_labels_{drain1}/item3/family_integral
24220:reading/region/component_3/rate_labels_{drain1}/item3/family_integral
24221:reading/region/component_3/rate_labels_{drain1}/item3/existence_condition
24223:reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand
24224:reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand
24226:reading/region/component_3/rate_labels_{drain1}/item3/object
24227:reading/region/component_3/rate_labels_{drain1}/item3/object
24228:reading/region/component_3/rate_labels_{drain1}/item3/status
```

```
$ for n in 21458 21461 22853 22856 24223 24226; do L=$(sed -n "${n}p" research/pde_ledger_v3/mathematica/out/S15_separation_departure.out); echo "line $n $(echo "$L" | cut -d: -f1)"; for f in '(2 + aa[1] + bb[1])' '(3 + aa[1] + bb[1])' '(4 + aa[1] + bb[1])' '(5 + aa[1] + bb[1])'; do echo "  literal $f: $(echo "$L" | grep -o -F "$f" | wc -l)"; done; echo '  every equality in the line:'; echo "$L" | grep -o -E '[^,{]{0,40}== -?[0-9/]+' | sort | uniq -c; done
line 21458 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand
  literal (2 + aa[1] + bb[1]): 4
  literal (3 + aa[1] + bb[1]): 4
  literal (4 + aa[1] + bb[1]): 2
  literal (5 + aa[1] + bb[1]): 2
  every equality in the line:
      6  -1/2 + aa[1]/2 + bb[1]/2 == -1
      2  1/2 + aa[1]/2 + bb[1]/2 == -1
      5  aa[1] + bb[1] == -1
      3  aa[1] + bb[1] == -3
line 21461 reading/region/component_1/rate_labels_{drain1}/item3/object
  literal (2 + aa[1] + bb[1]): 20
  literal (3 + aa[1] + bb[1]): 20
  literal (4 + aa[1] + bb[1]): 10
  literal (5 + aa[1] + bb[1]): 10
  every equality in the line:
     10  -1/2 + aa[1]/2 + bb[1]/2 == -1
      2  1/2 + aa[1]/2 + bb[1]/2 == -1
     18  aa[1] + bb[1] == -1
      9  aa[1] + bb[1] == -3
line 22853 reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand
  literal (2 + aa[1] + bb[1]): 4
  literal (3 + aa[1] + bb[1]): 4
  literal (4 + aa[1] + bb[1]): 2
  literal (5 + aa[1] + bb[1]): 2
  every equality in the line:
      6  -1/2 + aa[1]/2 + bb[1]/2 == -1
      2  1/2 + aa[1]/2 + bb[1]/2 == -1
      5  aa[1] + bb[1] == -1
      3  aa[1] + bb[1] == -3
line 22856 reading/region/component_2/rate_labels_{drain1}/item3/object
  literal (2 + aa[1] + bb[1]): 20
  literal (3 + aa[1] + bb[1]): 20
  literal (4 + aa[1] + bb[1]): 10
  literal (5 + aa[1] + bb[1]): 10
  every equality in the line:
     10  -1/2 + aa[1]/2 + bb[1]/2 == -1
      2  1/2 + aa[1]/2 + bb[1]/2 == -1
     18  aa[1] + bb[1] == -1
      9  aa[1] + bb[1] == -3
line 24223 reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand
  literal (2 + aa[1] + bb[1]): 4
  literal (3 + aa[1] + bb[1]): 4
  literal (4 + aa[1] + bb[1]): 2
  literal (5 + aa[1] + bb[1]): 2
  every equality in the line:
      6  -1/2 + aa[1]/2 + bb[1]/2 == -1
      2  1/2 + aa[1]/2 + bb[1]/2 == -1
      5  aa[1] + bb[1] == -1
      3  aa[1] + bb[1] == -3
line 24226 reading/region/component_3/rate_labels_{drain1}/item3/object
  literal (2 + aa[1] + bb[1]): 20
  literal (3 + aa[1] + bb[1]): 20
  literal (4 + aa[1] + bb[1]): 10
  literal (5 + aa[1] + bb[1]): 10
  every equality in the line:
     10  -1/2 + aa[1]/2 + bb[1]/2 == -1
      2  1/2 + aa[1]/2 + bb[1]/2 == -1
     18  aa[1] + bb[1] == -1
      9  aa[1] + bb[1] == -3
```

```
$ sed -n '21463p;22858p;24228p' research/pde_ledger_v3/mathematica/out/S15_separation_departure.out
reading/region/component_1/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular moments, distance substitution and finite-interval Appell integration, including logarithmic branches."
reading/region/component_2/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular moments, distance substitution and finite-interval Appell integration, including logarithmic branches."
reading/region/component_3/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular moments, distance substitution and finite-interval Appell integration, including logarithmic branches."
```

```
$ grep -E '^SUMMARY' /var/projects/_review_scratch/s15-reading-r1-claude/e10_scan_strata.stdout; grep -i 'NOT EVALUABLE' /var/projects/_review_scratch/s15-reading-r1-claude/e10_scan_strata.stdout | grep -v '^SUMMARY' | cut -c1-200
SUMMARY WL: strata tested 606; not evaluable 24
SUMMARY PY: strata tested 84; not evaluable 0
WL reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 5 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-93/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 4 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-73/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 4 + (2*(+ ((-1)/2) + (aa[1]/2) + (bb[1]/2)))) -> stratum a1=-7/20 b1=-53/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivi
WL reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 2 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-33/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_1/rate_labels_{drain1}/item3/object: divisor (+ 5 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-93/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_1/rate_labels_{drain1}/item3/object: divisor (+ 4 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-73/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_1/rate_labels_{drain1}/item3/object: divisor (+ 4 + (2*(+ ((-1)/2) + (aa[1]/2) + (bb[1]/2)))) -> stratum a1=-7/20 b1=-53/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_1/rate_labels_{drain1}/item3/object: divisor (+ 2 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-33/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 5 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-93/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 4 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-73/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 4 + (2*(+ ((-1)/2) + (aa[1]/2) + (bb[1]/2)))) -> stratum a1=-7/20 b1=-53/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivi
WL reading/region/component_2/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 2 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-33/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/object: divisor (+ 5 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-93/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/object: divisor (+ 4 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-73/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/object: divisor (+ 4 + (2*(+ ((-1)/2) + (aa[1]/2) + (bb[1]/2)))) -> stratum a1=-7/20 b1=-53/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_2/rate_labels_{drain1}/item3/object: divisor (+ 2 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-33/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 5 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-93/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 4 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-73/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 4 + (2*(+ ((-1)/2) + (aa[1]/2) + (bb[1]/2)))) -> stratum a1=-7/20 b1=-53/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivi
WL reading/region/component_3/rate_labels_{drain1}/item3/Q_closed_operand: divisor (+ 2 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-33/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/object: divisor (+ 5 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-93/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/object: divisor (+ 4 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-73/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/object: divisor (+ 4 + (2*(+ ((-1)/2) + (aa[1]/2) + (bb[1]/2)))) -> stratum a1=-7/20 b1=-53/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
WL reading/region/component_3/rate_labels_{drain1}/item3/object: divisor (+ 2 + aa[1] + bb[1]) -> stratum a1=-7/20 b1=-33/20 a2=-3/5 b2=2/5: NOT EVALUABLE ZeroDivisionError
```

```
$ head -6 /var/projects/_review_scratch/s15-reading-r1-claude/e05_locate_zero.stdout | cut -c1-200
E6 a1=-1 b1=-2 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: zero divisions 4; distinct 2
E6 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand:   divisor x2: ((+ 2 + aa[1] + bb[1])*(+ 3 + aa[1] + bb[1])*(+ dist + rin[t])*(+ rin[t] - rout[t])^2)
E6 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand:     inside: ['Piecewise default (conditions not met: (+ aa[1] + bb[1]) == (-1))']
E6 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand:   divisor x2: ((+ 2 + aa[1] + bb[1])*(+ 3 + aa[1] + bb[1])*(+ dist - rin[t])*(+ rin[t] - rout[t])^2)
E6 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand:     inside: ['Piecewise default (conditions not met: (+ aa[1] + bb[1]) == (-1))']
E7 a1=-2 b1=-3 reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: zero divisions 2; distinct 2
```

```
$ awk '{print $1, $2, $3, $4, $NF}' /var/projects/_review_scratch/s15-reading-r1-claude/e09_limit.stdout
F1 a1=-1/2 b1=-3/2+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/object: -3.186e-6
F1 a1=-1/2 b1=-3/2+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: -3.101e-9
F1 a1=-1/2 b1=-3/2+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 3.271e-9
E6 a1=-1 b1=-2+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/object: -7.993e-7
E6 a1=-1 b1=-2+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 2.093e-9
E6 a1=-1 b1=-2+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 3.697e-9
F2 a1=-3/2 b1=-5/2+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/object: 1.504e-7
F2 a1=-3/2 b1=-5/2+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 1.212e-8
F2 a1=-3/2 b1=-5/2+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 1.185e-8
F5 a1=-1 b1=-4+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/object: 4.715e-7
F5 a1=-1 b1=-4+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 3.302e-8
F5 a1=-1 b1=-4+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/object: 3.214e-8
F1 a1=-1/2 b1=-3/2+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -1.154e-6
F1 a1=-1/2 b1=-3/2+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -1.154e-9
F1 a1=-1/2 b1=-3/2+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: 1.154e-9
E6 a1=-1 b1=-2+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -6.893e-7
E6 a1=-1 b1=-2+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -6.893e-10
E6 a1=-1 b1=-2+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: 6.893e-10
F2 a1=-3/2 b1=-5/2+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -4.005e-7
F2 a1=-3/2 b1=-5/2+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -4.005e-10
F2 a1=-3/2 b1=-5/2+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: 4.005e-10
F5 a1=-1 b1=-4+(1/100000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -2.334e-7
F5 a1=-1 b1=-4+(1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: -2.334e-10
F5 a1=-1 b1=-4+(-1/100000000) reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand: 2.334e-10
```

```
$ grep -E '^SUMMARY' /var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout /var/projects/_review_scratch/s15-reading-r1-claude/e03b_compare_wl_final.stdout /var/projects/_review_scratch/s15-reading-r1-claude/e03b_compare_py.stdout /var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout
/var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout:SUMMARY WL-vs-WLPY: compared 286; max scaled residual 6.483e-16
/var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout:SUMMARY WLPY-vs-REF: compared 326; max scaled residual 8.812e-09; above 1e-7 or non-numeric: 4
/var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout:SUMMARY WLPY-vs-REF:   ('E6', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout:SUMMARY WLPY-vs-REF:   ('E7', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout:SUMMARY WLPY-vs-REF:   ('E6', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e06_compare_wlpy_final.stdout:SUMMARY WLPY-vs-REF:   ('E7', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e03b_compare_wl_final.stdout:SUMMARY PY: compared 330; max scaled residual 3.258e-08; above 1e-7 or non-numeric: 0
/var/projects/_review_scratch/s15-reading-r1-claude/e03b_compare_wl_final.stdout:SUMMARY WL: compared 286; max scaled residual 1.749e-15; above 1e-7 or non-numeric: 0
/var/projects/_review_scratch/s15-reading-r1-claude/e03b_compare_py.stdout:SUMMARY PY: compared 330; max scaled residual 3.258e-08; above 1e-7 or non-numeric: 0
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY PY: compared 210; max scaled residual 3.258e-08; above 1e-7 or non-numeric: 0
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY: compared 196; max scaled residual 7.879e-09; above 1e-7 or non-numeric: 14
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F1', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F2', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F3', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F4', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F5', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('E6', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('E7', 'reading/region/component_1/rate_labels_{drain1}/item3/object')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F1', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F2', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F3', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F4', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('F5', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('E6', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
/var/projects/_review_scratch/s15-reading-r1-claude/e08_compare_strata.stdout:SUMMARY WLPY:   ('E7', 'reading/region/component_1/rate_labels_{drain1}/item3/Q_closed_operand')
```

```
$ grep -F 'reading/region/component_1/rate_labels_{drain1}/item3/status' /var/projects/_review_scratch/s15-reading-r1-claude/runs/abl_W1/mathematica/out/S15_separation_departure.out
reading/region/component_1/rate_labels_{drain1}/item3/status: "Closed form from generating-polynomial angular moments, distance substitution and finite-interval Appell integration, including logarithmic branches."
```

```
$ head -12 /var/projects/_review_scratch/s15-reading-r1-claude/c05_abl_W1.stdout
W1: tags baseline 26588 ablated 26213
W1:         item1    same 336
W1:         item2    same 5712
W1:         item3   moved 96
W1:         item3    same 2208
W1:         item4   moved 342
W1:         item4    same 6282
W1:         item5    same 3360
W1:         item6    same 2016
W1: reading_other    same 672
W1:      toplevel missing 232
W1:      toplevel   moved 1
```

## M6. Items 1, 2, 5 and 6 against round 0

```
$ head -1 /var/projects/_review_scratch/s15-reading-r1-claude/c01_compare_tags.stdout; grep -E '^(PY|WL): +(item[0-9]|item3_kernel|inner|outer|dimension_argument|partition|family) ' /var/projects/_review_scratch/s15-reading-r1-claude/c01_compare_tags.stdout
PY: tags r0 17659 r1 27588 duplicates r0 1258 r1 6638
PY: dimension_argument     same 832
PY:         family     same 7
PY:          inner     same 570
PY:          item3   differ 213
PY:          item3  only_r1 288
PY:          item3     same 1803
PY:   item3_kernel  only_r1 2297
PY:          item4   differ 387
PY:          item4  only_r0 651
PY:          item4  only_r1 6043
PY:          item4     same 528
PY:          outer     same 570
PY:      partition     same 27
WL:         family     same 7
WL:          inner     same 570
WL:          item3   differ 198
WL:          item3  only_r1 288
WL:          item3     same 1818
WL:   item3_kernel  only_r1 2550
WL:          item4   differ 486
WL:          item4  only_r0 384
WL:          item4  only_r1 5424
WL:          item4     same 714
WL:          outer     same 570
```

```
$ grep -E '^(PY|WL): +[^ ]+ +(differ|only_r0|only_r1) ' /var/projects/_review_scratch/s15-reading-r1-claude/c01_compare_tags.stdout | grep -v -E 'WORKING_ON|#dup' | cut -c1-140
PY:          item3   differ 213
PY:          item3  only_r1 288
PY:   item3_kernel  only_r1 2297
PY:          item4   differ 387
PY:          item4  only_r0 651
PY:          item4  only_r1 6043
PY:  reading_other   differ 12
PY:       toplevel   differ 241
PY:       toplevel  only_r0 146
PY:       toplevel  only_r1 2098
WL:          item3   differ 198
WL:          item3  only_r1 288
WL:   item3_kernel  only_r1 2550
WL:          item4   differ 486
WL:          item4  only_r0 384
WL:          item4  only_r1 5424
WL:   item4_kernel  only_r1 14
WL:       toplevel   differ 120
WL:       toplevel  only_r0 282
WL:       toplevel  only_r1 1515
```

```
$ awk -F'; ' '{print $1"; "$NF}' /var/projects/_review_scratch/s15-reading-r1-claude/c02_pressure_value.stdout | cut -c1-160
reading/inner/pressure/component_1/whole_minus_expansion/E_geo: chars r0 3382 r1 3383; residual of integrands 0
reading/inner/pressure/component_1/whole_minus_expansion/E_resp: chars r0 6166 r1 6167; residual of integrands 0
reading/inner/pressure/component_2/whole_minus_expansion/E_geo: chars r0 3382 r1 3383; residual of integrands 0
reading/inner/pressure/component_2/whole_minus_expansion/E_resp: chars r0 6166 r1 6167; residual of integrands 0
reading/inner/pressure/component_3/whole_minus_expansion/E_geo: chars r0 3346 r1 3347; residual of integrands 0
reading/inner/pressure/component_3/whole_minus_expansion/E_resp: chars r0 6130 r1 6131; residual of integrands 0
reading/outer/pressure/component_1/whole_minus_expansion/E_geo: chars r0 3178 r1 3177; residual of integrands 0
reading/outer/pressure/component_1/whole_minus_expansion/E_resp: chars r0 5794 r1 5793; residual of integrands 0
reading/outer/pressure/component_2/whole_minus_expansion/E_geo: chars r0 3178 r1 3177; residual of integrands 0
reading/outer/pressure/component_2/whole_minus_expansion/E_resp: chars r0 5794 r1 5793; residual of integrands 0
reading/outer/pressure/component_3/whole_minus_expansion/E_geo: chars r0 3142 r1 3141; residual of integrands 0
reading/outer/pressure/component_3/whole_minus_expansion/E_resp: chars r0 5758 r1 5757; residual of integrands 0
```

```
$ grep -E '^(sympy|wolfram) ' /var/projects/_review_scratch/s15-reading-r1-grok/item1256_byte_compare.txt
sympy r0=6384 r1=6384 only_r0=0 only_r1=0 differ=12 same=6372
wolfram r0=7056 r1=7056 only_r0=0 only_r1=0 differ=0 same=7056
```

## M7. Controls

```
$ for a in P1 P3; do head -14 /var/projects/_review_scratch/s15-reading-r1-claude/c05_abl_$a.stdout; done
P1: tags baseline 27588 ablated 27648
P1:         item1    same 336
P1:         item2   moved 6
P1:         item2    same 5706
P1:         item3   moved 114
P1:         item3    same 2190
P1:         item4   moved 1340
P1:         item4     new 54
P1:         item4    same 5618
P1:         item5    same 2688
P1:         item6    same 2017
P1: reading_other    same 672
P1:      toplevel    same 4591
P1: moved/new tag stems 360 (first 60 by count)
P3: tags baseline 27588 ablated 35462
P3:         item1    same 336
P3:         item2   moved 825
P3:         item2    same 4887
P3:         item3   moved 330
P3:         item3    same 1974
P3:         item4   moved 1708
P3:         item4     new 6012
P3:         item4    same 5250
P3:         item5    same 2688
P3:         item6    same 2017
P3: reading_other    same 672
P3:      toplevel   moved 3
P3:      toplevel    same 4588
```

```
$ head -14 /var/projects/_review_scratch/s15-reading-r1-claude/c05_abl_P2_partial.stdout
P2_partial: tags baseline 27588 ablated 2766
P2_partial:         item1 missing 334
P2_partial:         item1    same 2
P2_partial:         item2 missing 5678
P2_partial:         item2    same 34
P2_partial:         item3 missing 2290
P2_partial:         item3   moved 2
P2_partial:         item3    same 12
P2_partial:         item4 missing 6928
P2_partial:         item4   moved 13
P2_partial:         item4    same 17
P2_partial:         item5 missing 2672
P2_partial:         item5    same 16
P2_partial:         item6 missing 2005
```

```
$ for d in grok-ablate-py grok-ablate-wl; do test -f _scratch/s11c/s15-reading-review-r1-guard/$d/outcome.json && echo "$d: $(tr -d ' \n' < _scratch/s11c/s15-reading-review-r1-guard/$d/outcome.json | cut -c1-120)"; done; ls _scratch/s11c/s15-reading-review-r1-guard | grep -i abl
grok-ablate-py: {"exitCode":0,"wallSeconds":2441.5913411640795,"unit":"s11c-guard-badf93357aeb","stderrBytes":0,"limitsVerified":true,"c
grok-ablate-wl: {"exitCode":0,"wallSeconds":234.10491017298773,"unit":"s11c-guard-5c25038cf25b","stderrBytes":0,"limitsVerified":true,"c
claude_abl_P1
claude_abl_P2
claude_abl_P2_b
claude_abl_P2_c
claude_abl_P3
claude_abl_W1
grok-ablate-py
grok-ablate-wl
```

## M8. Reproducibility

```
$ sha256sum /var/projects/_review_scratch/s15-reading-r1-claude/runs/rerun/scripts/out/S15_separation_departure.out /var/projects/_review_scratch/s15-reading-r1-claude/runs/rerun/mathematica/out/S15_separation_departure.out /var/projects/_review_scratch/s15-reading-r1-grok/rerun_py/out/S15_separation_departure.out
8115416bb4d21ac338101a3a2dcb8284b597d63e6d2c7d70bc94a287637db70a  /var/projects/_review_scratch/s15-reading-r1-claude/runs/rerun/scripts/out/S15_separation_departure.out
80265713a3146f82b1ecc9972712947b81bad845de6c4782026a349a5439efad  /var/projects/_review_scratch/s15-reading-r1-claude/runs/rerun/mathematica/out/S15_separation_departure.out
8115416bb4d21ac338101a3a2dcb8284b597d63e6d2c7d70bc94a287637db70a  /var/projects/_review_scratch/s15-reading-r1-grok/rerun_py/out/S15_separation_departure.out
```

```
$ head -3 /var/projects/_review_scratch/s15-reading-r1-claude/c03_rerun_py_tags.stdout; head -3 /var/projects/_review_scratch/s15-reading-r1-claude/c06_rerun_wl_tags.stdout
RERUN_PY: tags r0 27588 r1 27588 duplicates r0 6638 r1 6638
RERUN_PY:           P5_G     same 1
RERUN_PY:           P6_T     same 1
RERUN_WL: tags r0 26588 r1 26588 duplicates r0 4228 r1 4228
RERUN_WL: P5_momentum_density     same 1
RERUN_WL:      P6_stress     same 1
```

```
$ cat /var/projects/_review_scratch/s15-reading-r1-grok/d09_rerun_residual.stdout | cut -c1-160
reading/inner/pressure/component_1/whole_minus_expansion/E_geo: 0
reading/inner/pressure/component_1/whole_minus_expansion/E_resp: 0
reading/inner/pressure/component_2/whole_minus_expansion/E_geo: 0
reading/inner/pressure/component_2/whole_minus_expansion/E_resp: 0
reading/inner/pressure/component_3/whole_minus_expansion/E_geo: 0
reading/inner/pressure/component_3/whole_minus_expansion/E_resp: 0
reading/outer/pressure/component_1/whole_minus_expansion/E_geo: 0
reading/outer/pressure/component_1/whole_minus_expansion/E_resp: 0
reading/outer/pressure/component_2/whole_minus_expansion/E_geo: 0
reading/outer/pressure/component_2/whole_minus_expansion/E_resp: 0
reading/outer/pressure/component_3/whole_minus_expansion/E_geo: 0
reading/outer/pressure/component_3/whole_minus_expansion/E_resp: 0
DONE
```

```
$ for p in power_primitive_derivative kernel_definition_reconstruction rerun reproduc; do echo "$p in the run record: $(grep -c -i "$p" research/pde_ledger_v3/directives/_measurements/S15_separation_departure_run_record.md)"; done
power_primitive_derivative in the run record: 0
kernel_definition_reconstruction in the run record: 0
rerun in the run record: 0
reproduc in the run record: 0
```

```
$ echo "kernel_definition_reconstruction/residual lines: $(grep -c '/kernel_definition_reconstruction/residual: ' research/pde_ledger_v3/scripts/out/S15_separation_departure.out); equal to 0: $(grep -c '/kernel_definition_reconstruction/residual: 0$' research/pde_ledger_v3/scripts/out/S15_separation_departure.out)"
kernel_definition_reconstruction/residual lines: 258; equal to 0: 258
```

```
$ sed -n '499p;501,507p' research/pde_ledger_v3/mathematica/S15_separation_departure.wl | cut -c1-120; sed -n '776p' research/pde_ledger_v3/scripts/S15_separation_departure.py | cut -c1-140
  Block[{$MessageList={}},identityReturn=Quiet[trackedSolveAlways[query,dist]];messages=$MessageList];
  If[MatchQ[identityReturn,{{(_Rule|_RuleDelayed)...}...}] && messages==={},
   identityConditions=Or@@Map[And@@(Equal@@#& /@ (List@@@#))&,identityReturn];
   solo[prefix<>"/identity_candidate_conditions",identityConditions];
   (* SolveAlways supplies generic solutions; a returned list is not a proof
      of completeness on all exceptional parameter strata. Keep that gap. *)
   answer=identityReturn,
   answer=identityReturn]
            solver_operand,definitions=abbreviate_identity(coefficient,slots) if coefficient.has(d) else (coefficient,{})
```
