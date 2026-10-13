# Measurements — S9b comparator build review, round 1 (generated 2026-10-13 12:43)

Generator: `_scratch/s9b_build/gen/s9b_comparator_r1_lookups.sh`. Leg evidence: `_scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/` (copied from `/tmp/rev_r1/`)
and `_scratch/s9b_build/preserved/rev_r1_grok_evidence/` (copied from `/tmp/s9b-review-r1/` and `/tmp/s9b-enum-r1/`). Operand values
the legs computed stay in that evidence. Line numbers are those of the reviewed comparator `_scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py`.

## The reviewed artifacts
```
$ sha256sum -c _scratch/s9b_build/s9b_comparator_build_review_baseline_r1.sha256
_scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_r1/test_S9b_cross_engine_comparator.py: OK
_scratch/s9b_build/comparator_r1/S9b_comparator_build_report.md: OK
_scratch/s9b_build/comparator_output/repair-01/full-03/comparator_source.py: OK
_scratch/s9b_build/comparator_output/repair-01/full-03/accounting.jsonl: OK
_scratch/s9b_build/comparator_output/repair-01/full-03/output_summary.json: OK
_scratch/s9b_build/comparator_output/repair-01/full-03/not_formed_inventory.json: OK
_scratch/s9b_build/comparator_output/repair-01/full-03/guard/outcome.json: OK
_scratch/s9b_build/comparator_output/repair-01/tests-07/comparator_source.py: OK
_scratch/s9b_build/comparator_output/repair-01/tests-07/test_source.py: OK
_scratch/s9b_build/comparator_output/repair-01/tests-07/guard/outcome.json: OK
research/pde_ledger_v3/directives/S9b_comparator_build_directive.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out: OK
research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out: OK
_scratch/s9b_build/s9b_comparator_build_review_prompt_r1.md: OK
_scratch/s9b_build/comparator_output/repair-01/full-03/comparison.jsonl: OK
```

## Both verdicts
```
$ head -1 _scratch/s9b_build/s9b_comparator_build_review_r1_claude.md
## Verdict: NEEDS REVISION
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_comparator_build_review_r1_grok.txt
Verdict: CLEAR
```

```
$ cat _scratch/s9b_build/s9b_comparator_build_review_r1_grok.DONE
EXIT=0
```

## R1 (Claude 1): the residual's precision is judged against the residual's own size
```
$ sed -n '906,910p' _scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py
            return {'undefined':str(exc)}
        with mp.workdps(50):
            discrepancy=abs(high-low)
            if not mp.isfinite(high) or discrepancy>mp.mpf('1e-14')*abs(high):
                return {'unresolved':'independent quadrature precision check did not establish 12 printed digits'}
```

```
$ sed -n '1038,1042p;1103,1104p' _scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py
        return {'coverage':'unresolved coverage'}
    if residual is not None:
        if 'undefined' in residual: return {'coverage':'undefined residual; comparison not formed'}
        if 'unresolved' in residual: return {'coverage':'unresolved coverage'}
        return {'coverage':'residual formed','residual':residual}
        lv=numeric_value(left,point,profile); rv=numeric_value(right,point,profile)
        sampled_residual=None if boolean else numeric_value(residual,point,profile)
```

```
$ grep -c 'independent quadrature precision check did not establish 12 printed digits' _scratch/s9b_build/comparator_output/repair-01/full-03/comparison.jsonl
61
```

```
$ head -1 _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/rows_agree_unresolved.stdout
rows with agreeing numeric operands whose residual is printed unresolved by the precision rule: 61 samples 1513
```

```
$ sed -n 2p _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/rows_agree_unresolved.stdout | tr ' ' '\n' | sed -E 's/^(F_FLOW_|F_LIVE_|R_[A-Z]+_ONLY_)//; s/_G[0-9]+$/_G*/' | sort | uniq -c
     18 A_DEFLECTION_G*
      7 A_ONE_WAY_ER_G*
      7 A_ONE_WAY_RE_G*
     10 A_RADAR_LOG_G*
      7 A_ROUND_TRIP_G*
      6 B_RESIDUAL_DEFLECTION
      6 B_RESIDUAL_RADAR
```

```
$ grep -n 'comparator rule\|operand-scaled rule' _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/discarded_residual.stdout | head -2
3:   comparator rule |high-low| > 1e-14*|high| ? True   (=> prints 'unresolved' and discards both numbers)
4:   operand-scaled rule |high-low| <= 1e-14*max|operand| ? True
```

```
$ grep -o '^sub-[a-z0-9-]* A_DEFLECTION_G100 on/value: .*' _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/sub-c1-coef.detail
sub-base A_DEFLECTION_G100 on/value: {'unresolved coverage': 27, 'formed zero': 6}
sub-c1-coef A_DEFLECTION_G100 on/value: {'formed NONZERO': 24, 'unresolved coverage': 3, 'formed zero': 6}
```

Grok's own scan of A_DEFLECTION_G001 (operand-match counts only):
```
$ grep -E '^(samples|match_but_not_formed|formed_and_equal|differ|other) ' _scratch/s9b_build/preserved/rev_r1_grok_evidence/s9b-review-r1/g001_agree_stdout.txt
samples 33
match_but_not_formed 16
formed_and_equal 6
differ 8
other 3
```

The round-0 requirement this precision check implements:
```
$ grep -o 'A numeric value prints only where the integral exists and the quadrature result is shown reliable at the printed precision.' research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_r0_review_disposition.md
A numeric value prints only where the integral exists and the quadrature result is shown reliable at the printed precision.
```

## R2 (Claude 2): b and c_0 share a value at every sample point
```
$ sed -n '690,696p' _scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py
POINTS=tuple({k:sp.Rational(v) for k,v in row.items()} for row in [
    dict(s9b_r='11/3',b='7/3',Z_E='101/3',Z_R='137/3',c_0='7/3',GM='1/101',
         Phi='2/103',K='5/7',n='2',m='3/2',s='1/2',rho_0='11/7',s9b_gamma='2/3',s9b_b_min='1/7'),
    dict(s9b_r='17/5',b='11/5',Z_E='151/5',Z_R='173/5',c_0='11/5',GM='2/107',
         Phi='3/109',K='7/11',n='3',m='5/3',s='2',rho_0='13/11',s9b_gamma='3/5',s9b_b_min='1/11'),
    dict(s9b_r='23/7',b='13/7',Z_E='211/7',Z_R='263/7',c_0='13/7',GM='3/113',
         Phi='5/127',K='11/13',n='1/2',m='7/5',s='-1',rho_0='17/13',s9b_gamma='5/7',s9b_b_min='1/13')])
```

```
$ cat _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/point_coincidences.stdout
p0: {'7/3': ['b', 'c_0']}
p1: {'11/5': ['b', 'c_0']}
p2: {'13/7': ['b', 'c_0']}
```

```
$ grep -n 'samples identical' _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/sub-c4-b2c0.detail
5:samples identical (py, wl, residual, comparison): 33 differing: 0
```

```
$ grep -n 'generic rational points' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
56:   generic rational points where both operands are defined (E1; numeric-PIT practice).
```

## R3 (Claude 3): compared_leaves counts leaves whose comparison formed nothing
```
$ sed -n '1097p;1117p;1229p;1353,1355p' _scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py
        return result,True
    return result,True
                    representation=structural_delta(a.ref.node,b.ref.node,self.bindings))
        row={'name':key,'accounting':'joined',
             **{e:{'parsed_leaves':parsed[e],'compared_leaves':len(used[i]),
                    'uncompared_leaves':parsed[e]-len(used[i])} for i,e in enumerate(('py','wl'))}}
```

```
$ grep 'TOTAL' _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/leaf_formed.stdout
TOTAL {'py_nothing_formed': 1476377, 'py_formed': 6611, 'wl_nothing_formed': 793835, 'wl_formed': 6098}
```

```
$ grep -n 'count actually compared' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
107:   count actually compared. A row compared below its parsed count stays visible.
```

## R4 (Claude 4): the Wolfram function-name maps have no repoint control
```
$ sed -n '439,445p;474p;537,538p' _scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py
WL_BUILTINS = {'Plus':'Add','Times':'Mul','Power':'Pow','Equal':'Equality',
               'Unequal':'Unequality','Greater':'StrictGreaterThan','Less':'StrictLessThan',
               'GreaterEqual':'GreaterThan','LessEqual':'LessThan', 'List':'Sequence',
               'Integrate':'Integral', 'Inactive/Integrate':'Integral',
               'Inactive/Limit':'Limit','Element':'Contains',
               'Inactive/ForAll':'ForAll','Inactive/Exists':'Exists',
               'ForAll':'ForAll','Exists':'Exists'}
        if wh in WL_BUILTINS or wh in ('And','Or','Not','Piecewise','ConditionalExpression','Rational','Min','Max','Abs','Sin','Cos','Tan','Sec','ArcCos','Log','Exp','Sqrt','Re','Im','densityPower'):
    if h in ('sin','cos','tan','log','exp','sqrt','Abs','Min','Max','re','im','acos','Sin','Cos','Tan','Log','Exp','Sqrt','Re','Im','ArcCos','Sec'):
        function = {'Sin':sp.sin,'Cos':sp.cos,'Tan':sp.tan,'Log':sp.log,'Exp':sp.exp,'Sqrt':sp.sqrt,'Re':sp.re,'Im':sp.im,'ArcCos':sp.acos,'Sec':sp.sec}.get(h,getattr(sp,h,None))
```

```
$ grep -n -o 'Any map between the two engines. symbol and function names' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
41:Any map between the two engines' symbol and function names
```

```
$ grep -n -o 'A \*\*repoint control\*\* must show that pointing a mapped name at a different object moves the' research/pde_ledger_v3/directives/S9b_comparator_build_directive.md
45:A **repoint control** must show that pointing a mapped name at a different object moves the
```

```
$ tail -3 _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/tests_fa1.stderr
Ran 46 tests in 61.985s

OK
```

```
$ grep -c 'res=NONZERO' _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/sub-fa1.summary _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/sub-base.summary
_scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/sub-fa1.summary:5
_scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/sub-base.summary:0
```

## O1 (Grok, checked and found sound): the Subs form of a profile derivative
```
$ grep -E '^(prod_zero_class|form_zero_class|form_changes_printed_residual|prod_tests_exit|form_tests_exit|Ran )' _scratch/s9b_build/preserved/rev_r1_grok_evidence/s9b-review-r1/ablate_stdout.txt
prod_zero_class unresolved prod_residual -Subs(Derivative(s9b_xi(_profile_argument), _profile_argument), _profile_argument, s9b_r) + Derivative(s9b_xi(s9b_r), s9b_r)
form_zero_class zero form_residual 0
form_changes_printed_residual True
prod_tests_exit 0
Ran 46 tests in 63.557s
form_tests_exit 0
Ran 46 tests in 63.245s
```

## O2 (Claude, name table): the velocity integration-constant row
```
$ sed -n '282,283p' _scratch/s9b_build/comparator_r1/S9b_cross_engine_comparator.py
    ('s9b_velocity_integration_constant', 'velocityIntegrationConstant',
     'velocity reduction constant; repair item 6', 470, 524, 'parameter'),
```

## Coverage facts routed to the record (Claude)
```
$ head -1 _scratch/s9b_build/preserved/rev_r1_claude_evidence/rev_r1/never_evaluated_rows.stdout
symbolic comparisons with samples: 171  with NO sample where both operands evaluated (excluding zero-profile families): 36
```

