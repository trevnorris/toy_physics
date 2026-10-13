# Measurements — S9b comparator build review, round 0 (generated 2026-10-13 06:54)

Generator: `_scratch/s9b_build/gen/s9b_comparator_r0_lookups.sh`. Leg evidence (files under 2 MB): `_scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/` (from
`/tmp/s9b_cmp_review_r0_claudeleg/`) and `_scratch/s9b_build/s9b_comparator_build_review_r0_grok_evidence/` (from `/tmp/s9b_review_r0_*`). Operand values the legs computed
stay in that evidence.

## The reviewed artifacts
```
$ sha256sum -c _scratch/s9b_build/s9b_comparator_build_review_baseline_r0.sha256
research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py: OK
research/pde_ledger_v3/directives/_measurements/S9b_comparator_build_report.md: OK
_scratch/s9b_build/comparator_output/full-04/comparator_source.py: OK
_scratch/s9b_build/comparator_output/full-04/accounting.jsonl: OK
_scratch/s9b_build/comparator_output/full-04/output_summary.json: OK
_scratch/s9b_build/comparator_output/full-04/guard/outcome.json: OK
_scratch/s9b_build/comparator_output/tests-08/test_source.py: OK
_scratch/s9b_build/comparator_output/tests-08/outcome.json: OK
research/pde_ledger_v3/directives/S9b_comparator_build_directive.md: OK
research/pde_ledger_v3/directives/S9b_SHARED_PHYSICS.md: OK
research/pde_ledger_v3/scripts/S9b_light_bending_sympy_audit.py: OK
research/pde_ledger_v3/mathematica/S9b_light_bending_mathematica_audit.wl: OK
research/pde_ledger_v3/scripts/out/S9b_light_bending_sympy_audit.out: OK
research/pde_ledger_v3/mathematica/out/S9b_light_bending_mathematica_audit.out: OK
_scratch/s9b_build/s9b_comparator_build_review_prompt_r0.md: OK
_scratch/s9b_build/comparator_output/full-04/comparison.jsonl: OK
```

```
$ grep -o '"sha256": "[0-9a-f]*"' _scratch/s9b_build/comparator_output/full-04/output_summary.json
"sha256": "bf312dc2648ba75f2476838e935fb3859c04288535c341296f3125118310f501"
"sha256": "1bf34cbc891a4882151a596804be972a467c56ee4414cc2afc3e6770be1cabbe"
```

## Both verdicts
```
$ head -1 _scratch/s9b_build/s9b_comparator_build_review_r0_claude.md
**Verdict: NEEDS REVISION**
```

```
$ grep -o 'Verdict[*:]* *[A-Z][A-Z ]*' _scratch/s9b_build/s9b_comparator_build_review_r0_grok.txt
Verdict: NEEDS REVISION
```

```
$ cat _scratch/s9b_build/s9b_comparator_build_review_r0_grok.DONE
EXIT=0
```

## D1 (Claude 1, Grok 2): undefined samples on the G121 rows
```
$ grep -n 'Subs(\*args).doit(deep=False)\|has(sp.nan,sp.zoo)' research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
711:            return sp.Subs(*args).doit(deep=False)
741:    if value.has(sp.nan,sp.zoo): return {'undefined':math_text(value)}
```

```
$ grep -n "numeric operands printed" research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
829:    return {'coverage':'unresolved' if 'unresolved' in a or 'unresolved' in b else 'numeric operands printed'}
```

```
$ grep -o "^ *[0-9]* ('[A-Z_0-9]*G121', 'py undefined', 'wl approx')" _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/out/undefined_census.stdout
    30 ('A_ONE_WAY_ER_G121', 'py undefined', 'wl approx')
    30 ('A_ONE_WAY_RE_G121', 'py undefined', 'wl approx')
    30 ('A_ROUND_TRIP_G121', 'py undefined', 'wl approx')
    37 ('F_FLOW_A_ONE_WAY_ER_G121', 'py undefined', 'wl approx')
    37 ('F_FLOW_A_ONE_WAY_RE_G121', 'py undefined', 'wl approx')
    37 ('F_FLOW_A_ROUND_TRIP_G121', 'py undefined', 'wl approx')
    37 ('F_LIVE_A_ONE_WAY_ER_G121', 'py undefined', 'wl approx')
    37 ('F_LIVE_A_ONE_WAY_RE_G121', 'py undefined', 'wl approx')
    37 ('F_LIVE_A_ROUND_TRIP_G121', 'py undefined', 'wl approx')
```

```
$ grep -c 'zoo' _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/out/undefined_text.stdout
5
```

## D2 (Claude 2): quadrature with no convergence or singularity check
```
$ grep -n 'isfinite(high)' research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
757:        if not mp.isfinite(high): return {'undefined':str(high)}
```

```
$ grep -n 'samples with both operands finite approx' _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/out/undefined_census.stdout
12:brane_zero_at_point samples with both operands finite approx and |residual| > 1e-10: 18
```

```
$ grep -n "^ *[0-9]* ('brane_zero_at_point', '[a-z]*')" _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/out/precision_census.stdout
3:    78 ('brane_zero_at_point', 'wl')
4:    69 ('brane_zero_at_point', 'py')
6:    19 ('brane_zero_at_point', 'residual')
```

## D3 (Claude 3): mutation controls
```
$ sed -n '34,45p;61,71p;259,266p' research/pde_ledger_v3/scripts/test_S9b_cross_engine_comparator.py
def products(value):
    """Remove operands, sample values, diagnostics and raw serialization.

    Keep only mathematical products, comparison outcomes and child topology.
    No echoed operand (including a digest) can satisfy a mutation control.
    """
    if isinstance(value,list): return [products(x) for x in value]
    if not isinstance(value,dict): return value
    keys={'outcome','comparison','children','unpaired_metadata','residual','boolean_residual',
          'zero_class','samples','boolean_difference','same_exact_representation','coverage'}
    if 'outcome' in value or 'configuration' in value:
        return {k:products(v) for k,v in value.items() if k in keys}
    def mutation(self,py,wl,changed,*,channel='leaf'):
        before=self.fixture(py,wl); after=self.fixture(py,changed)
        self.assertEqual([r['py'] for r in before[2] if r['kind']=='operands'],
                         [r['py'] for r in after[2] if r['kind']=='operands'])
        self.assertNotEqual(products(before[0]),products(after[0]),
                            'mutation must move a computed comparison product')
        # An operand echo or parse change cannot rescue a disabled comparison.
        with patch.object(c.Comparator,channel,return_value=({'outcome':'not_formed'},False)):
            with self.assertRaises(AssertionError):
                self.assertNotEqual(products(self.fixture(py,wl)[0]),products(self.fixture(py,changed)[0]))
        return before,after
    def test_full_configuration_set_preserved(self):
        py='Integral(Function("s9b_delta")(Symbol("s9b_r")),Tuple(Symbol("s9b_r"),Integer(2),Integer(3)))'
        wl='Inactive[Integrate][deltaProfile[r],{r,2,3}]'
        comparison=self.fixture(py,wl)[0][0]
        rows=comparison['samples']
        expected={f'p{i}/{profile}' for i in range(len(c.POINTS)) for profile in c.PROFILE_NAMES}
        self.assertEqual({r['configuration'] for r in rows},expected)
        self.assertTrue(all('py' in r and 'wl' in r and 'residual' in r for r in rows))
```

```
$ tail -3 _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/ablate/copies/E1/tests.stderr
Ran 32 tests in 10.654s

FAILED (failures=4)
```

```
$ tail -3 _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/ablate/copies/F9/tests.stderr
Ran 32 tests in 17.214s

OK
```

```
$ tail -3 _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/ablate/copies/S1/tests.stderr
Ran 32 tests in 15.607s

OK
```

```
$ grep -n 'computed-products identical' _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/ablate/diff_ablated_corrupt.stdout
2:A_DEFLECTION_G000: computed-products identical; operand/residual-text echo identical; accounting DIFFERS
3:A_DEFLECTION_G001: computed-products identical; operand/residual-text echo differs; accounting same
70:BRANCH_EXISTENCE: computed-products identical; operand/residual-text echo identical; accounting same
72:B_CONDITION_DEFLECTION: computed-products identical; operand/residual-text echo identical; accounting same
75:B_GAMMA_DIFFERENCE: computed-products identical; operand/residual-text echo identical; accounting same
253:F_MASS_SOLUTION: computed-products identical; operand/residual-text echo identical; accounting same
339:A_DEFLECTION_G000: computed-products identical; operand/residual-text echo identical; accounting DIFFERS
408:BRANCH_EXISTENCE: computed-products identical; operand/residual-text echo identical; accounting same
410:B_CONDITION_DEFLECTION: computed-products identical; operand/residual-text echo identical; accounting same
413:B_GAMMA_DIFFERENCE: computed-products identical; operand/residual-text echo identical; accounting same
591:F_MASS_SOLUTION: computed-products identical; operand/residual-text echo identical; accounting same
```

## D4 (Claude 4): the B_GAMMA_DIFFERENCE gate
```
$ sed -n '108,111p;570,573p' research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
        if h == 'Tuple':
            if a and all(v.head == 'Sequence' and len(v.args) == 2
                         and v.args[0].head == 'Text' for v in a) and len({v.args[0].value for v in a}) == len(a):
                return record([(v.args[0].value, v.args[1]) for v in a])
                                  'applicability':extract(ref.child(1),engine)})
        if len(a)==2 and all(x.head=='Sequence' and len(x.args)==2 and attached(ref.child(i).child(0)) for i,x in enumerate(a)):
            return View(ref,role,{name: View(ref.child(i), 'case', {
                'domain':extract(ref.child(i).child(0),engine),
```

```
$ grep -c 'GATE NOT RECOGNIZED' _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/structure/gate_census.stdout
6
```

```
$ grep -o 'GATE NOT RECOGNIZED: [A-Z_]*' _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/structure/gate_census.stdout
GATE NOT RECOGNIZED: B_GAMMA_DIFFERENCE
GATE NOT RECOGNIZED: F_FLOW_B_GAMMA_DIFFERENCE
GATE NOT RECOGNIZED: F_LIVE_B_GAMMA_DIFFERENCE
GATE NOT RECOGNIZED: R_FLOW_ONLY_B_GAMMA_DIFFERENCE
GATE NOT RECOGNIZED: R_SPEED_ONLY_B_GAMMA_DIFFERENCE
GATE NOT RECOGNIZED: R_TILT_ONLY_B_GAMMA_DIFFERENCE
```

## D5 (Grok 1): strata coverage of an empty piece
```
$ sed -n '892,897p' research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
        for i,(_,_,piece) in enumerate(pieces):
            if i not in covered and numeric_value(piece,point,profile).get('truth') is True: covered.add(i)
    result['samples']=samples
    if pieces:
        result['strata']=[{'engine':eng,'piece':label,'condition':math_text(p),
                           'coverage':'sampled' if i in covered else 'unresolved coverage'} for i,(eng,label,p) in enumerate(pieces)]
```

```
$ grep -n 'PRINTED path strata\|SYNTHETIC sample0' _scratch/s9b_build/s9b_comparator_build_review_r0_grok_evidence/s9b_review_r0_recompute_pass2.stdout
51:PRINTED path strata [('py', 'condition', 'False', 'unresolved coverage'), ('py', 'complement', 'True', 'sampled'), ('wl', 'condition', 'False', 'unresolved coverage'), ('wl', 'complement', 'True', 'sampled')]
161:SYNTHETIC sample0 {'truth': False} {'boolean_difference': False}
```

## O1 (orchestrator): attached domains stated through a set or a quantifier are not compared
```
$ grep -n "raise Unsupported('unsupported constructor '+h)" research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
543:    raise Unsupported('unsupported constructor '+h)
```

```
$ grep -c "'ConditionSet'\|'ForAll'\|'Exists'" research/pde_ledger_v3/scripts/S9b_cross_engine_comparator.py
0
```

```
$ wc -l < _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/out/census_outcomes.stdout
348
```

```
$ grep -o "'nf:[^']*'" _scratch/s9b_build/s9b_comparator_build_review_r0_claude_evidence/out/census_outcomes.stdout | sort | uniq -c | sort -rn
    328 'nf:nested sibling has no bound counterpart'
    270 'nf:unsupported constructor ConditionSet'
     12 'nf:container/scalar representation differs'
      1 'nf:unsupported applied head Inactive/ForAll'
```

```
$ grep -A9 '"comparison_outcomes"' _scratch/s9b_build/comparator_output/full-04/output_summary.json
  "comparison_outcomes": {
    "container": 942,
    "not_formed": 1812,
    "wrapper": 408,
    "symbolic": 252,
    "marker": 265,
    "predicate": 45,
    "label": 44
  },
  "samples": {
```

