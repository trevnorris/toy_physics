# Measurements — O2 comparator build r4 review dispositions (generated 2026-10-07 20:29)

Generator: `_scratch/s9b_build/gen/o2_comparator_build_r4_lookups.sh` (sed/grep/sha256sum only; lines cut at 300 characters). Leg evidence quoted verbatim from the Claude leg's files, copied from `/tmp/o2r4_fc/` to `_scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/`.

```
$ sha256sum -c _scratch/s9b_build/o2_comparator_build_review_baseline_r4.sha256
research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_O2_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/ablate_O2_cross_engine_comparator.py: OK
_scratch/s11c/o2-comparator-build-r4/comparison.jsonl: OK
_scratch/s11c/o2-comparator-build-r4/accounting.jsonl: OK
_scratch/s11c/o2-comparator-build-r4/catalog.json: OK
_scratch/s11c/o2-comparator-build-r4/validation.json: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r4_report.md: OK
research/pde_ledger_v3/_measurements/O2_comparator_builder_r4_measurements.txt: OK
```

## W1 — the named-OPEN inventory on the structural-role path; WL provenance and held-operator names in the production output
```
$ sed -n '1054,1069p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
        if v.head=='Apply':
            key=action_key(v,engine)
            start=2 if engine=='wl' and v.args[0]==atom('Symbol','OpenAction') else 1
            if key is not None:
                own_head=sum(algebra_leaves(x) for x in v.args[:start]) if root else 0
                if root and key.startswith('native:'):
                    # OpenNativeField[OPEN[name, provenance]]: only these
                    # three name leaves identify the role, not its provenance.
                    native_head=v.args[0]
                    descriptor=native_head.args[1]
                    own_head=sum(algebra_leaves(x) for x in
                                 (native_head.args[0],descriptor.args[0],descriptor.args[1]))
                return own_head+sum(visit(x,scope) for x in v.args[start:])
            if engine=='wl' and v.args[0]==atom('Symbol','OPEN'):
                return visit(v.args[1],scope) if len(v.args)>1 else 0
            return sum(visit(x,scope,head=i==0) for i,x in enumerate(v.args))
```

```
$ sed -n '995p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
    'leaf_policy':'parsed terminal nodes and constructor options; consumed leaves construct a role/head, named operand/label, or complete live-object key; options, binder declarations outside a live key, provenance and other syntax are outside; no leaf is counted twice',
```

```
$ for t in '"wl::6"' '"wl::3,8"' '"wl::D"' '"wl::Map"'; do printf '%s lines=' "$t"; grep -cF "$t" _scratch/s11c/o2-comparator-build-r4/comparison.jsonl; done
"wl::6" lines=6
"wl::3,8" lines=6
"wl::D" lines=20
"wl::Map" lines=20
```

```
$ sed -n '266p' research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl
    "GeneralizedRates" -> operand[UnspecifiedRotationalNormalRates,"6","3,8"],
```

```
$ cat _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/generalized_rates_quote.stdout | head -12
operands line 373 py_path ['MATERIAL_POWER_PAIRING', 'unfixed_normal_generalized_rates'] wl_path ['B_E_STEADY', 'Entries', '2', 'Object', '@3', 'MaterialWork', '@3', 'GeneralizedRates']
  WL operand: ["Apply", "", [], [["Symbol", "OPEN", [], []], ["Symbol", "UnspecifiedRotationalNormalRates", [], []], ["Record", "", [], [["Entry", "Spec", [], [["Text", "6", [], []]]], ["Entry", "Contract", [], [["Text", "3,8", [], []]]]]]]]
role OPEN_UnfixedNormalGeneralizedRates
  wl named {'wl::UnspecifiedRotationalNormalRates': 1, 'wl::6': 1, 'wl::3,8': 1}
  wl limit {'repeated_live_objects': {}, 'parsed_leaves': 4, 'inventory_consumed_leaves': 4, 'outside_inventory_leaves': 0}
  differences.named [['B_A13', 1], ['H_core', 1], ['J_map', 1], ['N_br_live', 1], ['S12_boundary_domain', 1], ['S12_local_source_controller', 1], ['py::bulk_state', 1], ['py::mouth_core_data', 1], ['wl::3,8', -1], ['wl::6', -1], ['wl::UnspecifiedRotationalNormalRates', -1]]
role py::OPEN_OpenDomain
  wl named None
  wl limit None
  differences.named []
role py::OPEN_StateHistory
  wl named None
```

```
$ cat _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/head_leaf_quote.stdout | head -8
internal_force row, PY occurrence: parsed 3378 consumed 1616 outside 1762
hold_inplane[0] balance entry, PY: orientation -1 parsed 3379 consumed 1615 outside 1764
```

## W2 — balance-entry orientation comes from term_orientation; the harness orientation knives edit signed_actions
```
$ sed -n '1142,1161p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
def term_orientation(n):
    sign=1
    if n.head=='Number' and n.value.startswith('-'):
        return -1
    if n.head=='Mul':
        for v in n.args:
            sign *= term_orientation(v)
    elif n.head=='Derivative':
        sign *= term_orientation(n.args[0])
    elif n.head=='Lambda':
        sign *= term_orientation(n.args[1])
    elif n.head=='Apply':
        head=n.args[0]
        if head.value in ('OpenFirstVariation','Function','OPEN_SumOverAllNativeFaces'):
            body=1 if head.value=='OpenFirstVariation' else 2
            sign *= term_orientation(n.args[body])
        elif (head.head=='Apply' and head.args[0].value=='Inactive' and
              head.args[1].value in ('Total','Map','D')):
            sign *= term_orientation(n.args[1])
    return sign
```

```
$ sed -n '1185,1190p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
            fields,limit=inventory_details(term,bindings,engine)
            fields['head']=Counter({role:1})
            fields['orientation']=Counter({str(term_orientation(term)):1})
            fields['role']=Counter({role:1})
            entries.append({'role':role,'orientation':term_orientation(term),
                            'open_free':not is_open,'operand':data(mapped),'inventories':fields,'limit':limit})
```

```
$ grep -n 'signed_actions\|term_orientation' research/pde_ledger_v3/scripts/ablate_O2_cross_engine_comparator.py
85:      ('signed_actions(child,sign if i == body else 1,dependency or i != body)',
86:       'signed_actions(child,1,dependency or i != body)')),
```

```
$ cat _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/form_tests_summary.stdout
== F1_term_lambda_body guard_exit=0 OK
== F2_term_inactive_aggregate guard_exit=0 OK
== F3_term_held_body guard_exit=0 OK
== F4_term_derivative_body guard_exit=0 FAILED (failures=1)
FAIL: test_held_balance_entry_orientation (__main__.Controls)
== F5_term_negative_number guard_exit=0 FAILED (failures=3)
FAIL: test_balance_entry_field_differences (__main__.Controls) (field='orientation')
FAIL: test_closed_balance_entry (__main__.Controls)
FAIL: test_held_balance_entry_orientation (__main__.Controls)
== F6_symbol_pair_subtracted guard_exit=0 FAILED (failures=1)
FAIL: test_matching_text_and_names_not_zero_evidence (__main__.Controls)
```

```
$ grep -h 'differing lines' _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/diff_F1.stdout _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/diff_F2.stdout _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/diff_F3.stdout
differing lines 0
differing lines 0
differing lines 0
```

## W3 — the printed not_compared list; one-sided corruption D (held-derivative coordinate)
```
$ sed -n '988,996p' research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
OPEN_SCOPE = {
    'compared_inventories':['role_or_head_and_orientation','named_OPEN_operands_including_labels',
                           'live_profiles_and_derivatives_at_arguments'],
    'empty_difference_means':'only these inventories match',
    'not_compared':['where an object sits among OPEN arguments','how many times an object occurs',
                    'argument content outside these inventories, including profile-free velocity or metric algebra'],
    'occurrence_indices':'local engine enumeration only; no cross-engine argument or occurrence slot pairing',
    'leaf_policy':'parsed terminal nodes and constructor options; consumed leaves construct a role/head, named operand/label, or complete live-object key; options, binder declarations outside a live key, provenance and other syntax are outside; no leaf is counted twice',
    'raw_syntax_policy':'raw trees, censuses and serialization hashes are diagnostics, not additional OPEN comparisons',
```

```
$ tail -6 _scratch/s9b_build/o2_comparator_build_review_r4_claude_evidence/diff_D.stdout
   (operand text differs; printed operand reflects the input)
LINE 642 kind balance_entries row hold_inplane component [0]
   at ['entries', 'py', 5, 'operand', 3, 1, 3, 0, 1]
      base : x1
      other: x2
differing lines 5
```

## Note N1 — generalized-rates grain (accounting)
```
$ grep -n "generalized_rates\|UnspecifiedRotationalNormalRates'\|MATERIAL_POWER_PAIRING/rotational" research/pde_ledger_v3/scripts/O2_cross_engine_comparator.py
452:    J('generalized_rates','MATERIAL_POWER_PAIRING/unfixed_normal_generalized_rates','B_E_STEADY/Entries/2/Object/@3/MaterialWork/@3/GeneralizedRates',277,266,'§6'),
557:    ('OPEN_UnfixedNormalGeneralizedRates','operand:UnspecifiedRotationalNormalRates',277,266,'§6 generalized rates')))
763:    if name(head)=='OPEN' and len(n.args)>1 and name(n.args[1])=='UnspecifiedRotationalNormalRates':
764:        return 'operand:UnspecifiedRotationalNormalRates'
1317:        'MATERIAL_POWER_PAIRING/rotational':'WL retains rotational work within its joint material-work action, with no separately attributable rotational-work action (§6).',
```

```
$ sed -n '277,280p' research/pde_ledger_v3/scripts/O2_live_balance_sympy_audit.py
    normal_generalized_rates = action('UnfixedNormalGeneralizedRates',
                                      op['N_br_live'], op['J_map'], state, geom)
    rotational_power = action('RotationalGeneralizedWork', op['A_rot_live'],
                              material, stress_state, state)
```

```
$ sed -n '263,266p' research/pde_ledger_v3/mathematica/O2_live_balance_mathematica_audit.wl
stressWork = response[MaterialStressNormalRotationalWork,
  {stressInputs,inertia,normalResponse,rotation},
  <|"ForceAction" -> internalForce,"GraphVelocity" -> vMaterial,
    "GeneralizedRates" -> operand[UnspecifiedRotationalNormalRates,"6","3,8"],
```

```
$ grep -o '"path":\["MATERIAL_POWER_PAIRING","rotational"\],"reason":"[^"]*"' _scratch/s11c/o2-comparator-build-r4/accounting.jsonl
"path":["MATERIAL_POWER_PAIRING","rotational"],"reason":"WL retains rotational work within its joint material-work action, with no separately attributable rotational-work action (§6)."
```

