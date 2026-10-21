# Measurements: S14a part 1 comparator build review round 1 lookups (generated 2026-10-21 14:44)

Generator: `_scratch/s14a/gen/cmp_build_r1_lookups.sh` (grep/sed/sha256sum/diff/git/cat/sort/uniq only; this
file is written by it). Mechanical lookups only; every CAS result quoted here is a leg's saved stdout.

## E1. The round: repair, baseline, prompt, tree, reports

```
$ cat /var/projects/toy_physics/_scratch/s14a/cmp/repair1.DONE; grep -A1 'tokens used' /var/projects/toy_physics/_scratch/s14a/cmp/repair1.log | tail -2
exit 0
tokens used
615,948
```

```
$ cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s14a/cmp/review/s14a1_cmp_build_baseline_r1.sha256
research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_run_record.md: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_builder_report.md: OK
_scratch/s14a/cmp_work/r1_full_01.jsonl: OK
_scratch/s14a/cmp_work/r1_full_01_inventory.json: OK
_scratch/s14a/cmp_work/r1_full_01_unpaired_components.json: OK
_scratch/s14a/cmp_work/r1_tests_01.out: OK
_scratch/s14a/cmp_work/r1_run_inventory.json: OK
```

```
$ cd /var/projects/toy_physics && git status --porcelain | diff - /var/projects/toy_physics/_scratch/s14a/cmp/review/status_before_cmp_r1.txt && echo 'status unchanged'
6,7d5
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_build_review_r1_disposition.md
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_build_review_r1_lookups.md
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s14a/cmp/review/s14a1_cmp_build_review_prompt_r1.md
75eecf5f6987bffd4f6b2b74798e62f7a40ecb685565bff9e1792ea7affae9bf  /var/projects/toy_physics/_scratch/s14a/cmp/review/s14a1_cmp_build_review_prompt_r1.md
```

```
$ for r in r1_tests_01 r1_full_01; do printf '%s ' $r; grep -oE '"exitCode": *[0-9]+' /var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/$r/outcome.json | tr '\n' ' '; echo; done
r1_tests_01 "exitCode": 0 "exitCode": 0 
r1_full_01 "exitCode": 0 "exitCode": 0 
```

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r1-grok/DONE; grep -nE 'Verdict: ' /var/projects/_review_scratch/s14a1-cmp-r1-claude/REVIEW.md /var/projects/_review_scratch/s14a1-cmp-r1-grok/REVIEW.md
exit 0
/var/projects/_review_scratch/s14a1-cmp-r1-claude/REVIEW.md:3:**Verdict: CLEAR.** No findings. I found nothing that changes what the comparator computes or what may be claimed from
/var/projects/_review_scratch/s14a1-cmp-r1-grok/REVIEW.md:3:**Verdict: CLEAR**
```

```
$ grep -nE '^## Findings' -A2 /var/projects/_review_scratch/s14a1-cmp-r1-claude/REVIEW.md; grep -nE '^\*\*Findings:\*\*' /var/projects/_review_scratch/s14a1-cmp-r1-grok/REVIEW.md
12:## Findings
13-
14-None.
5:**Findings:** none.
```

```
$ cd /var/projects/toy_physics && git diff --stat 88bf7bf9 -- /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
 .../scripts/S14a_1_cross_engine_comparator.py      | 110 ++++++++++++---
 .../scripts/test_S14a_1_cross_engine_comparator.py | 149 +++++++++++++++++++++
 2 files changed, 242 insertions(+), 17 deletions(-)
```

## E2. H1: integral binders by nesting depth

```
$ sed -n '719,735p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
def alpha_integrals(expr):
    """Name integral binders by active lexical depth, without capture.

    Conversion first uses fresh Dummies. Only the completed expression is
    renamed: an inner integral cannot give an outer variable its own name.
    Flattened multi-limit integrals and nested integral syntax share the
    same outside-in binder order. Sibling scopes may reuse names safely.
    """
    def walk(value,env,depth):
        if isinstance(value,s.Integral):
            local=dict(env)
            for i,var in enumerate(reversed(value.variables)):
                local[var]=s.Symbol('_int'+str(depth+i),real=True)
            inner_depth=depth+len(value.limits)
            limits=[(local[var],walk(lo,local,inner_depth),walk(hi,local,inner_depth))
                    for var,lo,hi in value.limits]
            return s.Integral(walk(value.function,local,inner_depth),*limits)
```

```
$ sed -n '306,308p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
    def test_nested_integral_variable_identity(self):
        py='Integral(radius**2*Integral(radius*w,(w,0,1)),(radius,1,2))'
        wl='Inactive[Integrate][r^2*Inactive[Integrate][r*w,{w,0,1}],{r,1,2}]'
```

```
$ head -4 /var/projects/_review_scratch/s14a1-cmp-r1-claude/corrupt/C1_wl_nested_slot_order_lost/pdiff.stdout | cut -c1-200
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 1
  b169 [SphericalSpecialization/order lost [] | spherical specialization / shell/order lost []] changed items: 56
     B: as_computed_residual "-4*pi*Integral(_int0**2*Integral(SphericalSpecialization_inside_GammaB(_int0, _int1, t)*SphericalSpecialization_inside_n(_int0, _int1, t), (_int1, SphericalSpecialization
```

```
$ head -3 /var/projects/_review_scratch/s14a1-cmp-r1-claude/corrupt/C2_py_nested_slot_storage/pdiff.stdout | cut -c1-200
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 1
  b171 [SphericalSpecialization/storage [] | spherical specialization / shell/storage []] changed items: 38
```

```
$ grep -nE 'bindings changed in operands/values only' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/U1_undo_H1/pdiff.stdout | cut -c1-200
3:bindings changed in operands/values only: 5 ['b169', 'b171', 'b173', 'b174', 'b175']
```

```
$ grep -nE '^FAIL:|^FAILED' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/U1_undo_H1/tests.out | cut -c1-200
35:FAIL: test_nested_integral_variable_identity (__main__.Controls.test_nested_integral_variable_identity)
50:FAILED (failures=1)
```

```
$ grep -nE '_int1|_int0' /var/projects/_review_scratch/s14a1-cmp-r1-grok/ablations/extract.txt
1:baseline b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
2:undo_h1 b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
3:undo_h2 b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
4:undo_h3 b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
5:face b171 residual_len 429 has_int1 False has_inner_radius_deriv True nonzero_samples 18
7:face_h2 b171 residual_len 429 has_int1 False has_inner_radius_deriv True nonzero_samples 0
8:nested b171 residual_len 1005 has_int1 True has_inner_radius_deriv False nonzero_samples 36
10:raw b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
11:rerun b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
12:drop_condition b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
13:endpoint b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
22:baseline py_int1 2 wl_int1 2 py_int0 11
23:undo_h1 py_int1 0 wl_int1 0 py_int0 13
```

```
$ grep -nE '^nested |sample_equal baseline nested' -A1 /var/projects/_review_scratch/s14a1-cmp-r1-grok/ablations/extract.txt | cut -c1-120
8:nested b171 residual_len 1005 has_int1 True has_inner_radius_deriv False nonzero_samples 36
9-  first_nonzero [('affine', 'jump_lower', 'inside', 0, '-746458674039282257005226638767*pi/182308782128284134400000000
--
17:sample_equal baseline nested False
18-  differing_keys 72 by_form {'affine': 18, 'quadratic': 36, 'cubic': 18}
```

## E3. H2: a derivative taken at another point

```
$ sed -n '870,884p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
def argument_coupling(args,code):
    """Every first partial depends on every other argument in this form."""
    return (code+1)*sum((i+1)*(j+1)*args[i]*args[j]
                        for i in range(len(args)) for j in range(i+1,len(args)))/211


def basic_form(form,args,code):
    if form=='affine':
        return (code+1)/7+sum((i+1)*(code+1)*q/101 for i,q in enumerate(args))
    if form=='quadratic':
        return (code+1)/7+sum(((i+1)**2*(code+1)*q*q+(i+1)*q)/103 for i,q in enumerate(args))+argument_coupling(args,code)
    if form=='cubic':
        return (code+1)/7+sum(((i+1)**3*(code+1)*q**3+(i+1)**2*q)/107 for i,q in enumerate(args))
    raise InputError('unknown probe form '+form)

```

```
$ sed -n '333,338p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py
    def test_every_live_first_partial_at_another_argument_point(self):
        # H2 requires one discriminating form/stratum. The quadratic probes
        # cover every derivative slot crossed with every evaluation argument;
        # production still samples all forms, strata and rational points.
        options={'forms':('quadratic',),'strata':('generic',)}
        for field in (r for r in c.NAME_TABLE if r.arity):
```

```
$ for d in C3_wl_face_dt_point_storage C4_wl_field_dw_point_closed_operand_S C5_wl_field_dw_point_closed_operand_G C10_wl_face_dr_point_J; do echo "== $d"; sed -n '2,3p' /var/projects/_review_scratch/s14a1-cmp-r1-claude/corrupt/$d/pdiff.stdout | cut -c1-110; done
== C3_wl_face_dt_point_storage
bindings whose residual products moved: 1
  b171 [SphericalSpecialization/storage [] | spherical specialization / shell/storage []] changed items: 20
== C4_wl_field_dw_point_closed_operand_S
bindings whose residual products moved: 1
  b128 [SphericalSpecialization/material/closed integral of differential operand [] | spherical specialization
== C5_wl_field_dw_point_closed_operand_G
bindings whose residual products moved: 1
  b50 [General/material/closed integral of differential operand [] | general projection/material/closed integr
== C10_wl_face_dr_point_J
bindings whose residual products moved: 1
  b168 [SphericalSpecialization/J [] | spherical specialization / shell/J []] changed items: 20
```

```
$ grep -nE 'B: sample' /var/projects/_review_scratch/s14a1-cmp-r1-claude/corrupt/C3_wl_face_dt_point_storage/pdiff.stdout | head -1 | cut -c1-200
6:     B: sample {"form": "quadratic", "location": "inside", "point_index": 0, "residual": "-77477917886174994067*pi/473208382034853830400", "stratum": "generic"}
```

```
$ grep -nE '^FAILED' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/U2_undo_H2/tests.out; grep -m1 -nE '^FAIL: test_every_live_first_partial' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/U2_undo_H2/tests.out | cut -c1-200
10147:FAILED (failures=919)
953:FAIL: test_every_live_first_partial_at_another_argument_point (__main__.Controls.test_every_live_first_partial_at_another_argument_point) (pair='slot 0 at argument 1', function='General_face_lower
```

```
$ grep -nE '^(face|face_h2|baseline) b171|sample_equal (baseline face|undo_h2 face_h2)' -A1 /var/projects/_review_scratch/s14a1-cmp-r1-grok/ablations/extract.txt | cut -c1-140
1:baseline b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
2-undo_h1 b171 residual_len 422 has_int1 False has_inner_radius_deriv False nonzero_samples 0
--
5:face b171 residual_len 429 has_int1 False has_inner_radius_deriv True nonzero_samples 18
6-  first_nonzero [('quadratic', 'generic', 'inside', 0, '-77477917886174994067*pi/473208382034853830400', None), ('quadratic', 'generic', '
7:face_h2 b171 residual_len 429 has_int1 False has_inner_radius_deriv True nonzero_samples 0
8-nested b171 residual_len 1005 has_int1 True has_inner_radius_deriv False nonzero_samples 36
--
14:sample_equal baseline face False
15-  differing_keys 36 by_form {'quadratic': 36}
16:sample_equal undo_h2 face_h2 True
17-sample_equal baseline nested False
```

## E4. H3: unjoined reasons

```
$ sed -n '499,500p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
def is_label(tag):
    return tag in ('contract','engine','WORKING_ON') or is_comparison_note(tag) or tag.endswith('/representation')
```

```
$ grep -F '"tag":"transverse integration / comparison"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl | cut -c1-400
{"kind":"tag_accounting","engine":"wl","tag":"transverse integration / comparison","line":9,"occurrence":0,"parsed_leaves":1,"actually_compared_leaves":0,"accounting":"unjoined","reason":"Engine-local comparison note; this occurrence has no unique join under the active table.","classification":"engine-local comparison note","read_as_label":true,"reason_basis":"reader label classification and curre
```

```
$ grep '"kind":"tag_accounting"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl | grep -oE '"accounting":"[a-z_]+"' | sort | uniq -c
    327 "accounting":"bound"
    656 "accounting":"unjoined"
```

```
$ grep '"kind":"tag_accounting"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl | grep '"accounting":"unjoined"' | sed -E 's/.*"accounting":"unjoined","reason":"([^"(:;]{0,40}).*/\1/' | sort | uniq -c | sort -rn
    343 Reference-table/source curation 
    313 Engine-local comparison note
```

```
$ grep -c 'Reference-table/source curation (not computed active-table or counterpart-completeness evidence)' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl
343
```

```
$ grep -nE '^FAILED' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/U3_undo_H3/tests.out; grep -m1 -nE '^FAIL: test_freed_component' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/U3_undo_H3/tests.out | cut -c1-200
1422:FAILED (failures=126)
160:FAIL: test_freed_component_bindings_have_current_accounting (__main__.Controls.test_freed_component_bindings_have_current_accounting) (binding='b3')
```

```
$ grep -nE 'wl9' /var/projects/_review_scratch/s14a1-cmp-r1-grok/ablations/extract.txt | cut -c1-260
24:baseline wl9 {"tag": "transverse integration / comparison", "accounting": "unjoined", "reason": "Engine-local comparison note; this occurrence has no unique join under the active table.", "classification": "engine-local comparison note", "basis": "reader la
25:undo_h3 wl9 {"tag": "transverse integration / comparison", "accounting": "unjoined", "reason": "Wolfram-only fixed-slab coordinate construction; SymPy prints native-w integration objects (py171-185; wl64-99,108-125,143-155). Different printed objects; no co
```

## E5. What the repair changed on the real streams (round 0 output against round 1 output)

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r1-claude/misc/r0_vs_r1_accounting_status.stdout | cut -c1-200
tags r0: 983 tags r1: 983 keys only in one: 0
tags whose (accounting, compared leaves, parsed leaves) changed: 0
```

```
$ grep -nE '^bindings|^tag accounting|^summary|^below' /var/projects/_review_scratch/s14a1-cmp-r1-claude/misc/r0_full02_vs_r1_full01.pdiff.stdout | cut -c1-200
1:bindings in A: 227 in B: 227 | present in only one: 0 []
2:bindings whose residual products moved: 0
3:bindings changed in operands/values only: 190 ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15', 'b16', 'b17', 'b18', 'b19', 'b20', 'b21', 'b22', 'b23',
4:tag accounting changed (accounting/reason/compared leaves): 656
13:summary count changes: {}
14:below parsed count A/B: [] / []
```

## E6. FORM ablations (Claude leg), production effect and test outcome

```
$ for d in A1_raw_delta A2a_one_form A2b_nonseparating A3_same_test_function A4_no_range A5_cas_endpoint_weight A6_drop_condition A7_position_pairing A8_whole_tags_only A8b_whole_tag_refs A9_collapse_applied; do echo "== $d"; grep -E '^bindings whose residual|^bindings changed in operands|^bindings in A' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/$d/pdiff.stdout | cut -c1-140; grep -E '^(OK|FAILED)' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/$d/tests.out; done
== A1_raw_delta
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 0
bindings changed in operands/values only: 4 ['b48', 'b70', 'b126', 'b146']
FAILED (failures=1)
== A2a_one_form
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 190
bindings changed in operands/values only: 0 []
FAILED (failures=2)
== A2b_nonseparating
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 0
bindings changed in operands/values only: 190 ['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15
FAILED (failures=2)
== A3_same_test_function
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 0
bindings changed in operands/values only: 35 ['b105', 'b107', 'b110', 'b112', 'b126', 'b127', 'b128', 'b130', 'b131', 'b133', 'b134', 'b136'
FAILED (failures=1)
== A4_no_range
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 11
bindings changed in operands/values only: 30 ['b51', 'b52', 'b59', 'b60', 'b61', 'b73', 'b74', 'b79', 'b80', 'b81', 'b82', 'b129', 'b130', '
FAILED (failures=3)
== A5_cas_endpoint_weight
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 0
bindings changed in operands/values only: 0 []
FAILED (failures=1)
== A6_drop_condition
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 5
bindings changed in operands/values only: 0 []
FAILED (failures=4)
== A7_position_pairing
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 3
bindings changed in operands/values only: 0 []
FAILED (failures=4, errors=4)
== A8_whole_tags_only
bindings in A: 227 in B: 102 | present in only one: 125 ['b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15',
bindings whose residual products moved: 0
bindings changed in operands/values only: 0 []
FAILED (errors=1)
== A8b_whole_tag_refs
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 158
bindings changed in operands/values only: 0 []
FAILED (failures=1505, errors=1)
== A9_collapse_applied
bindings in A: 227 in B: 227 | present in only one: 0 []
bindings whose residual products moved: 141
bindings changed in operands/values only: 0 []
FAILED (failures=70, errors=1)
```

```
$ grep -nE '^(FAIL|ERROR):' /var/projects/_review_scratch/s14a1-cmp-r1-claude/abl/A4_no_range/tests.out | cut -c1-160
35:FAIL: test_condition_only_mutation_and_failure_stratum (__main__.Controls.test_condition_only_mutation_and_failure_stratum)
45:FAIL: test_equivalence_integral_alpha_and_closed_endpoint (__main__.Controls.test_equivalence_integral_alpha_and_closed_endpoint)
63:FAIL: test_integral_range_variable_integrand (__main__.Controls.test_integral_range_variable_integrand)
```

## E7. FORM ablations (Grok leg), first test failure and moved bindings

```
$ grep -nE '^CASE |^tests_exit|^first_failure|^flags' /var/projects/_review_scratch/s14a1-cmp-r1-grok/ablations/report.txt | cut -c1-200
1:CASE operational_existing_output
20:CASE undo_h1
21:tests_exit 1 tests_seconds 52.384
22:first_failure FAIL: test_nested_integral_variable_identity (__main__.Controls.test_nested_integral_variable_identity)
41:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 5}
51:CASE undo_h2
52:tests_exit 1 tests_seconds 32.439
53:first_failure FAIL: test_every_live_first_partial_at_another_argument_point (__main__.Controls.test_every_live_first_partial_at_another_argument_point) (pair='slot 0 at argument 1', function='Gener
69:flags {"tags": false, "test_functions": true, "forms": false, "summary_counts": false, "changed_bindings": 0}
72:CASE undo_h3
73:tests_exit 1 tests_seconds 50.713
74:first_failure FAIL: test_freed_component_bindings_have_current_accounting (__main__.Controls.test_freed_component_bindings_have_current_accounting) (binding='b3')
90:flags {"tags": true, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 0}
93:CASE raw_coefficient
94:tests_exit 1 tests_seconds 5.683
95:first_failure FAIL: test_delta_face_removal_normalization_and_regular_region (__main__.Controls.test_delta_face_removal_normalization_and_regular_region)
113:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 4}
123:CASE one_form
124:tests_exit 1 tests_seconds 1.692
125:first_failure FAIL: test_coverage_one_form_only (__main__.Controls.test_coverage_one_form_only)
141:flags {"tags": false, "test_functions": true, "forms": true, "summary_counts": false, "changed_bindings": 190}
151:CASE nonseparating_forms
152:tests_exit 1 tests_seconds 4.504
153:first_failure FAIL: test_coverage_nonseparating_forms_only (__main__.Controls.test_coverage_nonseparating_forms_only)
169:flags {"tags": false, "test_functions": true, "forms": false, "summary_counts": false, "changed_bindings": 0}
172:CASE same_test_function
173:tests_exit 1 tests_seconds 6.271
174:first_failure FAIL: test_distinct_functions_in_each_form_and_stratum (__main__.Controls.test_distinct_functions_in_each_form_and_stratum)
190:flags {"tags": false, "test_functions": true, "forms": false, "summary_counts": false, "changed_bindings": 4}
200:CASE drop_integral_range
201:tests_exit 1 tests_seconds 6.637
202:first_failure FAIL: test_condition_only_mutation_and_failure_stratum (__main__.Controls.test_condition_only_mutation_and_failure_stratum)
218:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 41}
228:CASE endpoint_half_weight
229:tests_exit 1 tests_seconds 33.614
230:first_failure FAIL: test_equivalence_integral_alpha_and_closed_endpoint (__main__.Controls.test_equivalence_integral_alpha_and_closed_endpoint)
256:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 0}
259:CASE drop_condition
260:tests_exit 1 tests_seconds 4.816
261:first_failure FAIL: test_condition_only_mutation_and_failure_stratum (__main__.Controls.test_condition_only_mutation_and_failure_stratum)
275:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 2}
285:CASE pair_by_position
286:tests_exit 1 tests_seconds 14.844
287:first_failure FAIL: test_each_production_component_correspondence (__main__.Controls.test_each_production_component_correspondence) (binding='b2', component='below')
303:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": true, "changed_bindings": 11}
313:CASE whole_tags
314:tests_exit 1 tests_seconds 0.330
315:first_failure FAIL: test_condition_only_mutation_and_failure_stratum (__main__.Controls.test_condition_only_mutation_and_failure_stratum)
334:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": true, "changed_bindings": 161}
344:CASE collapse_applied
345:tests_exit 1 tests_seconds 0.567
346:first_failure FAIL: test_condition_only_mutation_and_failure_stratum (__main__.Controls.test_condition_only_mutation_and_failure_stratum)
362:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": true, "changed_bindings": 141}
372:CASE corrupt_nested_slot
373:tests_exit 0 tests_seconds 77.329
374:first_failure no failure line
384:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 1}
394:CASE corrupt_face_point
395:tests_exit 0 tests_seconds 76.953
396:first_failure no failure line
406:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 1}
416:CASE corrupt_face_point_undo_h2
417:tests_exit 1 tests_seconds 32.647
418:first_failure FAIL: test_every_live_first_partial_at_another_argument_point (__main__.Controls.test_every_live_first_partial_at_another_argument_point) (pair='slot 0 at argument 1', function='Gene
434:flags {"tags": false, "test_functions": true, "forms": false, "summary_counts": false, "changed_bindings": 1}
444:CASE rerun_final
445:tests_exit 0 tests_seconds 77.644
446:first_failure no failure line
456:flags {"tags": false, "test_functions": false, "forms": false, "summary_counts": false, "changed_bindings": 0}
```

## E8. One-sided corruptions (Claude leg C6–C9)

```
$ for d in C6_py_form_core_channel C7_wl_delta_coeff_off_face_equiv C8_wl_delta_part_upper C9_wl_condition_only_volume_flux; do echo "== $d"; sed -n '2,3p' /var/projects/_review_scratch/s14a1-cmp-r1-claude/corrupt/$d/pdiff.stdout | cut -c1-120; grep -E '^bindings changed in operands' /var/projects/_review_scratch/s14a1-cmp-r1-claude/corrupt/$d/pdiff.stdout | cut -c1-120; done
== C6_py_form_core_channel
bindings whose residual products moved: 1
  b167 [SphericalSpecialization/core channel [] | spherical specialization / shell/core channel []] changed items: 56
bindings changed in operands/values only: 0 []
== C7_wl_delta_coeff_off_face_equiv
bindings whose residual products moved: 0
bindings changed in operands/values only: 0 []
bindings changed in operands/values only: 0 []
== C8_wl_delta_part_upper
bindings whose residual products moved: 1
  b126 [SphericalSpecialization/material/differentiated distribution [] | spherical specialization/material/distribution
bindings changed in operands/values only: 0 []
== C9_wl_condition_only_volume_flux
bindings whose residual products moved: 1
  b177 [SphericalSpecialization/outer volume flux [] | spherical specialization / shell/projected velocity volume flux [
bindings changed in operands/values only: 0 []
```

## E9. Parts whose exact residual as computed is not literally 0 or False on the real streams

```
$ grep '"kind":"residual_before_domain"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl | grep -vE '"as_computed_residual":"(0|False)"' | grep -oE '"binding":"b[0-9]+","part":"[^"]*"' | sort -u
"binding":"b128","part":"item 1: closed integral of material differential operand"
"binding":"b148","part":"item 1: closed integral of order differential operand"
"binding":"b165","part":"items 2-3: shell domain"
"binding":"b171","part":"items 2-3: storage"
"binding":"b173","part":"shell left"
"binding":"b174","part":"item 2: shell balance residual"
"binding":"b188","part":"external printed domain contexts"
"binding":"b189","part":"external printed domain contexts"
"binding":"b50","part":"item 1: closed integral of material differential operand"
"binding":"b72","part":"item 1: closed integral of order differential operand"
```

```
$ grep -c '"as_computed_residual"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl
521
```

```
$ grep '"kind":"operands"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl | grep -E '"binding":"b(188|189)"' | grep -oE '"binding":"b[0-9]+","py":\{"tag":"[^"]*"'
"binding":"b188","py":{"tag":"SphericalSpecialization/dimensions/projected velocity/dimension"
"binding":"b189","py":{"tag":"SphericalSpecialization/dimensions/volume flux/dimension"
```

```
$ grep -nE 'closed-endpoint' /var/projects/_review_scratch/s14a1-cmp-r1-grok/REVIEW.md | cut -c1-160
34:| Endpoint weight one half | Comparison fingerprint unchanged. These streams contain 0 closed-endpoint integrals. | `test_equivalence_integral_alpha_and_clos
```

```
$ grep '"kind":"residual_before_domain"' /var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl | grep -E '"binding":"b(188|189)","part":"external' | grep -oE '"binding":"b[0-9]+"|"(py_operand|wl_operand)":"[^"]{0,90}'
"binding":"b188"
"py_operand":"True
"wl_operand":"Ne(0, Integral(SphericalSpecialization_inside_n(outer_radius, _int0, t), (_int0, Spherical
"binding":"b189"
"py_operand":"True
"wl_operand":"Ne(0, Integral(SphericalSpecialization_inside_n(outer_radius, _int0, t), (_int0, Spherical
```

## E10. Reading, tokens, reproducibility

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r1-claude/misc/read_coverage.stdout | cut -c1-200 | head -20
Counter({'tag_accounting': 983, 'probe_domain': 540, 'residual_before_domain': 331, 'operands': 227, 'comparison': 227, 'binding_accounting': 227, 'progress_marker': 27, 'test_functions': 18, 'configu
{"kind": "summary", "counts": {"bound_parsed_leaves_py": 21714, "compared_leaves_py": 21714, "bound_parsed_leaves_wl": 24948, "compared_leaves_wl": 24948, "bindings": 227, "bindings_not_formed": 0, "b
py file lines 619 tag start lines 600 progress markers 13 remaining (continuation) lines 6
wl file lines 397 tag start lines 383 progress markers 14 remaining (continuation) lines 0
```

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r1-claude/misc/item6_tokens.stdout | cut -c1-200 | head -20
$ grep -c -E "\b(PASS|FAIL|AGREE|MATCH|VERDICT|STATUS)\b" r1_full_01.jsonl
0
$ grep -n -E "\b(PASS|FAIL|AGREE|MATCH|VERDICT|STATUS)\b" comparator.py
$ grep -n -i -E "expected|target" comparator.py
2:"""S14a part 1 native-stream comparator. No residual target was given.
77:    def take(self, expected=None):
79:        if expected is not None and v != expected:
80:            raise InputError(f'Wolfram expected {expected!r}, found {v!r}')
121:            raise InputError('unexpected Wolfram token ' + repr(tok))
590:                raise Unsupported('expected sequence at '+str(ref.path))
595:                raise Unsupported('expected conditioned object at '+str(ref.path))
600:                raise Unsupported('expected equation at '+str(ref.path))
1439:              'residual_target':'No residual target was given.'})
$ grep -n -E "sys.exit|exit\(|raise SystemExit" comparator.py
1603:    raise SystemExit(main())
```

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r1-claude/misc/repro.stdout | cut -c1-200
/var/projects/toy_physics/_scratch/s14a/cmp_work/r1_full_01.jsonl records(non-resources), resources, sha256(non-resources) = (2582, 1, '640b95d10248b2bedfd7d253c3163f25584628115df1066a5e1b9ffa12f8a2cb
abl/base/prod.jsonl records(non-resources), resources, sha256(non-resources) = (2582, 1, '640b95d10248b2bedfd7d253c3163f25584628115df1066a5e1b9ffa12f8a2cb')
tests output identical except timing lines
"exitCode":
"exitCode":
"guardReason":
"exitCode":
"exitCode":
"guardReason":
r1_full_01 .exitCode 0
r1_full_01 .childOutcome.exitCode 0
r1_full_01 .childOutcome.guardReason None
r1_tests_01 .exitCode 0
r1_tests_01 .childOutcome.exitCode 0
r1_tests_01 .childOutcome.guardReason None
```

```
$ grep -nE 'sample_equal baseline rerun|resources' /var/projects/_review_scratch/s14a1-cmp-r1-grok/ablations/extract.txt | cut -c1-200
20:sample_equal baseline rerun True
28:baseline resources {"kind": "resources", "runtime_seconds": 105.5439142069954, "peak_rss_kib": 200432}
29:rerun resources {"kind": "resources", "runtime_seconds": 105.81404200100224, "peak_rss_kib": 199760}
```
