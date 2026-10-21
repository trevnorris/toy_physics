# Measurements: S14a part 1 comparator build review round 0 lookups (generated 2026-10-21 12:52)

Generator: `_scratch/s14a/gen/cmp_build_r0_lookups.sh` (grep/sed/sha256sum/diff/git/cat only; this file is
written by it). Mechanical lookups only; every CAS result quoted here is a leg's saved stdout.

## E1. The round: build, baseline, prompt, tree, reports

```
$ cat /var/projects/toy_physics/_scratch/s14a/cmp/build.DONE; grep -A1 -m1 'tokens used' /var/projects/toy_physics/_scratch/s14a/cmp/build.log
exit 0
tokens used
302,055
```

```
$ sed -n '2,11p' /var/projects/toy_physics/_scratch/s14a/cmp/build.log
--------
workdir: /var/projects/toy_physics
model: gpt-6.1-sol
provider: openai
approval: never
sandbox: danger-full-access
reasoning effort: xhigh
reasoning summaries: none
session id: 01a15f56-1181-75c2-8f4a-742595b33523
--------
```

```
$ cd /var/projects/toy_physics && sha256sum -c /var/projects/toy_physics/_scratch/s14a/cmp/review/s14a1_cmp_build_baseline_r0.sha256
research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py: OK
research/pde_ledger_v3/scripts/test_S14a_1_cross_engine_comparator.py: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_run_record.md: OK
research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_builder_report.md: OK
_scratch/s14a/cmp_work/full_02.jsonl: OK
_scratch/s14a/cmp_work/full_02_inventory.json: OK
_scratch/s14a/cmp_work/full_02_unpaired_components.json: OK
_scratch/s14a/cmp_work/tests_04.out: OK
```

```
$ cd /var/projects/toy_physics && git status --porcelain | diff - /var/projects/toy_physics/_scratch/s14a/cmp/review/status_before_cmp_r0.txt && echo 'status unchanged'
2d1
< ?? research/pde_ledger_v3/directives/_measurements/S14a_1_comparator_build_review_r0_lookups.md
```

```
$ sha256sum /var/projects/toy_physics/_scratch/s14a/cmp/review/s14a1_cmp_build_review_prompt_r0.md
f0958ce83ce27a47b806b1708ec1eaf901a097fa69dfe8e7c51d50904d5cf98a  /var/projects/toy_physics/_scratch/s14a/cmp/review/s14a1_cmp_build_review_prompt_r0.md
```

```
$ for r in tests_04 full_02; do printf '%s ' $r; grep -oE '"exitCode": *[0-9]+' /var/projects/toy_physics/_scratch/s11c/s14a1_cmp_guard/$r/outcome.json; done
tests_04 "exitCode": 0
"exitCode": 0
full_02 "exitCode": 0
"exitCode": 0
```

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r0-grok/DONE; grep -nE 'Verdict: ' /var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md /var/projects/_review_scratch/s14a1-cmp-r0-grok/REVIEW.md | head -2
exit 0
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:5:## Verdict: NEEDS REVISION
/var/projects/_review_scratch/s14a1-cmp-r0-grok/REVIEW.md:3:**Verdict: NEEDS REVISION**
```

```
$ grep -nE '^### ' /var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md /var/projects/_review_scratch/s14a1-cmp-r0-grok/REVIEW.md
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:9:### F1. Nested integrals: the reader builds an object neither engine printed (item 1)
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:34:### prod_c2a_storage_nested_slot_swap
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:36:### prod_c2b_order_lost_nested_slot_swap
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:41:### fixed_c2a_storage_nested_slot_swap vs fixed_baseline
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:44:### fixed_c2b_order_lost_nested_slot_swap vs fixed_baseline
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:58:### F2. A face derivative taken at another point cannot be seen by evaluation, and for three parts evaluation is the only resolved check (items 1, 5)
/var/projects/_review_scratch/s14a1-cmp-r0-claude/REVIEW.md:89:### F3 (claim scope, minor). Unjoined reasons are fixed strings, not computed facts (item 2)
/var/projects/_review_scratch/s14a1-cmp-r0-grok/REVIEW.md:7:### 1. The transverse-integration comparison note is unjoined under the coordinate-construction reason
```

## E2. Nested integrals (Claude F1)

```
$ sed -n '788p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
        return integral.xreplace({var:s.Symbol('_int'+str(i),real=True) for i,var in enumerate(integral.variables)})
```

```
$ grep -nF 'as read' /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/probe_fixture_capture_and_separable.stdout | head -6
3:   py_operand as read: Integral(_int0**2*Integral(_int0**2, (_int0, 0, 1)), (_int0, 1, 2))
4:   wl_operand as read: Integral(_int0**2*Integral(_int0**2, (_int0, 0, 1)), (_int0, 1, 2))
7:   py_operand as read: Integral(_int0**2*Integral(_int0**2, (_int0, 0, 1)), (_int0, 1, 2))
8:   wl_operand as read: Integral(_int0**2*Integral(_int0**2, (_int0, 0, 1)), (_int0, 1, 2))
```

```
$ grep -nE 'sympy.integrate of' /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/probe_fixture_capture_and_separable.stdout
10:   sympy.integrate of the py fixture object: 15/8
11:   sympy.integrate of the second wl fixture object: 7/9
```

```
$ grep -nE '^### |bindings moved' /var/projects/_review_scratch/s14a1-cmp-r0-claude/corrupt/diff_prod.out
1:### prod_c1_core_channel_form
3:   bindings moved: 1 of 227
6:### prod_c2a_storage_nested_slot_swap
7:   bindings moved: 0 of 227
10:### prod_c2b_order_lost_nested_slot_swap
11:   bindings moved: 0 of 227
14:### prod_c3_storage_face_derivative_point
16:   bindings moved: 1 of 227
```

```
$ grep -nE '^### |bindings moved' /var/projects/_review_scratch/s14a1-cmp-r0-claude/corrupt/diff_fixed.out
1:### fixed_baseline vs production full_02
11:   bindings moved: 9 of 227
14:### fixed_c2a_storage_nested_slot_swap vs fixed_baseline
16:   bindings moved: 1 of 227
19:### fixed_c2b_order_lost_nested_slot_swap vs fixed_baseline
21:   bindings moved: 1 of 227
```

```
$ diff /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/fixed/S14a_1_cross_engine_comparator.py
788c788
<         return integral.xreplace({var:s.Symbol('_int'+str(i),real=True) for i,var in enumerate(integral.variables)})
---
>         return integral.xreplace({var:s.Symbol('_int'+str(len(scope))+'_'+str(i),real=True) for i,var in enumerate(integral.variables)})
```

```
$ grep -nE 'share one symbol' /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/probe_nested_all.out | cut -c1-160
1:b169 | py | SphericalSpecialization/order lost | 1 | outer and inner bound variable share one symbol
2:b169 | wl | spherical specialization / shell/order lost | 1 | outer and inner bound variable share one symbol
3:b171 | py | SphericalSpecialization/storage | 1 | outer and inner bound variable share one symbol
4:b171 | wl | spherical specialization / shell/storage | 1 | outer and inner bound variable share one symbol
5:b173 | py | SphericalSpecialization/shell balance/operands | 1 | outer and inner bound variable share one symbol
6:b173 | wl | spherical specialization / shell/balance | 1 | outer and inner bound variable share one symbol
7:b174 | py | SphericalSpecialization/shell balance/residual | 1 | outer and inner bound variable share one symbol
8:b174 | wl | spherical specialization / shell/balance residual | 1 | outer and inner bound variable share one symbol
9:b175 | py | SphericalSpecialization/face channel minus order lost/residual | 1 | outer and inner bound variable share one symbol
10:b175 | wl | spherical specialization / shell/face channel minus order lost | 1 | outer and inner bound variable share one symbol
```

```
$ grep -nE 'point 0:' /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/probe_capture_values.stdout | cut -c1-200
2:   point 0: py=-50197.6136518578 wl=-50197.6136518578 py-wl=0 | full_02 printed py=-50197.6136518578
6:   point 0: py=558.971447409312 wl=558.971447409312 py-wl=0 | full_02 printed py=558.971447409312
10:   point 0: py=161630.296988830 wl=161630.296988830 py-wl=0 | full_02 printed py=161630.296988830
14:   point 0: py=105317.232252639 wl=105317.232252639 py-wl=0 | full_02 printed py=105317.232252639
18:   point 0: py=-75997.2311562359 wl=-75997.2311562359 py-wl=0 | full_02 printed py=-50197.6136518578
22:   point 0: py=558.971447409312 wl=558.971447409312 py-wl=0 | full_02 printed py=558.971447409312
26:   point 0: py=161630.296988830 wl=161630.296988830 py-wl=0 | full_02 printed py=161630.296988830
30:   point 0: py=131116.849757017 wl=131116.849757017 py-wl=0 | full_02 printed py=105317.232252639
```

## E3. Separable test functions and the face-derivative point (Claude F2)

```
$ sed -n '837,844p;854,855p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
def basic_form(form,args,code):
    if form=='affine':
        return (code+1)/7+sum((i+1)*(code+1)*q/101 for i,q in enumerate(args))
    if form=='quadratic':
        return (code+1)/7+sum(((i+1)**2*(code+1)*q*q+(i+1)*q)/103 for i,q in enumerate(args))
    if form=='cubic':
        return (code+1)/7+sum(((i+1)**3*(code+1)*q**3+(i+1)**2*q)/107 for i,q in enumerate(args))
    raise InputError('unknown probe form '+form)
    lower=-1+basic_form(form,args,s.Integer(1))/100
    upper=lower+3+(1+sum((i+1)*q*q for i,q in enumerate(args)))/100
```

```
$ grep -nE "\('face'" /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/probe_fixture_capture_and_separable.stdout | head -20
31:   ('face', 'affine', 'density_zero', 'all mixed partials zero') 4
32:   ('face', 'affine', 'generic', 'all mixed partials zero') 4
33:   ('face', 'affine', 'jump_lower', 'all mixed partials zero') 4
34:   ('face', 'affine', 'jump_upper', 'all mixed partials zero') 4
35:   ('face', 'affine', 'order_zero_lower', 'all mixed partials zero') 4
36:   ('face', 'affine', 'order_zero_upper', 'all mixed partials zero') 4
37:   ('face', 'cubic', 'density_zero', 'all mixed partials zero') 4
38:   ('face', 'cubic', 'generic', 'all mixed partials zero') 4
39:   ('face', 'cubic', 'jump_lower', 'all mixed partials zero') 4
40:   ('face', 'cubic', 'jump_upper', 'all mixed partials zero') 4
41:   ('face', 'cubic', 'order_zero_lower', 'all mixed partials zero') 4
42:   ('face', 'cubic', 'order_zero_upper', 'all mixed partials zero') 4
43:   ('face', 'quadratic', 'density_zero', 'all mixed partials zero') 4
44:   ('face', 'quadratic', 'generic', 'all mixed partials zero') 4
45:   ('face', 'quadratic', 'jump_lower', 'all mixed partials zero') 4
46:   ('face', 'quadratic', 'jump_upper', 'all mixed partials zero') 4
47:   ('face', 'quadratic', 'order_zero_lower', 'all mixed partials zero') 4
48:   ('face', 'quadratic', 'order_zero_upper', 'all mixed partials zero') 4
```

```
$ grep -n '^==' /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/exact_nonzero_parts.out
1:== b50 part=item 1: closed integral of material differential operand
4:== b165 part=items 2-3: shell domain
7:== b173 part=shell left
10:== b174 part=item 2: shell balance residual
```

```
$ grep -nE 'c3|b171' /var/projects/_review_scratch/s14a1-cmp-r0-claude/corrupt/diff_prod.out | cut -c1-260
14:### prod_c3_storage_face_derivative_point
15:   b171 SphericalSpecialization/storage: moved [exact residual; also: on_supplied_domain_residual] baseline {'samples': 54, 'nonzero': 0, 'unresolved': 0, 'exact': 1, 'exact_nonzero': 1, 'not_formed': 0} -> variant {'samples': 54, 'nonzero': 0, 'unresolved'
```

```
$ grep -nE '^==|distinct sample residuals' /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/b171_exact_residuals.out | cut -c1-140
1:== production full_02: b171 component keys: ['as_computed_residual', 'on_supplied_domain_residual', 'part', 'py_operand', 'samples', 'wl_o
4:   distinct sample residuals: ['0'] n= 54
5:== production on c3: b171 component keys: ['as_computed_residual', 'on_supplied_domain_residual', 'part', 'py_operand', 'samples', 'wl_ope
8:   distinct sample residuals: ['0'] n= 54
```

```
$ sed -n 2p /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/b171_exact_residuals.out
   as_computed_residual: 4*pi*Integral(-_int0**2*SphericalSpecialization_inside_n(_int0, SphericalSpecialization_face_lower(_int0, t), t)*Derivative(SphericalSpecialization_face_lower(_int0, t), t), (_int0, inner_radius, outer_radius)) + 4*pi*Integral(_int0**2*SphericalSpecialization_inside_n(_int0, SphericalSpecialization_face_lower(_int0, t), t)*Derivative(SphericalSpecialization_face_lower(_int0, t), t), (_int0, inner_radius, outer_radius))
```

```
$ sed -n 7,12p /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/exact_nonzero_parts.out | cut -c1-400
== b173 part=shell left
   as_computed_residual: 4*pi*Integral(-_int0**2*SphericalSpecialization_inside_n(_int0, SphericalSpecialization_face_lower(_int0, t), t)*Derivative(SphericalSpecialization_face_lower(_int0, t), t), (_int0, inner_radius, outer_radius)) + 4*pi*Integral(_int0**2*SphericalSpecialization_inside_n(_int0, SphericalSpecialization_face_lower(_int0, t), t)*Derivative(SphericalSpecialization_face_lower(_int
   distinct sample residuals: ['0']
== b174 part=item 2: shell balance residual
   as_computed_residual: 4*pi*Integral(-_int0**2*SphericalSpecialization_inside_n(_int0, SphericalSpecialization_face_lower(_int0, t), t)*Derivative(SphericalSpecialization_face_lower(_int0, t), t), (_int0, inner_radius, outer_radius)) + 4*pi*Integral(_int0**2*SphericalSpecialization_inside_n(_int0, SphericalSpecialization_face_lower(_int0, t), t)*Derivative(SphericalSpecialization_face_lower(_int
   distinct sample residuals: ['0']
```

## E4. Unjoined reasons (Claude F3, Grok 1)

```
$ sed -n '1305,1339p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
def unjoined_reason(tag,engine):
    if tag.endswith('/comparison'):
        return 'Engine-local comparison note; no unused note about the same two counterparts, or split/grouped note has no injective note counterpart (py27-36,201-205; wl14-15,146,162).'
    if '/projection integration/' in tag or '/closed integration/' in tag:
        return 'Construction residuals subtract the two integration routes in opposite orders (py201-205; wl161-162,15); no residual or note binding.'
    if '/coordinate ' in tag or '/pulled current' in tag or '/transverse ' in tag or tag.startswith('transverse integration'):
        return 'Wolfram-only fixed-slab coordinate construction; SymPy prints native-w integration objects (py171-185; wl64-99,108-125,143-155). Different printed objects; no coordinate reconciliation.'
    if '/bulk derivative integral' in tag:
        return 'Native-w derivative integral or its endpoints (py171-176); Wolfram prints the pulled z-current derivative and relative endpoints (wl143,151-155). Different printed integrands/current objects.'
    if '/supplied material factor' in tag or '/supplied material face factor' in tag:
        return 'SymPy prints the one-sided trace w=face (py213-216); Wolfram prints unrestricted w-dependent factors (wl48,118). Different objects; faces are not substituted to reconcile them.'
    if 'dimension' in tag or '/sum ' in tag:
        return 'Engine-specific dimension audit occurrence: SymPy operand/sum/reference expression trees (py305-313,316-337), or Wolfram assignments/atoms/derivative/source/term ledgers and extra shell objects (wl244-256). No unused source-cited counterpart; term identity and position do not create bindings.'
    if tag.endswith('/projected velocity') and engine=='py':
        return 'SymPy prints velocity at free radius (py271,275); Wolfram prints it only at outer radius ro (wl218-219). The outer-radius conditional counterpart is bound separately.'
    if tag.endswith('/regular parts'):
        return 'Regional differential objects are already compared through supplied equation lhs components (py120-121,163,167; wl130,136); this redundant standalone occurrence has no unused counterpart.'
    if '/regional balance residuals' in tag or '/source law residual' in tag:
        return 'Wolfram construction checks subtract native supplied rhs (wl137,178); no SymPy check against that same counterpart is printed.'
    if '/encoded ' in tag:
        return 'Only SymPy prints the encoded global current before differentiation (py157-158,164-165); Wolfram prints regional currents and their differentiated distribution (wl129-138).'
    if tag.endswith('/projected differential expansion'):
        return 'Only SymPy prints the standalone projected differential expansion (py181,192); Wolfram prints it combined with face terms (wl158-170).'
    if tag.endswith('/face channel by face'):
        return 'Only SymPy prints shell-integrated transport by individual face (py251,263); Wolfram shell objects group interior/concentrated contributions across both faces (wl201-203,222).'
    if 'face channel minus order lost/operands' in tag:
        return 'Only SymPy prints the two operands separately (py270,33-34); Wolfram prints their difference (wl215,222), which is bound to the SymPy residual.'
    if tag.endswith('/projected source'):
        return 'Material source literal zero is already compared in the projected balance rhs (py188,198-200; wl163,173); this redundant occurrence has no unused counterpart.'
    if '/supplied/' in tag or '/supplied material face equation/' in tag or tag.endswith('/velocity'):
        return 'SymPy redundant supplied-field/equation occurrence (py109-148,197); corresponding current and regional equation components are bound where both are printed (wl129-136,172). No unused counterpart for this occurrence.'
    if engine=='wl':
        return 'Wolfram-only intermediate or redundant occurrence (wl114-123,134-180,199-222); no unused SymPy occurrence of this object. Named projected/shell objects are bound separately.'
    return 'SymPy-only geometry, premise, engine metadata or redundant construction occurrence (py220-238,127-217,249-277); no unused Wolfram occurrence of this object.'

```

```
$ sed -n '494,495p' /var/projects/toy_physics/research/pde_ledger_v3/scripts/S14a_1_cross_engine_comparator.py
def is_label(tag):
    return tag in ('contract','engine','WORKING_ON') or tag.endswith(('/comparison','/representation'))
```

```
$ head -40 /var/projects/_review_scratch/s14a1-cmp-r0-claude/probe/v08b_unjoined_reasons.out | cut -c1-260
46 ('py', 'SymPy redundant supplied-field/equation occurrence (py109-148,197); corresponding current and regional equation components are bound where both are printed (wl129-136,172). No unused counterpart for t') | e.g. tag General/supplied/below/material_den
16 ('py', 'SymPy-only geometry, premise, engine metadata or redundant construction occurrence (py220-238,127-217,249-277); no unused Wolfram occurrence of this object.') | e.g. tag General/material/face 0/interior leaving term
4 ('py', 'Construction residuals subtract the two integration routes in opposite orders (py201-205; wl161-162,15); no residual or note binding.') | e.g. tag General/material/projection integration/operands
19 ('py', 'Engine-specific dimension audit occurrence: SymPy operand/sum/reference expression trees (py305-313,316-337), or Wolfram assignments/atoms/derivative/source/term ledgers and extra shell objects (wl244') | e.g. tag General/dimensions/material/flux 0/
20 ('wl', 'Wolfram-only intermediate or redundant occurrence (wl114-123,134-180,199-222); no unused SymPy occurrence of this object. Named projected/shell objects are bound separately.') | e.g. tag general projection/supplied gamma identity
8 ('wl', 'Construction residuals subtract the two integration routes in opposite orders (py201-205; wl161-162,15); no residual or note binding.') | e.g. tag general projection/material/closed integration/operand_a
10 ('wl', 'Engine-specific dimension audit occurrence: SymPy operand/sum/reference expression trees (py305-313,316-337), or Wolfram assignments/atoms/derivative/source/term ledgers and extra shell objects (wl244') | e.g. tag general projection/dimensions/mater
```

```
$ cat /var/projects/_review_scratch/s14a1-cmp-r0-grok/check_endswith.stdout | cut -c1-300
SNIPPET_POINTS
['0x64', '0x65', '0x66', '0x20', '0x75', '0x6e', '0x6a', '0x6f', '0x69', '0x6e', '0x65', '0x64', '0x5f', '0x72', '0x65', '0x61', '0x73', '0x6f', '0x6e', '0x28', '0x74', '0x61', '0x67', '0x2c', '0x65', '0x6e', '0x67', '0x69', '0x6e', '0x65', '0x29', '0x3a', '0xa', '0x20', '0x20', '0x20', '0x20', '0x6
TAG_TAIL ['0x2f', '0x20', '0x63', '0x6f', '0x6d', '0x70', '0x61', '0x72', '0x69', '0x73', '0x6f', '0x6e']
SUFFIX ['0x2f', '0x63', '0x6f', '0x6d', '0x70', '0x61', '0x72', '0x69', '0x73', '0x6f', '0x6e']
ENDSWITH False
STARTSWITH True
FN Wolfram-only fixed-slab coordinate construction; SymPy prints native-w
```

```
$ grep -nF 'transverse integration / comparison' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | cut -c1-120
9:transverse integration / comparison: "There is nothing to compare this object with."
```

```
$ grep -cE ' / comparison: ' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out
1
```

## E5. Shell domain printed by each engine (Claude note; Grok checked list)

```
$ grep -nE '^SphericalSpecialization/shell domain: ' /var/projects/toy_physics/research/pde_ledger_v3/scripts/out/S14a_1_projection.out
336:SphericalSpecialization/shell domain: inner_radius < outer_radius
```

```
$ grep -nE '^spherical specialization / shell/domain: ' /var/projects/toy_physics/research/pde_ledger_v3/mathematica/out/S14a_1_projection.out | cut -c1-160
237:spherical specialization / shell/domain: Inactive[Less][0, ri, ro]
```

## E6. Reproducibility, both legs

```
$ grep -nE 'identical excluding|sha256 full_02 minus' /var/projects/_review_scratch/s14a1-cmp-r0-claude/repro/repro_compare.out
2:identical excluding the resources record: True
3:sha256 full_02 minus resources : 5cdef2382b84e2baf951c284ef118d82f95d95c65dcf83545b05ea593276f2ae
```

```
$ grep -nE 'mismatch|lines' /var/projects/_review_scratch/s14a1-cmp-r0-grok/diff_repro.stdout | head -5
1:lines 2583 2583
4:non_resource_mismatches 0
5:mismatch_kinds {}
```
