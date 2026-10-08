#!/usr/bin/env python3
"""Run isolated comparator-logic ablations, under the caller's pooled guard.

No production transcripts or engine imports. Every variant runs its designated
serialized synthetic control(s). Assertion failures are required; operational errors,
kills, and surviving variants stop the run. No retries or duration limits.
"""
import argparse
import difflib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


# One entry for each amendment-1 item-6 bullet, in directive order.
# The two further tests and the O4 grain repair have supplemental ablations.
# Each edit removes comparator logic; tests and their assertions are copied intact.
def A(bullet,label,required,*edits):
    return dict(bullet=bullet,label=label,required=(required,) if isinstance(required,str) else required,edits=edits)


def suppress(field):
    return ('field:bag_delta(left.get(field,{}),right.get(field,{}))',
            f'field:[] if field=={field!r} else bag_delta(left.get(field,{{}}),right.get(field,{{}}))')


_OPERANDS = "        emit(output,{'kind':'operands',**info,'py':data(a) if a else None,'wl':data(b) if b else None})\n"
_RESIDUAL = "        emit(output,{'kind':'residual','row':row.label,'comparison':result})\n"
ABLATIONS = (
    A('R1','subtraction_to_addition','test_residual_reconstructs_left_operand',
      ('delta = sp.cancel(left-right)','delta = sp.cancel(left+right)')),
    A('R2','operands_after_residual','test_operands_precede_residual',
      (_OPERANDS,''),(_RESIDUAL,_RESIDUAL+_OPERANDS)),
    A('R3','constant_residual','test_one_sided_operand',
      ('delta = sp.cancel(left-right)','delta = sp.Integer(0)')),
    A('R4','form_ignored','test_form_change',
      ('delta = sp.cancel(left-right)','delta = sp.Integer(len(left.free_symbols)-len(right.free_symbols))')),
    A('J1','repoints_ignored',('test_each_name_binding_repoint','test_each_join_binding_repoint'),
      ("return self.maps[engine].get(value, engine + '::' + value)",
       "return next((r.py for r in NAME_TABLE if getattr(r,engine)==value),engine+'::'+value)"),
      ('a,b = extract(py,row.py),extract(wl,row.wl)',
       'a,b = extract(py,row.py),extract(wl,next((r.wl for r in JOIN_TABLE if r.label==row.label),row.wl))')),
    A('J2','join_requires_same_tag','test_moved_tag_and_key_stays_joined',
      ('a,b = extract(py,row.py),extract(wl,row.wl)',
       'a,b = (extract(py,row.py),extract(wl,row.wl)) if row.py[0]==row.wl[0] else (None,None)')),
    A('J3','unsupported_as_zero','test_computed_and_open_same_role',
      ("except Unsupported as exc:\n        return {'outcome':'not_formed','reason':str(exc)},(0,0)",
       "except Unsupported as exc:\n        return {'outcome':'exact','residual':'0'},(0,0)")),
    A('J4','named_inventory_declared_only','test_unbound_open_name_is_in_difference',
      ("if v.head=='Text' or bindings.kinds.get(key) not in ('coordinate','parameter','profile','binder'):",
       "if key in bindings.kinds and (v.head=='Text' or bindings.kinds.get(key) not in ('coordinate','parameter','profile','binder')):")),
    A('P1','stripped_argument_restored','test_applied_argument_stripped',
      ("if bindings.kinds.get(key) not in ('coordinate','parameter'):",
       "if bindings.kinds.get(key)=='profile':\n            return sp.Function(key)(sp.Symbol('x1',real=True))\n        if bindings.kinds.get(key) not in ('coordinate','parameter'):")),
    A('P2','live_profile_frozen','test_live_profile_frozen',
      ('return sp.Function(key)(*(algebra(x,bindings,engine) for x in a[1:]))','return sp.Integer(11)')),
    A('P3','derivative_order_ignored','test_derivative_order',
      ('(var,int(orders[0].value))','(var,1)')),
    A('P4','evaluation_point_collapsed','test_derivative_evaluation_points',
      ('    return algebra(n,bindings,engine)\n\n\ndef layout',
       "    return sp.Symbol('fixed_point',real=True)\n\n\ndef layout")),
    A('P5','binder_scope_erased','test_contract_binder_structure',
      ("return scope.get(n.value,atom('Name',bindings.name(n.value,engine)))",
       "return atom('Name',bindings.name(n.value,engine))")),
    A('O1','head_difference_omitted','test_contract_changed_open_head',suppress('head')),
    A('O2','named_difference_omitted','test_contract_changed_named_operand',suppress('named_OPEN_operands')),
    A('O3','labels_omitted','test_contract_changed_label',
      ("if v.head in ('Symbol','FunctionName','Dummy','Text'):",
       "if v.head in ('Symbol','FunctionName','Dummy'):")),
    A('O4','live_difference_omitted','test_contract_changed_live_object',suppress('live_arguments')),
    A('O5','orientation_difference_omitted','test_contract_changed_orientation',suppress('orientation')),
    A('O6','canonical_key_removed','test_identical_open_objects_have_empty_deltas',
      ("return Node('ProfileDerivative',(seq(*(convert(v) for v in orders)),",
       "return Node('WolframDerivativeSyntax',(seq(*(convert(v) for v in orders)),")),
    A('O7','missing_sibling_as_exact','test_nested_sibling_removed',
      ("result[key] = {'outcome':'not_formed','reason':'nested sibling absent',",
       "result[key] = {'outcome':'exact','reason':'nested sibling absent',")),
    A('B1','plain_sign','test_plain_sum_orientation',('sign *= -1','sign *= 1')),
    A('B2','held_aggregate_sign','test_held_aggregate_orientation',
      ('            role = None\n','            linear_arguments = ()\n            role = None\n')),
    A('B3','sympy_body_sign','test_all_held_linear_orientations',
      ('signed_actions(child,sign if i == body else 1,dependency or i != body)',
       'signed_actions(child,1,dependency or i != body)')),
    A('B4','closed_part_first_term','test_each_closed_term_in_multiterm_balance',
      ("Node('Add',tuple(closed_terms)) if closed_terms","Node('Add',tuple(closed_terms[:1])) if closed_terms")),
    A('L1','transpose_dropped','test_transpose_layout_matches_explicit_transpose',
      ("if row.layout == 'transpose' and a.head == 'Matrix':",'if False:')),
    A('L2','relation_left_only','test_each_relation_operand_both_engines',
      ('for x,y in zip(a.args,b.args):','for x,y in zip(a.args[:1],b.args[:1]):')),
    A('L3','sqrt_to_cbrt','test_wolfram_sqrt_translation',
      ('return sp.sqrt(algebra(a[1], bindings, engine))','return sp.cbrt(algebra(a[1], bindings, engine))')),
    A('I1','repeated_live_omitted','test_limit_repeated_live_object',
      ("'repeated_live_objects':{k:v for k,v in live_counts.items() if v>1}","'repeated_live_objects':{}")),
    A('I2','outside_count_omitted','test_limit_outside_leaves',
      ("'outside_inventory_leaves':parsed-consumed","'outside_inventory_leaves':0")),
    A('T1','verdict_token_inserted','test_output_has_no_verdict_tokens',
      ("emit(output,{'kind':'comparison_scope','OPEN_limit':OPEN_SCOPE})",
       "emit(output,{'kind':'comparison_scope','OPEN_limit':OPEN_SCOPE,'VERDICT':'inserted'})")),
    A('extra-duplicate','duplicate_tables_admitted',('test_duplicate_join_rejected','test_duplicate_name_rejected_both_sides'),
      ('def validate_tables(joins, names):','def validate_tables(joins, names):\n    return')),
    A('extra-boolean','boolean_as_algebra','test_boolean_does_not_hide_algebraic_sibling',
      ("if a.head in ('Text','Boolean','Name','FunctionName') or b.head in ('Text','Boolean','Name','FunctionName'):",
       "if a.head in ('Text','Name','FunctionName') or b.head in ('Text','Name','FunctionName'):"),
      ('    h, a = n.head, n.args',"    h, a = n.head, n.args\n    if h=='Boolean':\n        return sp.Integer(n.value.lower()=='true')")),
    A('extra-O4','o4_container_grain','test_o4_components_have_same_role_joins',
      ("J('coupled_embedding_operand','COUPLED_INPUTS/embedding/0','B_HOLD_LIVE/O4Identity/@1',395,239,'§3.2 O4 operand'),\n    J('coupled_embedding_count','COUPLED_INPUTS/embedding/2','COUPLED_INPUTS_MODEL_POINT/O4EquationIdentityCount',396,340,'§3.2 unsettled equation identity/count'),",
       "J('coupled_embedding','COUPLED_INPUTS/embedding','B_HOLD_LIVE/O4Identity',395,239,'§3.2'),")),
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not os.environ.get('S11C_POOLED_GUARD_MANIFEST'):
        raise RuntimeError('run through s11c_guarded_run.py in pooled mode')
    args.output.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).with_name('O2_cross_engine_comparator.py').read_text()
    tests=Path(__file__).with_name('test_O2_cross_engine_comparator.py').read_text()
    summaries=[]
    for entry in ABLATIONS:
        label,required=entry['label'],entry['required']
        modified=source
        for before,after in entry['edits']:
            if modified.count(before)!=1:
                raise RuntimeError('ablation anchor is not unique: '+label)
            modified=modified.replace(before,after)
        folder=args.output/label
        folder.mkdir()
        (folder/'O2_cross_engine_comparator.py').write_text(modified)
        (folder/'test_O2_cross_engine_comparator.py').write_text(tests)
        (folder/'change.diff').write_text(''.join(difflib.unified_diff(
            source.splitlines(True),modified.splitlines(True),fromfile='baseline',tofile=label)))
        command=['/usr/bin/time','-f','elapsed_seconds=%e\npeak_rss_kib=%M',
                 '-o',str((folder/'resources.txt').resolve()),sys.executable,
                 str((folder/'test_O2_cross_engine_comparator.py').resolve()),
                 *('Controls.'+name for name in required)]
        with (folder/'stdout').open('w') as out,(folder/'stderr').open('w') as err:
            result=subprocess.run(command,stdout=out,stderr=err,
                                  env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
        transcript=(folder/'stderr').read_text()
        failures=re.findall(r'^FAIL: (.+)$',transcript,re.M)
        errors=re.findall(r'^ERROR: (.+)$',transcript,re.M)
        summary={'bullet':entry['bullet'],'ablation':label,'command':command,'exit_code':result.returncode,
                 'failing_tests':failures,'errors':errors,'required_tests':required,
                 'resources':(folder/'resources.txt').read_text()}
        summaries.append(summary)
        (args.output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
        print(json.dumps(summary),flush=True)
        if result.returncode!=1 or errors or not all(any(name in test for test in failures) for name in required):
            raise RuntimeError('ablation did not produce its required assertion failure: '+label)


if __name__=='__main__':
    main()
